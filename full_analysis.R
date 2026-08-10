# ==============================================================================
# BTC1859 Group Project - Sleep Disturbance and Quality of Life
# ==============================================================================

library(dplyr)
library(tidyr)
library(janitor)
library(gtsummary)
library(gt)
library(binom) # for Wilson 95% CIs
library(car)
library(ggplot2)
library(patchwork)
library(broom)
library(paletteer)
library(mice)

# Load data --------------------------------------------------------------------

df <- read.csv(
  "project_data.csv",
  na.strings = "NA"
)

# Explore data -----------------------------------------------------------------

str(df)
dim(df)
sum(is.na(df)) # 2191 NAs
colSums(is.na(df)) 

# Subset data ------------------------------------------------------------------

# Keeping important demographic and clinical variables 

df_updated <- df |>
  dplyr::select(
    Subject, Age, Gender, BMI, Time.from.transplant,
    Liver.Diagnosis, Recurrence.of.disease, Rejection.graft.dysfunction,
    Any.fibrosis, Renal.Failure, Depression, Corticoid,
    Pittsburgh.Sleep.Quality.Index.Score, Epworth.Sleepiness.Scale,
    Berlin.Sleepiness.Scale, Athens.Insomnia.Scale,
    SF36.PCS, SF36.MCS
  )

# Clean data -------------------------------------------------------------------

# Create copy of dataframe 
df_cleaned <- df_updated

# Convert column names to snake case (in line with Tidyverse)
df_cleaned <- clean_names(df_cleaned)

# Rename column names for conciseness
df_cleaned <- df_cleaned |> rename(
  psqi_score = pittsburgh_sleep_quality_index_score,
  ess_score = epworth_sleepiness_scale,
  bss_score = berlin_sleepiness_scale,
  ais_score = athens_insomnia_scale
)

# Convert select columns to factors (from int)
df_cleaned <- df_cleaned |>
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2)
    ),
    
    
    liver_diagnosis = factor(
      liver_diagnosis,
      levels = c(1, 2, 3, 4, 5)
    ),
    
    across(
      c(
        recurrence_of_disease,
        rejection_graft_dysfunction,
        any_fibrosis,
        renal_failure,
        depression,
        corticoid
      ),
      ~ factor(
        .x,
        levels = c(0, 1)
      )
    )
  )

# Verify that age entries are valid
min(df_cleaned$age, na.rm = TRUE)
max(df_cleaned$age, na.rm = TRUE)

# Verify that BMI entries are valid
min(df_cleaned$bmi, na.rm = TRUE)
max(df_cleaned$bmi, na.rm = TRUE)

# Verify that liver diagnosis entries are valid (i.e. 1-5)
unique(df_cleaned$liver_diagnosis)

# Verify that sleep disturbance measure scores are valid (based on their coding
# scheme)
invalid_sleep_scores <- df_cleaned |>
  select(
    subject,
    psqi_score,
    ess_score,
    ais_score
  ) |>
  pivot_longer(
    cols = c(psqi_score, ess_score, ais_score),
    names_to = "sleep_measure",
    values_to = "score"
  ) |>
  filter(
    case_when(
      sleep_measure == "psqi_score" ~ !is.na(score) & !between(score, 0, 21),
      sleep_measure == "ess_score" ~ !is.na(score) & !between(score, 0, 24),
      sleep_measure == "ais_score" ~ !is.na(score) & !between(score, 0, 24)
    )
  )

# Correct the invalid ESS score by marking it as NA
df_cleaned <- df_cleaned |>
  mutate(
    ess_score = if_else(
      ess_score > 24,
      NA,
      ess_score
    )
  )

# Create categorical variables for sleep disturbance instruments ---------------

df_cleaned <- df_cleaned |>
  mutate(
    psqi_binary = case_when(
      is.na(psqi_score) ~ NA_integer_,
      psqi_score > 4    ~ 1L,
      TRUE        ~ 0L
    ),
    
    ess_binary = case_when(
      is.na(ess_score) ~ NA_integer_,
      ess_score > 10   ~ 1L,
      TRUE       ~ 0L
    ),
    
    ais_binary = case_when(
      is.na(ais_score) ~ NA_integer_,
      ais_score > 5    ~ 1L,
      TRUE       ~ 0L
    )
  ) |>
  relocate(psqi_binary, .after = psqi_score) |>
  relocate(ess_binary, .after = ess_score) |>
  relocate(ais_binary, .after = ais_score)

# Define reusable column vectors -----------------------------------------------

sleep_measure_scores <- c(
  "psqi_score",
  "ess_score",
  "ais_score"
)

sleep_measure_binary <- c(
  "psqi_binary",
  "ess_binary",
  "ais_binary",
  "bss_score"
)

clinical_variables <- c(
  "liver_diagnosis",
  "recurrence_of_disease",
  "rejection_graft_dysfunction",
  "any_fibrosis",
  "renal_failure",
  "depression",
  "corticoid",
  "time_from_transplant"
)

binary_clinical_variables <- c(
  "recurrence_of_disease",
  "rejection_graft_dysfunction",
  "any_fibrosis",
  "renal_failure",
  "depression",
  "corticoid"
)

demographic_variables <- c(
  "age",
  "gender",
  "bmi"
)

# Create descriptive summary table ---------------------------------------------

# Create descriptive statistics table (overall, ungrouped)
descriptive_table <- df_cleaned |>
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    liver_diagnosis = factor(
      liver_diagnosis,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "Hepatitis C",
        "Hepatitis B",
        "PSC/PBC/AIH",
        "Alcohol-related",
        "Other"
      )
    ),
    
    bss_score = factor(
      bss_score,
      levels = c(0, 1),
      labels = c(
        "Low likelihood",
        "High likelihood"
      )
    ),
    
    across(all_of(
      binary_clinical_variables
    ),
    ~ factor(
      .x,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    )
  ) |>
  tbl_summary(
    include = all_of(c(
      demographic_variables,
      clinical_variables,
      sleep_measure_scores,
      "bss_score"
    )
    ),
    
    type = bss_score ~ "dichotomous",
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    
    value = list(
      all_of(binary_clinical_variables) ~ "Yes",
      bss_score ~ "High likelihood"
    ),
    
    label = list(
      age = "Age",
      gender = "Gender",
      bmi = "Body Mass Index (BMI)",
      liver_diagnosis = "Liver diagnosis",
      recurrence_of_disease = "Recurrence of disease",
      rejection_graft_dysfunction = "Rejection or graft dysfunction",
      any_fibrosis = "Any fibrosis (grade A2 and higher)",
      renal_failure = "Renal failure",
      depression = "Depression",
      corticoid = "Corticosteroid use",
      time_from_transplant = "Time from transplant (years)",
      psqi_score = "Pittsburgh Sleep Quality Index score",
      ess_score = "Epworth Sleepiness Scale score",
      ais_score = "Athens Insomnia Scale score",
      bss_score = "High likelihood of sleep-disordered breathing (Berlin Sleepiness Scale)"
    ),
    
    missing = "no"
  )

descriptive_table |>
  as_gt() |>
  gtsave("tbl_characteristics.png")

# Create descriptive statistics table stratified by PSQI status ----------------

psqi_table <- df_cleaned |>
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    liver_diagnosis = factor(
      liver_diagnosis,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "Hepatitis C",
        "Hepatitis B",
        "PSC/PBC/AIH",
        "Alcohol-related",
        "Other"
      )
    ),
    
    psqi_binary = factor(
      psqi_binary,
      levels = c(0, 1),
      labels = c(
        "No sleep disturbance",
        "Yes sleep disturbance"
      )
    ),
    
    across(all_of(
      binary_clinical_variables
    ),
    ~ factor(
      .x,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    )
  ) |>
  tbl_summary(
    by = psqi_binary,
    
    include = all_of(c(
      demographic_variables,
      clinical_variables
    )
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    
    value = list(
      all_of(binary_clinical_variables) ~ "Yes"
    ),
    
    label = list(
      age = "Age",
      gender = "Gender",
      bmi = "Body Mass Index (BMI)",
      liver_diagnosis = "Liver diagnosis",
      recurrence_of_disease = "Recurrence of disease",
      rejection_graft_dysfunction = "Rejection or graft dysfunction",
      any_fibrosis = "Any fibrosis (grade A2 and higher)",
      renal_failure = "Renal failure",
      depression = "Depression",
      corticoid = "Corticosteroid use",
      time_from_transplant = "Time from transplant (years)"
    ),
    
    missing = "no"
  )

# Create descriptive statistics table stratified by ESS status -----------------
# (NOTE: this used to be a copy-paste of the PSQI table above, stratifying by
# psqi_binary and overwriting the psqi_table object. Fixed below to actually
# stratify by ess_binary and save to its own object, ess_table.)

ess_table <- df_cleaned |>
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    liver_diagnosis = factor(
      liver_diagnosis,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "Hepatitis C",
        "Hepatitis B",
        "PSC/PBC/AIH",
        "Alcohol-related",
        "Other"
      )
    ),
    
    ess_binary = factor(
      ess_binary,
      levels = c(0, 1),
      labels = c(
        "No sleep disturbance",
        "Yes sleep disturbance"
      )
    ),
    
    across(all_of(
      binary_clinical_variables
    ),
    ~ factor(
      .x,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    )
  ) |>
  tbl_summary(
    by = ess_binary,
    
    include = all_of(c(
      demographic_variables,
      clinical_variables
    )
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    
    value = list(
      all_of(binary_clinical_variables) ~ "Yes"
    ),
    
    label = list(
      age = "Age",
      gender = "Gender",
      bmi = "Body Mass Index (BMI)",
      liver_diagnosis = "Liver diagnosis",
      recurrence_of_disease = "Recurrence of disease",
      rejection_graft_dysfunction = "Rejection or graft dysfunction",
      any_fibrosis = "Any fibrosis (grade A2 and higher)",
      renal_failure = "Renal failure",
      depression = "Depression",
      corticoid = "Corticosteroid use",
      time_from_transplant = "Time from transplant (years)"
    ),
    
    missing = "no"
  )

# Create descriptive statistics table stratified by AIS status -----------------

ais_table <- df_cleaned |>
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    liver_diagnosis = factor(
      liver_diagnosis,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "Hepatitis C",
        "Hepatitis B",
        "PSC/PBC/AIH",
        "Alcohol-related",
        "Other"
      )
    ),
    
    ais_binary = factor(
      ais_binary,
      levels = c(0, 1),
      labels = c(
        "No sleep disturbance",
        "Yes sleep disturbance"
      )
    ),
    
    across(all_of(
      binary_clinical_variables
    ),
    ~ factor(
      .x,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    )
  ) |>
  tbl_summary(
    by = ais_binary,
    
    include = all_of(c(
      demographic_variables,
      clinical_variables
    )
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    
    value = list(
      all_of(binary_clinical_variables) ~ "Yes"
    ),
    
    label = list(
      age = "Age",
      gender = "Gender",
      bmi = "Body Mass Index (BMI)",
      liver_diagnosis = "Liver diagnosis",
      recurrence_of_disease = "Recurrence of disease",
      rejection_graft_dysfunction = "Rejection or graft dysfunction",
      any_fibrosis = "Any fibrosis (grade A2 and higher)",
      renal_failure = "Renal failure",
      depression = "Depression",
      corticoid = "Corticosteroid use",
      time_from_transplant = "Time from transplant (years)"
    ),
    
    missing = "no"
  )

# Create descriptive statistics table stratified by Berlin (BSS) status --------

bss_table <- df_cleaned |>
  mutate(
    gender = factor(
      gender,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    liver_diagnosis = factor(
      liver_diagnosis,
      levels = c(1, 2, 3, 4, 5),
      labels = c(
        "Hepatitis C",
        "Hepatitis B",
        "PSC/PBC/AIH",
        "Alcohol-related",
        "Other"
      )
    ),
    
    bss_score = factor(
      bss_score,
      levels = c(0, 1),
      labels = c(
        "Low likelihood",
        "High likelihood"
      )
    ),
    
    across(all_of(
      binary_clinical_variables
    ),
    ~ factor(
      .x,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    )
  ) |>
  tbl_summary(
    by = bss_score,
    
    include = all_of(c(
      demographic_variables,
      clinical_variables
    )
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    
    value = list(
      all_of(binary_clinical_variables) ~ "Yes"
    ),
    
    label = list(
      age = "Age",
      gender = "Gender",
      bmi = "Body Mass Index (BMI)",
      liver_diagnosis = "Liver diagnosis",
      recurrence_of_disease = "Recurrence of disease",
      rejection_graft_dysfunction = "Rejection or graft dysfunction",
      any_fibrosis = "Any fibrosis (grade A2 and higher)",
      renal_failure = "Renal failure",
      depression = "Depression",
      corticoid = "Corticosteroid use",
      time_from_transplant = "Time from transplant (years)"
    ),
    
    missing = "no"
  )

# Merge stratified descriptive statistics tables -------------------------------

stratified_table <- tbl_merge(
  tbls = list(psqi_table, ess_table, ais_table, bss_table),
  tab_spanner = c(
    "**PSQI**",
    "**ESS**",
    "**AIS**",
    "**BSS**"
  )
)

stratified_gt <- stratified_table |>
  as_gt() |>
  gt::tab_options(
    table.layout = "auto",
    table.font.size = gt::px(12)
  )

gtsave(
  stratified_gt,
  filename = "tbl_characteristics_stratified.png",
  vwidth = 2400,
  vheight = 1800,
  zoom = 1.5,
  expand = 20
)

# Create sleep disturbance prevalence statistics tables ------------------------

sleep_disturbance_table <- df_cleaned |>
  mutate(
    across(all_of(
      sleep_measure_binary
    ),
    ~ factor(
      .x,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    )
  ) |>
  tbl_summary(

    include = all_of(
      sleep_measure_binary
    ),
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = all_continuous() ~ 2,
    
    value = list(
      all_of(sleep_measure_binary) ~ "Yes"
    ),
    
    label = list(
      psqi_binary = "Pittsburgh Sleep Quality Index",
      ess_binary = "Epworth Sleepiness Scale",
      ais_binary = "Athens Insomnia Scale",
      bss_score = "Berlin Sleepiness Scale"
    ),
    
    missing = "no"
  ) |>
  
  modify_header(
    stat_0 ~ "**Prevalence**"
  ) |>
  
  add_n() |>
  
  add_ci(
    pattern = NULL,
    statistic = all_categorical() ~ "[{conf.low}%, {conf.high}%]",
    style_fun = list(
      all_categorical() ~ label_style_number(
        digits = 1,
        scale = 100
      )
    )
  ) |>
  
  modify_footnote_body(
    columns = "label",
    rows = label == "Pittsburgh Sleep Quality Index",
    footnote = "PSQI > 4"
  ) |>
  
  modify_footnote_body(
    columns = "label",
    rows = label == "Epworth Sleepiness Scale",
    footnote = "ESS > 10"
  ) |>
  
  modify_footnote_body(
    columns = "label",
    rows = label == "Athens Insomnia Scale",
    footnote = "AIS > 5"
  ) |>
  
  modify_footnote_body(
    columns = "label",
    rows = label == "Berlin Sleepiness Scale",
    footnote = "BSS: high likelihood of sleep-disordered breathing"
  )
  
sleep_disturbance_table |>
  as_gt() |>
  gtsave("tbl_sleep_disturbance.png")

# ==============================================================================
# RQ1A - Prevalence
# ==============================================================================

# Prevalence with proper 95% CIs (Wilson score interval) -----------------------

prevalence_ci <- function(x_vec, label) {
  x <- sum(x_vec, na.rm = TRUE)
  n <- sum(!is.na(x_vec))
  ci <- binom.confint(x, n, methods = "wilson")
  data.frame(
    instrument = label, positive = x, n = n,
    prevalence_percentage = round(100 * ci$mean, 1),
    ci_lower = round(100 * ci$lower, 1),
    ci_upper = round(100 * ci$upper, 1)
  )
}

prevalence_table <- bind_rows(
  prevalence_ci(df_cleaned$ess_binary, "ESS > 10 (daytime sleepiness)"),
  prevalence_ci(df_cleaned$psqi_binary, "PSQI > 4 (poor sleep quality)"),
  prevalence_ci(df_cleaned$ais_binary, "AIS > 5 (insomnia)"),
  prevalence_ci(df_cleaned$bss_score, "Berlin (high SDB risk)")
)
print(prevalence_table)
write.csv(prevalence_table, "prevalence_estimates.csv", row.names = FALSE)

# Prevalence by sleep instrument plot ------------------------------------------

prevalence_table <- prevalence_table |>
  mutate(
    instrument = factor(
      instrument,
      levels = instrument[order(-prevalence_percentage)]
    )
  )

fig1_prevalence <- ggplot(
  prevalence_table,
  aes(
    x = instrument,
    y = prevalence_percentage,
    fill = instrument
  )
) +
  
  geom_col(
    width = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.15,
    linewidth = 0.3,
    colour = "grey40"
  ) +

  scale_fill_paletteer_d(
    "fishualize::Epinephelus_lanceolatus"
  ) +
  
  scale_x_discrete(
    labels = c(
      "PSQI > 4 (poor sleep quality)" = "PSQI",
      "AIS > 5 (insomnia)" = "AIS",
      "Berlin (high SDB risk)" = "BSS",
      "ESS > 10 (daytime sleepiness)" = "ESS"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    x = "Sleep Instrument",
    y = "Prevalence (%)",
    caption = paste0(
      "Cut-offs: PSQI > 4; ESS > 10; AIS > 5.\n",
      "BSS = high likelihood of sleep-disordered breathing.\n",
      "Error bars represent 95% Wilson score confidence intervals"
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(
      size = 10,
      colour = "grey20"
    ),
    axis.text = element_text(
      size = 8,
      colour = "grey30"
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 7,
      colour = "grey30",
      margin = margin (t = 8)
    ),
    plot.caption.position = "plot",
    legend.position = "none"
  )

print(fig1_prevalence)

ggsave(
  filename = "fig1_prevalence.png",
  plot = fig1_prevalence,
  width = 5.8,
  height = 3.6,
  units = "in",
  dpi = 300
)

# ==============================================================================
# Univariate measure screen
# ==============================================================================

# Count missing values per sleep instrument
df_cleaned |>
  summarise(across(
    all_of(sleep_measure_binary),
    ~ sum(is.na(.))
  ))

df_cleaned |> tbl_summary(
  include = sleep_measure_scores,
  statistic = list(all_continuous() ~ "{mean} ({sd})")
)

# ==============================================================================
# Binary vs. continuous outcomes
# ==============================================================================

# For RQ1B and RQ2, continuous scores may preserve more information and power
# Compare continuous vs. binary models empirically using model fit

# Compare PSQI representations -------------------------------------------------

m_cont <- lm(sf36_pcs ~ psqi_score, data = df_cleaned)
m_bin  <- lm(sf36_pcs ~ psqi_binary, data = df_cleaned)

summary(m_cont)$r.squared
summary(m_bin)$r.squared

# Lower AIC indicates better fit when models use the same data and outcome
AIC(m_cont, m_bin)

# Load MASS for stepwise selection. Loaded here to avoid MASS::select() masking
# earlier dplyr::select() calls
library(MASS)

# Create complete-case datasets ------------------------------------------------

# StepAIC() requires models being compared to use the same observations. Lock in complete cases for the outcomes and predictors before fitting models.

complete_case_data <- function(formula, data) {
  vars <- all.vars(formula)
  data |> dplyr::select(all_of(vars)) |> drop_na()
}

# ==============================================================================
# RQ1B: selecting predictors for four sleep disturbance instruments
# ==============================================================================

# Multivariable regression: logistic and linear models -------------------------

# Fit multivariable logistic and linear regression models for the four sleep 
# instruments. Keep the full domain-knowledge model as the primary model. Use
# stepAIC-reduced models as a comparison rather than the sole analysis. 
# Agreement between models strengthens evidence; differences should be
# discussed.

# ------------------------------------------------------------------------------
# PSQI 
# ------------------------------------------------------------------------------

# ---- PSQI: complete logistic regression model --------------------------------

# Create logistic regression model for PSQI
psqi_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "psqi_binary"
)

df_psqi_cc <- complete_case_data(psqi_formula, df_cleaned)
m_log_psqi <- glm(psqi_formula, data = df_psqi_cc, family = "binomial")

# Read the raw model results (coefficients on the log-odds scale, p-values)
summary(m_log_psqi)

# Convert log-odds coefficients into odds ratios + 95% CIs (interpretable scale:
# "patients with X are __ times more likely to have PSQI-defined sleep
# disturbance")
psqi_or <- exp(cbind(OR = coef(m_log_psqi), confint(m_log_psqi)))
print(psqi_or)

# Check for multicollinearity among predictors (VIF > ~5 is a concern)
m_psqi_vif <- vif(m_log_psqi)
print(m_psqi_vif)

# Residual diagnostics: look for systematic patterns (points should scatter
# randomly around 0, no obvious trend/curve)
res_psqi_p <- residuals(m_log_psqi, type = "pearson")

plot(fitted(m_log_psqi), res_psqi_p,
  main = "PSQI model: Pearson residuals",
  xlab = "Fitted values",
  ylab = "Pearson residuals"
)

abline(h = 0, lty = 2)

res_psqi_d <- residuals(m_log_psqi, type = "deviance")

plot(fitted(m_log_psqi), res_psqi_d,
  main = "PSQI model: Deviance residuals",
  xlab = "Fitted values",
  ylab = "Deviance residuals"
)

abline(h = 0, lty = 2)

# ---- PSQI: stepAIC logistic regression model ---------------------------------

m_log_psqi_step <- stepAIC(m_log_psqi, direction = "backward", trace = FALSE)
summary(m_log_psqi_step)

# Resulting model: psqi_binary ~ age + gender + recurrence_of_disease +
# renal_failure + depression
exp(cbind(OR = coef(m_log_psqi_step), confint(m_log_psqi_step)))

# ---- PSQI: linear regression model -------------------------------------------

psqi_formula_cont <- reformulate(
  termlabels = c(demographic_variables, clinical_variables),
  response = "psqi_score"
)

df_psqi_cont_cc <- complete_case_data(psqi_formula_cont, df_cleaned)

m_lin_psqi <- lm(psqi_formula_cont, data = df_psqi_cont_cc)

summary(m_lin_psqi)

vif(m_lin_psqi)

par(mfrow = c(2, 2))
plot(m_lin_psqi)
par(mfrow = c(1, 1))

# ---- PSQI: stepAIC linear regression model -----------------------------------

m_lin_psqi_step <- stepAIC(m_lin_psqi, direction = "backward", trace = FALSE)

summary(m_lin_psqi_step)

# Resulting model: psqi_score ~ age + gender + bmi + recurrence_of_disease +
# rejection_graft_dysfunction + any_fibrosis + depression + corticoid

# ---- PSQI: model comparison --------------------------------------------------

# Did the same predictors survive selection on the continuous score vs. the 
# binary threshold for PSQI?

names(coef(m_log_psqi_step))
names(coef(m_lin_psqi_step))

# The shared predictors are age, gender, recurrence of disease, and depression

# ------------------------------------------------------------------------------
# ESS 
# ------------------------------------------------------------------------------

# ---- ESS: complete logistic regression model ---------------------------------

# Create logistic regression model for ESS

ess_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "ess_binary"
)

df_ess_cc <- complete_case_data(ess_formula, df_cleaned)
m_log_ess <- glm(ess_formula, data = df_ess_cc, family = "binomial")

summary(m_log_ess)

ess_or <- exp(cbind(OR = coef(m_log_ess), confint(m_log_ess)))
print(ess_or)

m_ess_vif <- vif(m_log_ess)
print(m_ess_vif)

res_ess_p <- residuals(m_log_ess, type = "pearson")

plot(fitted(m_log_ess), res_ess_p,
  main = "ESS model: Pearson residuals",
  xlab = "Fitted values",
  ylab = "Pearson residuals"
)

abline(h = 0, lty = 2)

res_ess_d <- residuals(m_log_ess, type = "deviance")

plot(fitted(m_log_ess), res_ess_d,
  main = "ESS model: Deviance residuals",
  xlab = "Fitted values",
  ylab = "Deviance residuals"
)

abline(h = 0, lty = 2)

# ---- ESS: stepAIC logistic regression model ----------------------------------

m_log_ess_step <- stepAIC(m_log_ess, direction = "backward", trace = FALSE)

summary(m_log_ess_step)

# Resulting model: ess_binary ~ liver_diagnosis + rejection_graft_dysfunction + 
# renal_failure + corticoid

exp(cbind(OR = coef(m_log_ess_step), confint(m_log_ess_step)))

# ---- ESS: complete linear regression model -----------------------------------

ess_formula_cont <- reformulate(
  termlabels = c(demographic_variables, clinical_variables),
  response = "ess_score"
)

df_ess_cont_cc <- complete_case_data(ess_formula_cont, df_cleaned)

m_lin_ess <- lm(ess_formula_cont, data = df_ess_cont_cc)

summary(m_lin_ess)

vif(m_lin_ess)

par(mfrow = c(2, 2)); plot(m_lin_ess); par(mfrow = c(1, 1))

# ---- ESS: stepAIC linear regression model ------------------------------------

m_lin_ess_step <- stepAIC(m_lin_ess, direction = "backward", trace = FALSE)

summary(m_lin_ess_step)

# Resulting model: ess_score ~ gender + rejection_graft_dysfunction + depression + corticoid

# ---- ESS: model comparison ---------------------------------------------------

names(coef(m_log_ess_step))

names(coef(m_lin_ess_step))

# The shared predictors are rejection graft dysfunction and corticoid

# ------------------------------------------------------------------------------
# BSS 
# ------------------------------------------------------------------------------

# No continuous scores exist for the BSS so only a logistic regression model
# will be used

# ---- BSS: complete logistic regression model ---------------------------------

# Create logistic regression model for BSS
bss_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "bss_score"
)

df_bss_cc <- complete_case_data(bss_formula, df_cleaned)
m_log_bss <- glm(bss_formula, data = df_bss_cc, family = "binomial")

summary(m_log_bss)

bss_or <- exp(cbind(OR = coef(m_log_bss), confint(m_log_bss)))
print(bss_or)

m_bss_vif <- vif(m_log_bss)
print(m_bss_vif)

res_bss_p <- residuals(m_log_bss, type = "pearson")

plot(fitted(m_log_bss), res_bss_p,
     main = "BSS model: Pearson residuals",
     xlab = "Fitted values", ylab = "Pearson residuals"
)

abline(h = 0, lty = 2)

res_bss_d <- residuals(m_log_bss, type = "deviance")

plot(fitted(m_log_bss), res_bss_d,
     main = "BSS model: Deviance residuals",
     xlab = "Fitted values", ylab = "Deviance residuals"
)

abline(h = 0, lty = 2)

# ---- BSS: stepAIC logistic regression model ----------------------------------

m_log_bss_step <- stepAIC(m_log_bss, direction = "backward", trace = FALSE)

summary(m_log_bss_step)

# Resulting model: bss_score ~ bmi + recurrence_of_disease + renal_failure
exp(cbind(OR = coef(m_log_bss_step), confint(m_log_bss_step)))

# ------------------------------------------------------------------------------
# AIS 
# ------------------------------------------------------------------------------

# ---- AIS: complete logistic regression model ---------------------------------

# Create logistic regression model for AIS
ais_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "ais_binary"
)

df_ais_cc <- complete_case_data(ais_formula, df_cleaned)
m_log_ais <- glm(ais_formula, data = df_ais_cc, family = "binomial")

summary(m_log_ais)

ais_or <- exp(cbind(OR = coef(m_log_ais), confint(m_log_ais)))
print(ais_or)

m_ais_vif <- vif(m_log_ais)
print(m_ais_vif)

res_ais_p <- residuals(m_log_ais, type = "pearson")
plot(fitted(m_log_ais), res_ais_p,
     main = "AIS model: Pearson residuals",
     xlab = "Fitted values", ylab = "Pearson residuals"
)
abline(h = 0, lty = 2)

res_ais_d <- residuals(m_log_ais, type = "deviance")
plot(fitted(m_log_ais), res_ais_d,
     main = "AIS model: Deviance residuals",
     xlab = "Fitted values", ylab = "Deviance residuals"
)
abline(h = 0, lty = 2)

# ---- AIS: stepAIC logistic regression model ----------------------------------

m_log_ais_step <- stepAIC(m_log_ais, direction = "backward", trace = FALSE)

summary(m_log_ais_step)

# Resulting model: is_binary ~ age + recurrence_of_disease + depression +
# corticoid
exp(cbind(OR = coef(m_log_ais_step), confint(m_log_ais_step)))

# ---- BSS: complete linear regression model -----------------------------------

ais_formula_cont <- reformulate(
  termlabels = c(demographic_variables, clinical_variables),
  response = "ais_score"
)
df_ais_cont_cc <- complete_case_data(ais_formula_cont, df_cleaned)

m_lin_ais <- lm(ais_formula_cont, data = df_ais_cont_cc)

summary(m_lin_ais)

vif(m_lin_ais)

par(mfrow = c(2, 2)); plot(m_lin_ais); par(mfrow = c(1, 1))

# ---- AIS: stepAIC linear regression model ------------------------------------

m_lin_ais_step <- stepAIC(m_lin_ais, direction = "backward", trace = FALSE)
summary(m_lin_ais_step)

# Resulting model: ais_score ~ age + bmi + recurrence_of_disease + corticoid

# ---- AIS: model comparison ---------------------------------------------------

names(coef(m_log_ais_step))
names(coef(m_lin_ais_step))

# The shared predictors are age, reccurence of disease, and corticoid

# Which predictors survived AIC-based selection vs. the full model, per
# instrument? Divergence here is exactly what the assignment prompt
# anticipates ("some predictors better associated with some measures").
AIC(m_log_psqi, m_log_psqi_step)
AIC(m_log_ess,  m_log_ess_step)
AIC(m_log_bss,  m_log_bss_step)
AIC(m_log_ais,  m_log_ais_step)

# All the stepwise selected models have less predictors and smaller AIC values

# ==============================================================================
# RQ2: Relationship between sleep disturbance and quality of life
# ==============================================================================

# ------------------------------------------------------------------------------
# SF36-PCS: physical quality of life
# ------------------------------------------------------------------------------

# ---- PCS: unadjusted multivariable linear regression models ------------------

# Examine unadjusted associations between sleep scores and SoL. These models do
# not account for potential confounding patient characteristics

# Create linear regression model for physical quality of life
m_lin_pcs <- lm(sf36_pcs ~ psqi_score + ess_score + bss_score + ais_score,
                data = df_cleaned
)
summary(m_lin_pcs)

pcs_formula_adj <- reformulate(
  termlabels = c(
    "psqi_score", "ess_score", "bss_score", "ais_score",
    demographic_variables,
    clinical_variables
  ),
  response = "sf36_pcs"
)

df_pcs_adj_cc <- complete_case_data(pcs_formula_adj, df_cleaned)

m_lin_pcs_adj <- lm(pcs_formula_adj, data = df_pcs_adj_cc)

summary(m_lin_pcs_adj)

m_pcs_adj_vif <- vif(m_lin_pcs_adj)

print(m_pcs_adj_vif)

# ---- PCS: stepAIC multivariable linear regression models ---------------------

m_lin_pcs_adj_step <- stepAIC(m_lin_pcs_adj, direction = "backward", trace = FALSE)

summary(m_lin_pcs_adj_step)

# resulting model: sf36_pcs ~ ess_score + ais_score + age + bmi + 
# recurrence_of_disease

# ------------------------------------------------------------------------------
# SF36-MCS: mental quality of life
# ------------------------------------------------------------------------------

# ---- MCS: unadjusted multivariable linear regression models ------------------

# Create linear regression model for mental quality of life
m_lin_mcs <- lm(sf36_mcs ~ psqi_score + ess_score + bss_score + ais_score,
                data = df_cleaned
)
summary(m_lin_mcs)

par(mfrow = c(2, 2)); plot(m_lin_pcs_adj); par(mfrow = c(1, 1))

mcs_formula_adj <- reformulate(
  termlabels = c(
    "psqi_score", "ess_score", "bss_score", "ais_score",
    demographic_variables,
    clinical_variables
  ),
  response = "sf36_mcs"
)

df_mcs_adj_cc <- complete_case_data(mcs_formula_adj, df_cleaned)

m_lin_mcs_adj <- lm(mcs_formula_adj, data = df_mcs_adj_cc)

summary(m_lin_mcs_adj)

m_mcs_adj_vif <- vif(m_lin_mcs_adj)

print(m_mcs_adj_vif)

# ---- MCS: stepAIC multivariable linear regression models ---------------------

m_lin_mcs_adj_step <- stepAIC(m_lin_mcs_adj, direction = "backward", trace = FALSE)

summary(m_lin_mcs_adj_step)

# Resulting model: sf36_mcs ~ psqi_score + ais_score + age + liver_diagnosis + 
# rejection_graft_dysfunction + any_fibrosis + depression + time_from_transplant

# ------------------------------------------------------------------------------
# SF36-PCS and SF36-MCS: residual diagnostics
# ------------------------------------------------------------------------------

par(mfrow = c(2, 2)); plot(m_lin_mcs_adj); par(mfrow = c(1, 1))

# Compare full and stepwise models
AIC(m_lin_pcs_adj, m_lin_pcs_adj_step)
AIC(m_lin_mcs_adj, m_lin_mcs_adj_step)

# Lower AIC indicates better relative fit. Check which sleep scores remain after
# stepAIC selection. If only AIS remains, the other sleep measures may add
# little independent explanatory value once the remaining predictors are
# included.

# ------------------------------------------------------------------------------
# SF36-PCS and SF36-MCS: unadjusted vs. adjusted models
# ------------------------------------------------------------------------------

summary(m_lin_pcs)$r.squared #0.207
summary(m_lin_pcs_adj)$r.squared #0.292
summary(m_lin_pcs_adj_step)$r.squared #0.245

summary(m_lin_mcs)$r.squared #0.319
summary(m_lin_mcs_adj)$r.squared #0.474
summary(m_lin_mcs_adj_step)$r.squared #0.470


# ==============================================================================
# Sensitivity analysis: multiple imputation vs. complete-case for PSQI
# ==============================================================================

# PSQI has substantial missingness (85/268, 32%). Compare multiple-imputation
# results with complete-case results to assess whether PSQI missingness
# meaningfully affects the conclusions.

# This sensitivity analysis uses continuous PSQI scores only. Binary PSQI would
# require rederiving PSQI_binary after imputation.

# Prepare data for imputation --------------------------------------------------

imp_vars <- c(
  demographic_variables,
  clinical_variables,
  "psqi_score",
  "ess_score",
  "ais_score",
  "bss_score",
  "sf36_pcs",
  "sf36_mcs"
)

df_for_imputation <- df_cleaned |>
  dplyr::select(
    all_of(
      imp_vars
    )
  )

# Perform multiple imputation --------------------------------------------------

# Use 20 imputations and predictive mean matching (PMM) for numeric variables.
# PMM uses observed donor values, helping keep scores within plausible ranges.

imp <- mice(
  df_for_imputation,
  m = 20,
  method = "pmm",
  seed = 431859,
  print = FALSE
)

# RQ1B: PSQI as a continuous outcome -------------------------------------------

# Fit the same model in each inputed dataset and pool estimates using Rubin's
# rules.
fits_psqi_mi <- lapply(
  seq_len(imp$m),
  function(i) {
    lm(psqi_formula_cont,
       data = complete(imp, i))
  }
)

pooled_psqi <- pool(as.mira(fits_psqi_mi))

summary(pooled_psqi)

# Compare PSQI complete-case vs. imputed reults --------------------------------
summary(m_lin_psqi)$coefficients
summary(pooled_psqi)

# RQ2: PSQI as a predictor of physical quality of life -------------------------

fits_pcs_mi <- lapply(seq_len(imp$m), function(i) {
  lm(pcs_formula_adj, data = complete(imp, i))
})

pooled_pcs <- pool(as.mira(fits_pcs_mi))

summary(pooled_pcs)

# RQ2: PSQI as a predictor of mental quality of life ---------------------------

fits_mcs_mi <- lapply(seq_len(imp$m), function(i) {
  lm(mcs_formula_adj, data = complete(imp, i))
})

pooled_mcs <- pool(as.mira(fits_mcs_mi))

summary(pooled_mcs)

# Compare complete-case vs. imputed results ------------------------------------

# Similar coefficients, CIs and conclusions would support robustness of the 
# complete-case findings to PSQI missingness

summary(m_lin_pcs_adj)$coefficients
summary(pooled_pcs)

summary(m_lin_mcs_adj)$coefficients
summary(pooled_mcs)

# ---- Conclusion from this sensitivity check ----
#
# RQ2 (sleep-QoL relationship) is robust to how PSQI's missingness is handled
#   - psqi_score -> MCS: coefficient values are similar with both being significant
#     (CC: beta=-0.83, p=0.007; MI: beta=-0.79, p=0.005).
#   - psqi_score -> PCS: not significant in either analysis (CC: p=0.62;
#     MI: p=0.94) 
#   Conclusion: PSQI's relationship with mental QoL but not physical QoL
#   is not an artifact of dropping the 85 PSQI-missing patients. Report
#   complete-case as primary for RQ2, and cite this MI check as
#   confirming robustness.
#
# RQ1b (predictors of PSQI itself) is NOT robust -- several predictors
# change conclusions between complete-case and MI:
#   - rejection_graft_dysfunction: significant in CC (p=0.018,
#     beta=-2.52), not significant under MI (p=0.37, beta=-0.82 --
#     shrinks to about a third).
#   - corticoid: flips from significant (CC p=0.023) to borderline
#     (MI p=0.056).
#   - depression: flips the other way, from borderline (CC p=0.063) to
#     significant (MI p=0.031).


# ==============================================================================
# Result tables for the regression analyses (RQ1b + RQ2 + sensitivity check)
# ==============================================================================
# Append this to the end of full_analysis.R and run it AFTER the rest of that
# script -- it does not refit any of the original models or change any of
# their estimates. It only builds publication-style tables from the model
# objects already created above (m_log_psqi, m_log_psqi_step, m_lin_pcs_adj,
# pooled_psqi, etc.).
#
# One exception: for display purposes only, RQ1b/RQ2 models are refit on a
# copy of the data with factor levels LABELED (e.g. gender "1"/"2" ->
# "Male"/"Female", matching the labels already used in the descriptive
# tables earlier in the script). This changes nothing about the numbers --
# same rows, same predictors, same coefficients -- it just makes the tables
# readable instead of showing raw numeric codes.
#
# NOTE: this was written and proofread without an R session available to
# execute it (this sandbox has no R installation and no network access to
# install one) -- run it in your own R session and flag anything that
# errors so it can be fixed.

# ---- Shared predictor labels, reused across every table --------------------

var_labels <- list(
  age                          ~ "Age",
  gender                       ~ "Gender",
  bmi                          ~ "Body Mass Index (BMI)",
  liver_diagnosis              ~ "Liver diagnosis",
  recurrence_of_disease        ~ "Recurrence of disease",
  rejection_graft_dysfunction  ~ "Rejection or graft dysfunction",
  any_fibrosis                 ~ "Any fibrosis (grade A2 and higher)",
  renal_failure                ~ "Renal failure",
  depression                   ~ "Depression",
  corticoid                    ~ "Corticosteroid use",
  time_from_transplant         ~ "Time from transplant",
  psqi_score                   ~ "PSQI score",
  ess_score                    ~ "ESS score",
  bss_score                    ~ "BSS (high SDB risk)",
  ais_score                    ~ "AIS score"
)

# tbl_regression()'s `label` argument errors (rather than silently ignoring
# entries) if it's given a label for a variable that isn't actually a
# PREDICTOR in the model -- e.g. passing the full var_labels list (which
# includes psqi_score, ess_score, bss_score, ais_score as predictor labels
# for the RQ2 models) to a model where one of those is instead the RESPONSE
# (m_lin_psqi/m_lin_ess/m_lin_ais/m_log_bss all have one of these four score
# variables as their outcome, not a predictor). Using all.vars(formula(model))
# alone doesn't catch this, since that includes the response too --
# delete.response() strips it first, so only genuine predictor names are
# checked against var_labels.
labels_for_model <- function(model) {
  predictor_vars <- all.vars(delete.response(terms(model)))
  Filter(function(f) all.vars(f)[1] %in% predictor_vars, var_labels)
}

# Logistic models -> odds ratios; linear models -> beta coefficients
or_tbl <- function(model) {
  tbl_regression(model, exponentiate = TRUE, label = labels_for_model(model)) |>
    bold_p() |>
    modify_header(label = "**Predictor**")
}

beta_tbl <- function(model) {
  tbl_regression(model, label = labels_for_model(model)) |>
    bold_p() |>
    modify_header(label = "**Predictor**")
}

# ---- Labeled copy of the data, for table display only ----------------------
# Coefficients don't change when a factor's *labels* change (only the
# reference level would) -- this just makes tables print "Female" /
# "Hepatitis B" / "Yes" instead of raw codes "2" / "2" / "1".

label_clinical_vars <- function(data) {
  data |>
    mutate(
      gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
      liver_diagnosis = factor(
        liver_diagnosis,
        levels = c(1, 2, 3, 4, 5),
        labels = c("Hepatitis C", "Hepatitis B", "PSC/PBC/AIH", "Alcohol-related", "Other")
      ),
      across(
        all_of(binary_clinical_variables),
        ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))
      )
    )
}

df_cleaned_lbl <- label_clinical_vars(df_cleaned)

# Refits a formula's full model and its already-derived stepwise model on ONE
# locked-in complete-case subset of the labeled data -- same logic as
# complete_case_data() above, so the stepwise refit isn't accidentally fit on
# a different (larger) row subset than the full model it's being compared to.
refit_pair_for_table <- function(full_formula, step_model, data = df_cleaned_lbl,
                                 family = NULL) {
  cc <- complete_case_data(full_formula, data)
  fit_one <- function(f) {
    if (is.null(family)) lm(f, data = cc) else glm(f, data = cc, family = family)
  }
  list(full = fit_one(full_formula), step = fit_one(formula(step_model)))
}

merge_pair <- function(pair, caption, exponentiate) {
  tbl_fn <- if (exponentiate) or_tbl else beta_tbl
  tbl_merge(
    tbls = list(tbl_fn(pair$full), tbl_fn(pair$step)),
    tab_spanner = c("**Full model**", "**Stepwise-selected model**"),
    quiet = TRUE
  ) |>
    modify_caption(caption)
}

# ==============================================================================
# RQ1b tables: predictors of each sleep-disturbance instrument
# (full domain-knowledge model vs. AIC stepwise-selected model)
# ==============================================================================

psqi_log <- refit_pair_for_table(psqi_formula, m_log_psqi_step, family = "binomial")
psqi_lin <- refit_pair_for_table(psqi_formula_cont, m_lin_psqi_step)
ess_log  <- refit_pair_for_table(ess_formula, m_log_ess_step, family = "binomial")
ess_lin  <- refit_pair_for_table(ess_formula_cont, m_lin_ess_step)
bss_log  <- refit_pair_for_table(bss_formula, m_log_bss_step, family = "binomial")
ais_log  <- refit_pair_for_table(ais_formula, m_log_ais_step, family = "binomial")
ais_lin  <- refit_pair_for_table(ais_formula_cont, m_lin_ais_step)

tbl_psqi_log <- merge_pair(
  psqi_log,
  "**Table. Predictors of PSQI-defined sleep disturbance (logistic regression, OR [95% CI])**",
  exponentiate = TRUE
)
tbl_psqi_lin <- merge_pair(
  psqi_lin,
  "**Table. Predictors of PSQI score (linear regression, β [95% CI])**",
  exponentiate = FALSE
)
tbl_ess_log <- merge_pair(
  ess_log,
  "**Table. Predictors of ESS-defined excessive daytime sleepiness (logistic regression, OR [95% CI])**",
  exponentiate = TRUE
)
tbl_ess_lin <- merge_pair(
  ess_lin,
  "**Table. Predictors of ESS score (linear regression, β [95% CI])**",
  exponentiate = FALSE
)
tbl_bss_log <- merge_pair(
  bss_log,
  "**Table. Predictors of high-risk sleep-disordered breathing, Berlin Questionnaire (logistic regression, OR [95% CI])**",
  exponentiate = TRUE
)
tbl_ais_log <- merge_pair(
  ais_log,
  "**Table. Predictors of AIS-defined insomnia (logistic regression, OR [95% CI])**",
  exponentiate = TRUE
)
tbl_ais_lin <- merge_pair(
  ais_lin,
  "**Table. Predictors of AIS score (linear regression, β [95% CI])**",
  exponentiate = FALSE
)

tbl_psqi_log %>%
  as_gt() %>%
  gtsave("tbl_psqi_log.png")
tbl_psqi_lin %>%
  as_gt() %>%
  gtsave("tbl_psqi_lin.png")
tbl_ess_log %>%
  as_gt() %>%
  gtsave("tbl_ess_log.png")
tbl_ess_lin %>%
  as_gt() %>%
  gtsave("tbl_ess_lin.png")
tbl_bss_log %>%
  as_gt() %>%
  gtsave("tbl_bss_log.png")
tbl_ais_log %>%
  as_gt() %>%
  gtsave("tbl_ais_log.png")
tbl_ais_lin %>%
  as_gt() %>%
  gtsave("tbl_ais_lin.png")

# ==============================================================================
# RQ2 tables: sleep disturbance and quality of life
# (crude vs. adjusted vs. adjusted+stepwise, for PCS and MCS)
# ==============================================================================

df_pcs_adj_cc_lbl <- complete_case_data(pcs_formula_adj, df_cleaned_lbl)
df_mcs_adj_cc_lbl <- complete_case_data(mcs_formula_adj, df_cleaned_lbl)

m_lin_pcs_crude_tbl <- lm(sf36_pcs ~ psqi_score + ess_score + bss_score + ais_score,
                          data = df_cleaned_lbl)
m_lin_pcs_adj_tbl   <- lm(pcs_formula_adj, data = df_pcs_adj_cc_lbl)
m_lin_pcs_step_tbl  <- lm(formula(m_lin_pcs_adj_step), data = df_pcs_adj_cc_lbl)

m_lin_mcs_crude_tbl <- lm(sf36_mcs ~ psqi_score + ess_score + bss_score + ais_score,
                          data = df_cleaned_lbl)
m_lin_mcs_adj_tbl   <- lm(mcs_formula_adj, data = df_mcs_adj_cc_lbl)
m_lin_mcs_step_tbl  <- lm(formula(m_lin_mcs_adj_step), data = df_mcs_adj_cc_lbl)

tbl_pcs <- tbl_merge(
  tbls = list(
    beta_tbl(m_lin_pcs_crude_tbl),
    beta_tbl(m_lin_pcs_adj_tbl),
    beta_tbl(m_lin_pcs_step_tbl)
  ),
  tab_spanner = c("**Crude**", "**Adjusted**", "**Adjusted, stepwise**"),
  quiet = TRUE
) |>
  modify_caption("**Table. Sleep disturbance and physical quality of life (SF-36 PCS, β [95% CI])**")

tbl_mcs <- tbl_merge(
  tbls = list(
    beta_tbl(m_lin_mcs_crude_tbl),
    beta_tbl(m_lin_mcs_adj_tbl),
    beta_tbl(m_lin_mcs_step_tbl)
  ),
  tab_spanner = c("**Crude**", "**Adjusted**", "**Adjusted, stepwise**"),
  quiet = TRUE
) |>
  modify_caption("**Table. Sleep disturbance and mental quality of life (SF-36 MCS, β [95% CI])**")

tbl_pcs %>%
  as_gt() %>%
  gtsave("tbl_pcs.png")
tbl_mcs %>%
  as_gt() %>%
  gtsave("tbl_mcs.png")

# ==============================================================================
# Sensitivity tables: complete-case vs. multiple imputation, for PSQI
# ==============================================================================
# Pooled mice::mipo objects (pooled_psqi, pooled_pcs, pooled_mcs, built in the
# MI section above) don't carry the same categorical-variable metadata
# tbl_regression()/tbl_merge() rely on, so these comparisons are built
# directly from broom::tidy() (complete-case) and summary(..., conf.int =
# TRUE) (pooled) instead. term_labels below just prettifies the raw dummy-
# variable names (e.g. "gender2" -> "Gender: Female") using the same coding
# scheme as the un-labeled models used in the MI section (df_cleaned, not
# df_cleaned_lbl -- the imputation was run on the original coding).

term_labels <- c(
  age = "Age", gender2 = "Gender: Female", bmi = "BMI",
  liver_diagnosis2 = "Liver diagnosis: Hepatitis B",
  liver_diagnosis3 = "Liver diagnosis: PSC/PBC/AIH",
  liver_diagnosis4 = "Liver diagnosis: Alcohol-related",
  liver_diagnosis5 = "Liver diagnosis: Other",
  recurrence_of_disease1 = "Recurrence of disease: Yes",
  rejection_graft_dysfunction1 = "Rejection/graft dysfunction: Yes",
  any_fibrosis1 = "Any fibrosis: Yes",
  renal_failure1 = "Renal failure: Yes",
  depression1 = "Depression: Yes",
  corticoid1 = "Corticosteroid use: Yes",
  time_from_transplant = "Time from transplant",
  psqi_score = "PSQI score", ess_score = "ESS score",
  bss_score = "BSS (high SDB risk)", ais_score = "AIS score"
)

pretty_term <- function(term) {
  dplyr::if_else(term %in% names(term_labels), unname(term_labels[term]), term)
}

cc_mi_compare <- function(cc_model, pooled_model, digits = 2) {
  cc <- broom::tidy(cc_model, conf.int = TRUE) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::transmute(
      term = pretty_term(term),
      `CC beta` = round(estimate, digits),
      `CC 95% CI` = paste0("(", round(conf.low, digits), ", ", round(conf.high, digits), ")"),
      `CC p` = signif(p.value, 3)
    )
  
  mi <- summary(pooled_model, conf.int = TRUE) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::transmute(
      term = pretty_term(term),
      `MI beta` = round(estimate, digits),
      `MI 95% CI` = paste0("(", round(`2.5 %`, digits), ", ", round(`97.5 %`, digits), ")"),
      `MI p` = signif(p.value, 3)
    )
  
  dplyr::full_join(cc, mi, by = "term")
}

tbl_sens_psqi <- cc_mi_compare(m_lin_psqi, pooled_psqi) |>
  gt::gt() |>
  gt::tab_header(title = "Sensitivity check: PSQI score predictors, complete-case vs. MI")

tbl_sens_pcs <- cc_mi_compare(m_lin_pcs_adj, pooled_pcs) |>
  gt::gt() |>
  gt::tab_header(title = "Sensitivity check: adjusted PCS model, complete-case vs. MI")

tbl_sens_mcs <- cc_mi_compare(m_lin_mcs_adj, pooled_mcs) |>
  gt::gt() |>
  gt::tab_header(title = "Sensitivity check: adjusted MCS model, complete-case vs. MI")

tbl_sens_psqi %>%
  gtsave("tbl_sens_psqi.png")
tbl_sens_pcs %>%
  gtsave("tbl_sens_pcs.png")
tbl_sens_mcs %>%
  gtsave("tbl_sens_mcs.png")

# ==============================================================================
# Figure 2: sleep scores by significant clinical predictors (boxplots)
# ==============================================================================

box_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11))

p1 <- ggplot(df_cleaned_lbl, aes(x = factor(bss_score, labels = c("Low risk", "High risk")), y = bmi, fill = bss_score == 1)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("FALSE" = "#B0B7BF", "TRUE" = "#2E5A87")) +
  labs(x = "Berlin Sleepiness Scale (BSS)", y = "BMI (kg/m²)", title = "A. BMI by Berlin sleep apnea risk") +
  box_theme

p2 <- ggplot(df_cleaned_lbl, aes(x = depression, y = psqi_score, fill = depression)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("No" = "#B0B7BF", "Yes" = "#2E5A87")) +
  labs(x = "Depression", y = "PSQI score", title = "B. PSQI score by depression status") +
  box_theme

p3 <- ggplot(df_cleaned_lbl, aes(x = corticoid, y = ess_score, fill = corticoid)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("No" = "#B0B7BF", "Yes" = "#2E5A87")) +
  labs(x = "Corticosteroid use", y = "ESS score", title = "C. ESS score by corticosteroid use") +
  box_theme

p4 <- ggplot(df_cleaned_lbl, aes(x = recurrence_of_disease, y = ais_score, fill = recurrence_of_disease)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("No" = "#B0B7BF", "Yes" = "#2E5A87")) +
  labs(x = "Recurrence of disease", y = "AIS score", title = "D. AIS score by disease recurrence") +
  box_theme

figure2 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Figure 2. Distribution of sleep instrument scores by significant clinical predictors",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave("figure2_clinical_predictors.png", figure2, width = 10, height = 8, dpi = 300)


# ==============================================================================
# Figure 3: forest plot of odds ratios, full logistic models, all 4 instruments
# ==============================================================================

# pull OR + CI from each full model, tag with instrument name, keep only the
# 6 shared demographic/clinical predictors (renal_failure dropped -- separation)
keep_terms <- c("age", "gender2", "bmi", "recurrence_of_disease1", "depression1", "corticoid1")
term_labels <- c(
  age = "Age", gender2 = "Gender: Female", bmi = "BMI",
  recurrence_of_disease1 = "Recurrence of disease",
  depression1 = "Depression", corticoid1 = "Corticosteroid use"
)

get_or <- function(model, instrument) {
  tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term %in% keep_terms) |>
    mutate(instrument = instrument, term = term_labels[term])
}

forest_data <- bind_rows(
  get_or(m_log_psqi, "PSQI"),
  get_or(m_log_ess, "ESS"),
  get_or(m_log_ais, "AIS"),
  get_or(m_log_bss, "BSS")
) |>
  mutate(
    term = factor(term, levels = rev(unname(term_labels))),
    instrument = factor(instrument, levels = c("PSQI", "ESS", "AIS", "BSS")),
    significant = p.value < 0.05
  )

figure3 <- ggplot(forest_data, aes(x = estimate, y = instrument, xmin = conf.low, xmax = conf.high, color = instrument)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(size = significant), fatten = 3) +
  scale_size_manual(values = c("TRUE" = 1.1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = c(PSQI = "#2E5A87", ESS = "#5B9BD5", AIS = "#8FBFE0", BSS = "#B0B7BF")) +
  scale_x_log10() +
  facet_grid(rows = vars(term), scales = "free_y", switch = "y") +
  labs(
    x = "Odds ratio (log scale), 95% CI", y = NULL, color = NULL,
    title = "Figure 3. Odds ratios for key predictors of sleep disturbance across four instruments",
    subtitle = "Larger points indicate p < 0.05. Renal failure excluded (unstable due to separation)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.placement = "outside",
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

ggsave("figure3_forest_plot_predictors.png", figure3, width = 8.5, height = 8, dpi = 300)

