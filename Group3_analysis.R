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
  dplyr::select(
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

sleep_measure_scores_all <- c(
  "psqi_score",
  "ess_score",
  "ais_score",
  "bss_score"
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

# Identify predictors for logistic model tables --------------------------------

# Predictors that can be included in the logistic models. Renal failure is
# disregarded
logistic_predictors <- c(
  demographic_variables,
  setdiff(
    clinical_variables,
    "renal_failure"
  )
)

# ------------------------------------------------------------------------------
# PSQI 
# ------------------------------------------------------------------------------

# ---- PSQI: complete logistic regression model --------------------------------

# Create logistic regression model for PSQI
psqi_formula <- reformulate(
  termlabels = logistic_predictors,
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
  termlabels = logistic_predictors,
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
  termlabels = logistic_predictors,
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
  termlabels = logistic_predictors,
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

# ---- AIS: complete linear regression model -----------------------------------

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

# ---- PCS: sleep-measure-only multivariable linear regression models ----------

# Examine sleep-measure-only associations between sleep scores and QoL. These
# models do not account for potential confounding patient characteristics

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

# ---- MCS: sleep-measure-only multivariable linear regression models ----------

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

# ==============================================================================
# Main regression tables
# ==============================================================================

# Define labels for tables -----------------------------------------------------

main_predictor_labels <- list(
  age = "Age",
  gender = "Gender",
  bmi = "Body mass index (BMI)",
  liver_diagnosis = "Liver diagnosis",
  recurrence_of_disease = "Recurrence of disease",
  rejection_graft_dysfunction = "Rejection or graft dysfunction",
  any_fibrosis = "Any fibrosis (grade A2 and higher)",
  renal_failure = "Renal failure",
  depression = "Depression",
  corticoid = "Corticosteroid use",
  time_from_transplant = "Time from transplant (years)"
)

# Create reusable liver function relabeler function ----------------------------

relabel_liver_diagnosis <- function(table) {
  table |>
    modify_table_body(
      ~ .x |>
        mutate(
          label = case_when(
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "1" ~ "Hepatitis C",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "2" ~ "Hepatitis B",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "3" ~ "PSC/PBC/AIH",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "4" ~ "Alcohol-related",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "5" ~ "Other",
            
            TRUE ~ label
          )
        )
    )
}

# Create reusable table function -----------------------------------------------

make_main_regression_table <- function(
    model,
    model_type = c("logistic", "linear")
) {
  
  model_type <- match.arg(model_type)
  
  is_logistic <- model_type == "logistic"
  
  regression_table <- model |>
    tbl_regression(
      exponentiate = is_logistic,
      
      label = main_predictor_labels,
      
      show_single_row = any_of(
        c(
          "gender",
          binary_clinical_variables
        )
      ),
      
      estimate_fun = if (is_logistic) {
        label_style_ratio(digits = 2)
      } else {
        label_style_number(digits = 2)
      },
      
      pvalue_fun = label_style_pvalue(digits = 3)
    ) |>
    
    # Replace numeric liver-diagnosis codes
    modify_table_body(
      ~ .x |>
        mutate(
          label = case_when(
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "1" ~ "Hepatitis C",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "2" ~ "Hepatitis B",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "3" ~ "PSC/PBC/AIH",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "4" ~ "Alcohol-related",
            
            variable == "liver_diagnosis" &
              row_type == "level" &
              label == "5" ~ "Other",
            
            TRUE ~ label
          )
        )
    ) |>
    
    # Merge only rows that have an estimate
    modify_column_merge(
      pattern = "{estimate} ({conf.low}, {conf.high})",
      rows = !is.na(estimate)
    ) |>
    
    # Remove the p-value column
    modify_column_hide(
      columns = p.value
    ) |>
    
    modify_header(
      label = "**Predictor**",
      estimate = if (is_logistic) {
        "**aOR (95% CI)**"
      } else {
        "**B (95% CI)**"
      }
    )
  
  regression_table
}

# ------------------------------------------------------------------------------
# Predictors of sleep disturbance
# ------------------------------------------------------------------------------

# Logistic regression model table ----------------------------------------------

tbl_psqi_logistic_main <- make_main_regression_table(
  model = m_log_psqi,
  model_type = "logistic"
)

tbl_ess_logistic_main <- make_main_regression_table(
  model = m_log_ess,
  model_type = "logistic"
)

tbl_bss_logistic_main <- make_main_regression_table(
  model = m_log_bss,
  model_type = "logistic"
)

tbl_ais_logistic_main <- make_main_regression_table(
  model = m_log_ais,
  model_type = "logistic"
)

# Merge and format combined logistic regression model table

tbl_logistic_combined <- tbl_merge(
  tbls = list(
    tbl_psqi_logistic_main,
    tbl_ess_logistic_main,
    tbl_ais_logistic_main,
    tbl_bss_logistic_main
  ),
  
  tab_spanner = c(
    paste0(
      "**PSQI**  \nN = ",
      nobs(m_log_psqi)
    ),
    paste0(
      "**ESS**  \nN = ",
      nobs(m_log_ess)
    ),
    paste0(
      "**AIS**  \nN = ",
      nobs(m_log_ais)
    ),
    paste0(
      "**BSS**  \nN = ",
      nobs(m_log_bss)
    )
  )
) |>
  
  relabel_liver_diagnosis() |>
  
  # Combine each aOR and its CI
  modify_column_merge(
    pattern = "{estimate_1} ({conf.low_1}, {conf.high_1})",
    rows = !is.na(estimate_1)
  ) |>
  
  modify_column_merge(
    pattern = "{estimate_2} ({conf.low_2}, {conf.high_2})",
    rows = !is.na(estimate_2)
  ) |>
  
  modify_column_merge(
    pattern = "{estimate_3} ({conf.low_3}, {conf.high_3})",
    rows = !is.na(estimate_3)
  ) |>
  
  modify_column_merge(
    pattern = "{estimate_4} ({conf.low_4}, {conf.high_4})",
    rows = !is.na(estimate_4)
  ) |>
  
  # Remove all p-value columns
  modify_column_hide(
    columns = starts_with("p.value_")
  ) |>
  
  modify_header(
    label = "**Predictor**",
    estimate_1 = "**aOR (95% CI)**",
    estimate_2 = "**aOR (95% CI)**",
    estimate_3 = "**aOR (95% CI)**",
    estimate_4 = "**aOR (95% CI)**"
  )

tbl_logistic_combined

tbl_logistic_combined |>
  as_gt() |>
  opt_horizontal_padding(scale = 3) |>
  gtsave(
    filename = "tbl_log_models.png",
    zoom = 2,
    expand = 20
  )

# Linear regression model table ------------------------------------------------

tbl_psqi_linear_main <- make_main_regression_table(
  model = m_lin_psqi,
  model_type = "linear"
)

tbl_ess_linear_main <- make_main_regression_table(
  model = m_lin_ess,
  model_type = "linear"
)

tbl_ais_linear_main <- make_main_regression_table(
  model = m_lin_ais,
  model_type = "linear"
)

# Merge and format combine linear regression model table

# Combine the linear regression tables
tbl_linear_combined <- tbl_merge(
  tbls = list(
    tbl_psqi_linear_main,
    tbl_ess_linear_main,
    tbl_ais_linear_main
  ),
  tab_spanner = c(
    paste0(
      "**PSQI score**  \nN = ",
      nobs(m_lin_psqi)
    ),
    paste0(
      "**ESS score**  \nN = ",
      nobs(m_lin_ess)
    ),
    paste0(
      "**AIS score**  \nN = ",
      nobs(m_lin_ais)
    )
  )
) |>
  modify_header(
    label = "**Predictor**",
    estimate_1 = "**B (95% CI)**",
    estimate_2 = "**B (95% CI)**",
    estimate_3 = "**B (95% CI)**"
  )

tbl_linear_combined |>
  as_gt() |>
  opt_horizontal_padding(scale = 3) |>
  gtsave(
    filename = "tbl_linear_models.png",
    zoom = 2,
    expand = 20
  )

# ==============================================================================
# RQ2: Combined SF-36 PCS and MCS regression table
# ==============================================================================

qol_sleep_labels <- list(
  psqi_score = "PSQI score (per 1-point increase)",
  ess_score = "ESS score (per 1-point increase)",
  bss_score = "BSS high likelihood (vs. low likelihood)",
  ais_score = "AIS score (per 1-point increase)"
)

# Create the individual model tables
make_qol_model_table <- function(model) {
  model |>
    tbl_regression(
      include = all_of(sleep_measure_scores_all),
      intercept = FALSE,
      conf.level = 0.95,
      label = qol_sleep_labels,
      estimate_fun = label_style_number(digits = 2)
    )
}

tbl_pcs_sleep <- make_qol_model_table(
  m_lin_pcs
)

tbl_pcs_adjusted <- make_qol_model_table(
  m_lin_pcs_adj
)

tbl_mcs_sleep <- make_qol_model_table(
  m_lin_mcs
)

tbl_mcs_adjusted <- make_qol_model_table(
  m_lin_mcs_adj
)

# Combine the four model tables
tbl_qol_models <- tbl_merge(
  tbls = list(
    tbl_pcs_sleep,
    tbl_pcs_adjusted,
    tbl_mcs_sleep,
    tbl_mcs_adjusted
  ),
  tab_spanner = c(
    paste0(
      "**Sleep measures-only model**  \n",
      "N = ", nobs(m_lin_pcs),
      "; R² = ",
      sprintf("%.2f", summary(m_lin_pcs)$r.squared)
    ),
    paste0(
      "**Adjusted model**  \n",
      "N = ", nobs(m_lin_pcs_adj),
      "; R² = ",
      sprintf("%.2f", summary(m_lin_pcs_adj)$r.squared)
    ),
    paste0(
      "**Sleep-measure model**  \n",
      "N = ", nobs(m_lin_mcs),
      "; R² = ",
      sprintf("%.2f", summary(m_lin_mcs)$r.squared)
    ),
    paste0(
      "**Adjusted model**  \n",
      "N = ", nobs(m_lin_mcs_adj),
      "; R² = ",
      sprintf("%.2f", summary(m_lin_mcs_adj)$r.squared)
    )
  )
) |>
  
  # Combine coefficients and confidence intervals
  modify_column_merge(
    pattern = "{estimate_1} ({conf.low_1}, {conf.high_1})",
    rows = !is.na(estimate_1)
  ) |>
  modify_column_merge(
    pattern = "{estimate_2} ({conf.low_2}, {conf.high_2})",
    rows = !is.na(estimate_2)
  ) |>
  modify_column_merge(
    pattern = "{estimate_3} ({conf.low_3}, {conf.high_3})",
    rows = !is.na(estimate_3)
  ) |>
  modify_column_merge(
    pattern = "{estimate_4} ({conf.low_4}, {conf.high_4})",
    rows = !is.na(estimate_4)
  ) |>
  
  # Remove p-value columns
  modify_column_hide(
    columns = starts_with("p.value_")
  ) |>
  
  # Column headings
  modify_header(
    label = "**Sleep measure**",
    estimate_1 = "**B (95% CI)**",
    estimate_2 = "**B (95% CI)**",
    estimate_3 = "**B (95% CI)**",
    estimate_4 = "**B (95% CI)**"
  ) |>
  
  # Higher-level spanning headers
  modify_spanning_header(
    c(estimate_1, estimate_2) ~ "**PCS**",
    c(estimate_3, estimate_4) ~ "**MCS**",
    level = 2
  )

tbl_qol_models <- tbl_qol_models |>
  as_gt() |>
  tab_source_note(
    source_note = md(
      paste(
        "**Note.** B = unstandardized regression coefficient;",
        "CI = confidence interval.",
        "Sleep-measure models include PSQI, ESS, BSS, and AIS",
        "simultaneously. Adjusted models additionally include age,",
        "gender, BMI, liver diagnosis, recurrence of disease,",
        "rejection or graft dysfunction, fibrosis, renal failure,",
        "depression, corticosteroid use, and time from transplant."
      )
    )
  ) |>
  opt_horizontal_padding(scale = 1.5) |>
  tab_options(
    table.font.size = px(12)
  )

gtsave(
  data = tbl_qol_models_gt,
  filename = "tbl_sf36_qol_models.png",
  zoom = 2,
  expand = 20
)

# ==============================================================================
# Figure 2: sleep scores by significant clinical predictors (boxplots)
# ==============================================================================

# Prepare plotting data
figure2_data <- df_cleaned_lbl |>
  mutate(
    bss_risk = factor(
      bss_score,
      levels = c(0, 1),
      labels = c("Low risk", "High risk")
    )
  )

# Extract colours from the Epinephelus lanceolatus palette
epinephelus_palette <- paletteer_d(
  "fishualize::Epinephelus_lanceolatus"
) |>
  as.character()

# Use the pale grey-blue for No/Low and dark blue for Yes/High
binary_fill <- setNames(
  unname(epinephelus_palette[c(3, 1)]),
  c("No", "Yes")
)

bss_fill <- setNames(
  unname(epinephelus_palette[c(3, 1)]),
  c("Low risk", "High risk")
)

# Shared theme
box_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 12,
      margin = margin(b = 4)
    )
  )

# Reusable boxplot function
make_boxplot <- function(
    data,
    x,
    y,
    x_label,
    y_label,
    title
) {
  ggplot(
    data = data,
    mapping = aes(
      x = {{ x }},
      y = {{ y }},
      fill = {{ x }}
    )
  ) +
    geom_boxplot(
      width = 0.65,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = unname(
        as.character(
          paletteer::paletteer_d(
            "fishualize::Epinephelus_lanceolatus"
          )
        )[c(2, 3)]
      )
    ) +
    labs(
      x = x_label,
      y = y_label,
      title = title
    ) +
    box_theme
}

p1 <- make_boxplot(
  data = figure2_data |>
    filter(!is.na(bss_risk)),
  x = bss_risk,
  y = bmi,
  x_label = "Berlin Sleepiness Scale (BSS)",
  y_label = "BMI (kg/m²)",
  title = "A. BMI by Berlin sleep apnea risk"
)

p2 <- make_boxplot(
  data = figure2_data,
  x = depression,
  y = psqi_score,
  x_label = "Depression",
  y_label = "PSQI score",
  title = "B. PSQI score by depression status"
)

p3 <- make_boxplot(
  data = figure2_data,
  x = corticoid,
  y = ess_score,
  x_label = "Corticosteroid use",
  y_label = "ESS score",
  title = "C. ESS score by corticosteroid use"
)

p4 <- make_boxplot(
  data = figure2_data,
  x = recurrence_of_disease,
  y = ais_score,
  x_label = "Recurrence of disease",
  y_label = "AIS score",
  title = "D. AIS score by disease recurrence"
)

top_row <- p1 | p2
bottom_row <- p3 | p4

fig2_box <- top_row /
  plot_spacer() /
  bottom_row +
  plot_layout(
    heights = c(1, 0.02, 1)
  )

ggsave(
  filename = "figure2_clinical_predictors.png",
  plot = fig2_box,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# Figure 3: Forest plot of adjusted odds ratios
# ==============================================================================

# Terms included in the forest plot
keep_terms <- c(
  "age",
  "gender2",
  "bmi",
  "recurrence_of_disease1",
  "depression1",
  "corticoid1"
)

# Labels displayed in the plot
term_labels <- c(
  age = "Age",
  gender2 = "Gender: Female",
  bmi = "BMI",
  recurrence_of_disease1 = "Recurrence of disease",
  depression1 = "Depression",
  corticoid1 = "Corticosteroid use"
)

# Extract odds ratios and confidence intervals
get_or <- function(model, instrument) {
  model |>
    broom::tidy(
      exponentiate = TRUE,
      conf.int = TRUE
    ) |>
    filter(
      term %in% keep_terms
    ) |>
    mutate(
      instrument = instrument,
      term = unname(term_labels[term])
    )
}

# Combine results from the four logistic regression models
forest_data <- bind_rows(
  get_or(m_log_psqi, "PSQI"),
  get_or(m_log_ess, "ESS"),
  get_or(m_log_ais, "AIS"),
  get_or(m_log_bss, "BSS")
) |>
  mutate(
    term = factor(
      term,
      levels = rev(unname(term_labels))
    ),
    instrument = factor(
      instrument,
      levels = c(
        "PSQI",
        "ESS",
        "AIS",
        "BSS"
      )
    ),
    significant = p.value < 0.05
  )

# Create the forest plot
fig3_predictors <- ggplot(
  data = forest_data,
  mapping = aes(
    x = estimate,
    y = instrument
  )
) +
  
  # Reference line at an odds ratio of 1
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  # Confidence intervals
  geom_linerange(
    mapping = aes(
      xmin = conf.low,
      xmax = conf.high,
      colour = instrument
    ),
    linewidth = 0.8,
    show.legend = TRUE
  ) +
  
  # Odds-ratio estimates
  geom_point(
    mapping = aes(
      colour = instrument,
      size = significant
    ),
    shape = 16,
    show.legend = TRUE
  ) +
  
  # Make statistically significant estimates more noticeable
  scale_size_manual(
    values = c(
      "FALSE" = 2.5,
      "TRUE" = 5
    ),
    guide = "none"
  ) +
  
  # Original fishualize palette
  scale_colour_manual(
    values = {
      fish_colours <- as.character(
        paletteer_d(
          "fishualize::Epinephelus_lanceolatus"
        )
      )
      
      c(
        "PSQI" = unname(fish_colours[[1]]),
        "ESS" = unname(fish_colours[[2]]),
        "AIS" = "#6B7477",
        "BSS" = unname(fish_colours[[5]])
      )
    },
    name = NULL
  ) +
  
  scale_x_log10() +
  
  facet_grid(
    rows = vars(term),
    scales = "free_y",
    switch = "y"
  ) +
  
  labs(
    x = "Adjusted odds ratio (95% CI; log scale)",
    y = NULL,
    title = paste(
      "Figure 3. Odds ratios for key predictors of sleep",
      "disturbance across four instruments"
    ),
    caption = paste(
      "Note. Larger points indicate p < .05.",
      "Renal failure was excluded because of sparse data",
      "and complete separation."
    )
  ) +
  
  # Ensure the legend contains a line and circle
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 16,
        size = 3,
        linewidth = 0.8
      )
    )
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold"
    ),
    strip.placement = "outside",
    
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      margin = margin(t = 8)
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    legend.justification = "center",
    legend.box.just = "center",
    legend.direction = "vertical"
  )

# Save the forest plot
ggsave(
  filename = "fig3_predictors.png",
  plot = fig3_predictors,
  width = 8.5,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# Sensitivity analyses: complete-case versus multiple imputation
# ==============================================================================

# Labels for model terms --------------------------------------------------------

term_labels <- c(
  age = "Age",
  gender2 = "Gender: Female",
  bmi = "Body mass index (BMI)",
  liver_diagnosis2 = "Liver diagnosis: Hepatitis B",
  liver_diagnosis3 = "Liver diagnosis: PSC/PBC/AIH",
  liver_diagnosis4 = "Liver diagnosis: Alcohol-related",
  liver_diagnosis5 = "Liver diagnosis: Other",
  recurrence_of_disease1 = "Recurrence of disease: Yes",
  rejection_graft_dysfunction1 = "Rejection or graft dysfunction: Yes",
  any_fibrosis1 = "Any fibrosis: Yes",
  renal_failure1 = "Renal failure: Yes",
  depression1 = "Depression: Yes",
  corticoid1 = "Corticosteroid use: Yes",
  time_from_transplant = "Time from transplant (years)",
  psqi_score = "PSQI score",
  ess_score = "ESS score",
  bss_score = "BSS high likelihood",
  ais_score = "AIS score"
)

# Relabel model terms -----------------------------------------------------------

pretty_term <- function(term) {
  if_else(
    term %in% names(term_labels),
    unname(term_labels[term]),
    term
  )
}

# Format coefficients and confidence intervals ---------------------------------

format_estimate_ci <- function(
    estimate,
    conf_low,
    conf_high,
    digits = 2
) {
  paste0(
    formatC(
      estimate,
      format = "f",
      digits = digits
    ),
    " (",
    formatC(
      conf_low,
      format = "f",
      digits = digits
    ),
    ", ",
    formatC(
      conf_high,
      format = "f",
      digits = digits
    ),
    ")"
  )
}

# Format p-values using APA style -----------------------------------------------

format_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ NA_character_,
    p_value < 0.001 ~ "<.001",
    TRUE ~ sub(
      pattern = "^0",
      replacement = "",
      x = sprintf("%.3f", p_value)
    )
  )
}

# Compare complete-case and multiple-imputation results -------------------------

cc_mi_compare <- function(
    cc_model,
    pooled_model,
    digits = 2
) {
  cc_results <- cc_model |>
    tidy(
      conf.int = TRUE
    ) |>
    filter(
      term != "(Intercept)"
    ) |>
    transmute(
      predictor = pretty_term(term),
      cc_estimate = format_estimate_ci(
        estimate,
        conf.low,
        conf.high,
        digits
      ),
      cc_p = format_p_value(p.value)
    )
  
  mi_results <- pooled_model |>
    summary(
      conf.int = TRUE
    ) |>
    filter(
      term != "(Intercept)"
    ) |>
    transmute(
      predictor = pretty_term(term),
      mi_estimate = format_estimate_ci(
        estimate,
        `2.5 %`,
        `97.5 %`,
        digits
      ),
      mi_p = format_p_value(p.value)
    )
  
  cc_results |>
    full_join(
      mi_results,
      by = "predictor"
    )
}

# Create a formatted sensitivity table -----------------------------------------

make_sensitivity_table <- function(
    cc_model,
    pooled_model,
    table_title
) {
  cc_mi_compare(
    cc_model = cc_model,
    pooled_model = pooled_model
  ) |>
    gt() |>
    tab_header(
      title = md(
        paste0(
          "**",
          table_title,
          "**"
        )
      )
    ) |>
    tab_spanner(
      label = md("**Complete-case analysis**"),
      columns = c(
        cc_estimate,
        cc_p
      )
    ) |>
    tab_spanner(
      label = md("**Multiple imputation**"),
      columns = c(
        mi_estimate,
        mi_p
      )
    ) |>
    cols_label(
      predictor = md("**Predictor**"),
      cc_estimate = md("**B (95% CI)**"),
      cc_p = md("***p***"),
      mi_estimate = md("**B (95% CI)**"),
      mi_p = md("***p***")
    ) |>
    cols_align(
      align = "left",
      columns = predictor
    ) |>
    cols_align(
      align = "center",
      columns = c(
        cc_estimate,
        cc_p,
        mi_estimate,
        mi_p
      )
    ) |>
    tab_source_note(
      source_note = md(
        paste(
          "**Note.** B = unstandardized regression coefficient;",
          "CI = confidence interval."
        )
      )
    ) |>
    opt_horizontal_padding(
      scale = 1.5
    ) |>
    tab_options(
      table.font.size = px(12),
      data_row.padding = px(5)
    )
}

# Create sensitivity-analysis tables -------------------------------------------

tbl_sens_psqi <- make_sensitivity_table(
  cc_model = m_lin_psqi,
  pooled_model = pooled_psqi,
  table_title = paste(
    "Sensitivity analysis of predictors of PSQI score:",
    "complete-case analysis versus multiple imputation"
  )
)

tbl_sens_pcs <- make_sensitivity_table(
  cc_model = m_lin_pcs_adj,
  pooled_model = pooled_pcs,
  table_title = paste(
    "Sensitivity analysis of the adjusted SF-36 PCS model:",
    "complete-case analysis versus multiple imputation"
  )
)

tbl_sens_mcs <- make_sensitivity_table(
  cc_model = m_lin_mcs_adj,
  pooled_model = pooled_mcs,
  table_title = paste(
    "Sensitivity analysis of the adjusted SF-36 MCS model:",
    "complete-case analysis versus multiple imputation"
  )
)

# Display tables ---------------------------------------------------------------

tbl_sens_psqi
tbl_sens_pcs
tbl_sens_mcs

# Save tables ------------------------------------------------------------------

tbl_sens_psqi |>
  gtsave(
    filename = "tbl_sens_psqi.png",
    zoom = 2,
    expand = 20
  )

tbl_sens_pcs |>
  gtsave(
    filename = "tbl_sens_pcs.png",
    zoom = 2,
    expand = 20
  )

tbl_sens_mcs |>
  gtsave(
    filename = "tbl_sens_mcs.png",
    zoom = 2,
    expand = 20
  )

# ==============================================================================
# Heatmap for logistic regression models
# ==============================================================================
# ---- readable labels for each model term ----
clean_labels <- c(
  age = "Age",
  gender2 = "Gender: Female",
  bmi = "BMI",
  liver_diagnosis2 = "Liver diagnosis: Hep B",
  liver_diagnosis3 = "Liver diagnosis: PSC/PBC/AIH",
  liver_diagnosis4 = "Liver diagnosis: Alcohol",
  liver_diagnosis5 = "Liver diagnosis: Other",
  recurrence_of_disease1 = "Recurrence of disease",
  rejection_graft_dysfunction1 = "Rejection/graft dysfunction",
  any_fibrosis1 = "Fibrosis",
  renal_failure1 = "Renal failure",
  depression1 = "Depression",
  corticoid1 = "Corticosteroid use",
  time_from_transplant = "Time from transplant"
)

# ---- pull tidy coefficients from each stepAIC-reduced model ----
tidy_model <- function(model, instrument, exponentiate) {
  broom::tidy(model, exponentiate = exponentiate, conf.int = FALSE) |>
    filter(term != "(Intercept)") |>
    mutate(instrument = instrument,
           predictor = dplyr::coalesce(unname(clean_labels[term]), term))
}

binary_coefs <- bind_rows(
  tidy_model(m_log_psqi_step, "PSQI", exponentiate = TRUE),
  tidy_model(m_log_ess_step,  "ESS",  exponentiate = TRUE),
  tidy_model(m_log_ais_step,  "AIS",  exponentiate = TRUE),
  tidy_model(m_log_bss_step,  "BSS",  exponentiate = TRUE)
)

# 1. Drop renal failure — same exclusion rationale as your forest plot
binary_coefs <- binary_coefs |> filter(!str_starts(term, "renal_failure"))

# 2. Collapse liver diagnosis to one row (keep the most significant contrast)
liver_row <- binary_coefs |>
  filter(str_starts(term, "liver_diagnosis")) |>
  slice_min(p.value, n = 1, by = instrument) |>
  mutate(predictor = "Liver diagnosis*")

binary_coefs <- binary_coefs |>
  filter(!str_starts(term, "liver_diagnosis")) |>
  bind_rows(liver_row)

# 3. Rebuild grid, add log(OR) for a color scale that's actually legible
grid_binary <- expand_grid(
  predictor  = unique(binary_coefs$predictor),
  instrument = c("PSQI", "ESS", "AIS", "BSS")
) |>
  left_join(binary_coefs |> dplyr::select(predictor, instrument, estimate),
            by = c("predictor", "instrument")) |>
  mutate(log_estimate = log(estimate))

# order rows: most consistently-selected predictors on top
predictor_order <- grid_binary |>
  filter(!is.na(estimate)) |>
  count(predictor) |>
  arrange(n) |>
  pull(predictor)

grid_binary <- grid_binary |> mutate(predictor = factor(predictor, levels = predictor_order))

ggplot(grid_binary, aes(x = instrument, y = predictor, fill = log_estimate)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
            size = 4.2, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    na.value = "grey92", name = "Odds\nRatio",
    breaks = log(c(0.25, 0.5, 1, 2, 4)), labels = c("0.25", "0.5", "1", "2", "4")
  ) +
  labs(title = "AIC-Selected Predictors — Binary (Logistic) Models",
       subtitle = "Blank = not retained by AIC. Renal failure excluded (unstable due to separation, n=4)",
       caption = "*Liver diagnosis shows the contrast with the largest effect (ESS model only); other levels omitted for clarity.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold", size = 12))

ggsave("heatmap_binary_predictors.png", width = 7.5, height = 5.5, dpi = 300)

# ==============================================================================
# Heatmap for linear regression models
# ==============================================================================
cont_coefs <- bind_rows(
  tidy_model(m_lin_psqi_step, "PSQI", exponentiate = FALSE),
  tidy_model(m_lin_ess_step,  "ESS",  exponentiate = FALSE),
  tidy_model(m_lin_ais_step,  "AIS",  exponentiate = FALSE)
)

grid_cont <- expand_grid(
  predictor  = unique(cont_coefs$predictor),
  instrument = c("PSQI", "ESS", "AIS")
) |>
  left_join(cont_coefs |> dplyr::select(predictor, instrument, estimate),
            by = c("predictor", "instrument"))

ggplot(grid_cont, aes(x = instrument, y = fct_rev(predictor), fill = estimate)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = ifelse(is.na(estimate), "", sprintf("%.2f", estimate))),
            size = 4.2, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    na.value = "grey92", name = "β\ncoefficient"
  ) +
  labs(title = "AIC-Selected Predictors — Continuous (Linear) Models",
       subtitle = "Blank cells = not retained by stepwise AIC selection for that instrument",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold", size = 12))

ggsave("heatmap_continuous_predictors.png", width = 7.5, height = 5.5, dpi = 300)

# ==============================================================================
# SF-36 PCS and MCS coefficient plots
# ==============================================================================

# Sleep-measure labels ---------------------------------------------------------

qol_sleep_labels <- c(
  psqi_score = "PSQI score",
  ess_score = "ESS score",
  ais_score = "AIS score",
  bss_score = "BSS: high vs. low likelihood"
)

# Extract sleep-measure estimates from a model ---------------------------------

extract_qol_estimates <- function(
    model,
    model_name
) {
  model |>
    tidy(
      conf.int = TRUE
    ) |>
    filter(
      term %in% names(qol_sleep_labels)
    ) |>
    transmute(
      sleep_measure = unname(
        qol_sleep_labels[term]
      ),
      model = model_name,
      estimate,
      conf.low,
      conf.high
    )
}

# Prepare PCS results ----------------------------------------------------------

pcs_plot_data <- bind_rows(
  extract_qol_estimates(
    model = m_lin_pcs,
    model_name = "Sleep measures only"
  ),
  extract_qol_estimates(
    model = m_lin_pcs_adj,
    model_name = "Adjusted"
  )
)

# Prepare MCS results ----------------------------------------------------------

mcs_plot_data <- bind_rows(
  extract_qol_estimates(
    model = m_lin_mcs,
    model_name = "Sleep measures only"
  ),
  extract_qol_estimates(
    model = m_lin_mcs_adj,
    model_name = "Adjusted"
  )
)

# Apply consistent factor ordering ---------------------------------------------

pcs_plot_data <- pcs_plot_data |>
  mutate(
    sleep_measure = factor(
      sleep_measure,
      levels = rev(
        unname(qol_sleep_labels)
      )
    ),
    model = factor(
      model,
      levels = c(
        "Sleep measures only",
        "Adjusted"
      )
    )
  )

mcs_plot_data <- mcs_plot_data |>
  mutate(
    sleep_measure = factor(
      sleep_measure,
      levels = rev(
        unname(qol_sleep_labels)
      )
    ),
    model = factor(
      model,
      levels = c(
        "Sleep measures only",
        "Adjusted"
      )
    )
  )

# Use the same horizontal scale on both slides ---------------------------------

qol_limits <- range(
  c(
    pcs_plot_data$conf.low,
    pcs_plot_data$conf.high,
    mcs_plot_data$conf.low,
    mcs_plot_data$conf.high
  ),
  na.rm = TRUE
)

qol_padding <- diff(qol_limits) * 0.08

qol_limits <- qol_limits +
  c(
    -qol_padding,
    qol_padding
  )

# Reusable plotting function ---------------------------------------------------

make_qol_coefficient_plot <- function(
    plot_data,
    outcome,
    title,
    n_sleep,
    n_adjusted
) {
  dodge_position <- position_dodge(
    width = 0.5
  )
  
  ggplot(
    data = plot_data,
    mapping = aes(
      x = estimate,
      y = sleep_measure,
      colour = model,
      group = model
    )
  ) +
    
    # No-association reference line
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey35",
      linewidth = 0.8
    ) +
    
    # Confidence intervals
    geom_linerange(
      mapping = aes(
        xmin = conf.low,
        xmax = conf.high
      ),
      position = dodge_position,
      orientation = "y",
      linewidth = 1.1
    ) +
    
    # Coefficient estimates
    geom_point(
      position = dodge_position,
      shape = 16,
      size = 4
    ) +
    
    # Epinephelus_lanceolatus palette
    scale_colour_manual(
      values = c(
        "Sleep measures only" = unname(
          as.character(
            paletteer_d(
              "fishualize::Epinephelus_lanceolatus"
            )[[2]]
          )
        ),
        "Adjusted" = unname(
          as.character(
            paletteer_d(
              "fishualize::Epinephelus_lanceolatus"
            )[[1]]
          )
        )
      ),
      name = NULL
    ) +
    
    # Add major and minor vertical grid lines
    scale_x_continuous(
      breaks = qol_major_breaks,
      minor_breaks = qol_minor_breaks
    ) +
    
    coord_cartesian(
      xlim = qol_limits
    ) +
    
    labs(
      x = paste0(
        "Change in SF-36 ",
        outcome,
        " score, B (95% CI)"
      ),
      y = NULL,
      title = title
    ) +
    
    theme_minimal(
      base_size = 14
    ) +
    
    theme(
      panel.grid.major.x = element_line(
        colour = "grey78",
        linewidth = 0.5
      ),
      panel.grid.minor.x = element_line(
        colour = "grey88",
        linewidth = 0.35
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      plot.title.position = "plot",
      plot.caption.position = "plot",
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 18,
        lineheight = 1.1,
        margin = margin(
          b = 8
        )
      ),
      
      plot.subtitle = element_text(
        size = 13,
        lineheight = 1.1,
        margin = margin(
          b = 12
        )
      ),
      
      plot.caption = element_text(
        size = 10,
        hjust = 0,
        lineheight = 1.15,
        margin = margin(
          t = 12
        )
      ),
      
      axis.title.x = element_text(
        size = 13,
        margin = margin(
          t = 8
        )
      ),
      
      axis.text = element_text(
        size = 13
      ),
      
      legend.position = "bottom",
      legend.justification = "left",
      legend.text = element_text(
        size = 13
      ),
      
      plot.margin = margin(
        t = 20,
        r = 35,
        b = 20,
        l = 25
      )
    )
}

# PCS slide --------------------------------------------------------------------

figure_pcs_qol <- make_qol_coefficient_plot(
  plot_data = pcs_plot_data,
  outcome = "PCS",
  title = paste(
    "Sleep Disturbance and Physical Quality of Life"
  ),
  n_sleep = nobs(m_lin_pcs),
  n_adjusted = nobs(m_lin_pcs_adj)
)

figure_pcs_qol

# MCS slide --------------------------------------------------------------------

figure_mcs_qol <- make_qol_coefficient_plot(
  plot_data = mcs_plot_data,
  outcome = "MCS",
  title = paste(
    "Sleep Disturbance and Mental Quality of Life"
  ),
  n_sleep = nobs(m_lin_mcs),
  n_adjusted = nobs(m_lin_mcs_adj)
)

figure_mcs_qol

# Save plots -------------------------------------------------------------------

ggsave(
  filename = "figure_pcs_qol_slide.png",
  plot = figure_pcs_qol,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "figure_mcs_qol_slide.png",
  plot = figure_mcs_qol,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)