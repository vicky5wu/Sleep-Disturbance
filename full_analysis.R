# ==============================================================================
# BTC1859 Group Project - Sleep Disturbance and Quality of Life
# ==============================================================================

library(dplyr)
library(tidyr)
library(janitor)
library(gtsummary)
library(binom) # for Wilson 95% CIs
library(car)

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
    
    digits = all_continuous() ~ 2,
    
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
      psqi_score = "Pittsburgh Sleep Quality Index score",
      ess_score = "Epworth Sleepiness Scale score",
      ais_score = "Athens Insomnia Scale score",
      bss_score = "High likelihood of sleep-disordered breathing (Berlin Sleepiness Scale)"
    ),
    
    missing = "no"
  )

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
    
    digits = all_continuous() ~ 2,
    
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
      corticoid = "Corticosteroid use"
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
    
    digits = all_continuous() ~ 2,
    
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
      corticoid = "Corticosteroid use"
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
    
    digits = all_continuous() ~ 2,
    
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
      corticoid = "Corticosteroid use"
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
    
    digits = all_continuous() ~ 2,
    
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
      corticoid = "Corticosteroid use"
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

# RQ1a: Prevalence with proper 95% CIs (Wilson score interval) -----------------

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


# ============================================================
# Binary vs. continuous: which should we use downstream?
# ============================================================
# RQ1a itself ("what is the PREVALENCE") requires a binary case
# definition by construction -- prevalence is a %, so the table above
# stands regardless. The real choice is for RQ1b (predictors) and RQ2
# (relationship with QoL), where the outcome could be modeled either
# way. There's no single hypothesis test that declares one "correct" --
# it's a modeling decision -- but here's a principled way to handle it:
#
# (A) The general statistical argument (goes in your Discussion):
#     dichotomizing a continuous score at a cutoff throws away
#     information. Two patients on opposite sides of PSQI = 4 (say 3
#     vs. 5) get treated as maximally different, while PSQI = 5 and
#     PSQI = 25 get treated as identical. This is a well-documented way
#     to lose statistical power and distort effect estimates (see e.g.
#     Royston, Altman & Sauerbrei, "Dichotomizing continuous predictors
#     in multiple regression: a bad idea," Statistics in Medicine,
#     2006 -- confirm this citation yourselves before using it).
#
# (B) An empirical check on YOUR data: for a given outcome, fit the
#     same relationship once with the continuous score and once with
#     the binarized version, and compare fit. Example below uses SF36
#     PCS as the outcome and PSQI as the predictor:
#
# NOTE: previously this referenced df_updated (raw, pre-cleaning data
# frame) and PSQI_binary, which does not exist as a column anywhere.
# Fixed to use df_cleaned with the correct snake_case column names.

m_cont <- lm(sf36_pcs ~ psqi_score, data = df_cleaned)
m_bin  <- lm(sf36_pcs ~ psqi_binary, data = df_cleaned)

summary(m_cont)$r.squared
summary(m_bin)$r.squared
AIC(m_cont, m_bin)   # lower AIC = better fit, same data/outcome so comparable

# Univariate measure screen ----------------------------------------------------
# counts missing values per sleep instrument
df_cleaned |>
  summarise(across(
    all_of(sleep_measure_binary),
    ~ sum(is.na(.))
  ))

df_cleaned |> tbl_summary(
  include = sleep_measure_scores,
  statistic = list(all_continuous() ~ "{mean} ({sd})")
)

library(MASS)   # for stepAIC() -- loaded here (not at the top), because
# MASS::select() would otherwise mask the bare select() call used earlier
# in the data-validation section above

# Helper: build a locked-in complete-case subset before fitting -----------------
# stepAIC() errors out ("number of rows in use has changed") if dropping a
# predictor changes which rows have complete data -- which happens here
# because bmi (23 missing) and age (2 missing) are candidate predictors, and
# each sleep/QoL outcome has its own missingness on top. glm()/lm() silently
# drop incomplete rows per-model, so the full model and a stepwise-reduced
# model can end up fit on different samples, making their AIC incomparable.
# Fix: build the complete-case data ONCE per model (outcome + predictors in
# that model's formula only), then fit glm()/lm() directly against that named
# object -- NOT inside a wrapper function. (Fitting inside a helper function
# would break stepAIC()/update() afterwards: they re-evaluate the model's
# stored call in the environment attached to its formula, which is wherever
# psqi_formula etc. were originally created, i.e. this script's top level --
# not a function's local environment, which no longer exists once the
# function has returned. So the complete-case data frame needs to be a
# regular named object here, not a local variable inside a function.)

complete_case_data <- function(formula, data) {
  vars <- all.vars(formula)
  # dplyr::select() used explicitly here -- library(MASS) is now attached,
  # and MASS::select() would otherwise mask dplyr's, same issue as above
  data |> dplyr::select(all_of(vars)) |> drop_na()
}

# RQ1b: selecting predictors for four sleep disturbance instruments

# Multivariable logistic regression --------------------------------------------

# Variable selection using AIC-based stepwise selection (stepAIC)
# ================================================================
# stepwise selection is data-driven, and once a model is "chosen" by an automated
# procedure, p-values/CIs from that same model are no longer strictly
# valid (they're conditional on the selection having happened). The
# lecture explicitly states that choosing predictors from a-priori
# domain knowledge is preferred over pure stepwise selection. Our
# approach: keep both the full domain-knowledge model (already built
# above, using every clinically relevant covariate named in the
# assignment) and the stepAIC-reduced model, and compare/discuss them
# rather than silently picking one. Where they agree, that's your
# strongest evidence; where they disagree, that's worth a sentence in
# Results and Discussion.



# ---- PSQI ----------------------------------------------------------------

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
# "patients with X are __ times more likely to have PSQI-defined sleep disturbance")
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
     xlab = "Fitted values", ylab = "Pearson residuals"
)
abline(h = 0, lty = 2)

res_psqi_d <- residuals(m_log_psqi, type = "deviance")
plot(fitted(m_log_psqi), res_psqi_d,
     main = "PSQI model: Deviance residuals",
     xlab = "Fitted values", ylab = "Deviance residuals"
)
abline(h = 0, lty = 2)

# -- AIC-based selection (stepAIC), PSQI --
m_log_psqi_step <- stepAIC(m_log_psqi, direction = "backward", trace = FALSE)
summary(m_log_psqi_step)
# result model: psqi_binary ~ age + gender + recurrence_of_disease + 
# renal_failure + depression
exp(cbind(OR = coef(m_log_psqi_step), confint(m_log_psqi_step)))

# -- Continuous version: PSQI score as a continuous outcome, same predictors --
psqi_formula_cont <- reformulate(
  termlabels = c(demographic_variables, clinical_variables),
  response = "psqi_score"
)
df_psqi_cont_cc <- complete_case_data(psqi_formula_cont, df_cleaned)
m_lin_psqi <- lm(psqi_formula_cont, data = df_psqi_cont_cc)
summary(m_lin_psqi)
vif(m_lin_psqi)
par(mfrow = c(2, 2)); plot(m_lin_psqi); par(mfrow = c(1, 1))

m_lin_psqi_step <- stepAIC(m_lin_psqi, direction = "backward", trace = FALSE)
summary(m_lin_psqi_step)
# result model: psqi_score ~ age + gender + bmi + recurrence_of_disease + 
# rejection_graft_dysfunction + any_fibrosis + depression + corticoid

# -- Compare: did the same predictors survive selection on the continuous
# score vs. the binary threshold for PSQI? --
names(coef(m_log_psqi_step))
names(coef(m_lin_psqi_step))

# the shared predictors are age, gender, recurrence of disease, and depression

# ---- ESS -----------------------------------------------------------------

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
     xlab = "Fitted values", ylab = "Pearson residuals"
)
abline(h = 0, lty = 2)

res_ess_d <- residuals(m_log_ess, type = "deviance")
plot(fitted(m_log_ess), res_ess_d,
     main = "ESS model: Deviance residuals",
     xlab = "Fitted values", ylab = "Deviance residuals"
)
abline(h = 0, lty = 2)

# -- AIC-based selection (stepAIC), ESS --
m_log_ess_step <- stepAIC(m_log_ess, direction = "backward", trace = FALSE)
summary(m_log_ess_step)
# result model: ess_binary ~ liver_diagnosis + rejection_graft_dysfunction + 
# renal_failure + corticoid
exp(cbind(OR = coef(m_log_ess_step), confint(m_log_ess_step)))

# -- Continuous version: ESS score as a continuous outcome, same predictors --
ess_formula_cont <- reformulate(
  termlabels = c(demographic_variables, clinical_variables),
  response = "ess_score"
)
df_ess_cont_cc <- complete_case_data(ess_formula_cont, df_cleaned)
m_lin_ess <- lm(ess_formula_cont, data = df_ess_cont_cc)
summary(m_lin_ess)
vif(m_lin_ess)
par(mfrow = c(2, 2)); plot(m_lin_ess); par(mfrow = c(1, 1))

m_lin_ess_step <- stepAIC(m_lin_ess, direction = "backward", trace = FALSE)
summary(m_lin_ess_step)
# result model: ess_score ~ gender + rejection_graft_dysfunction + depression + corticoid

# -- Compare: continuous vs. binary predictor sets, ESS --
names(coef(m_log_ess_step))
names(coef(m_lin_ess_step))

# the shared predictors are rejection graft dysfunction and corticoid

# ---- BSS (binary only, no continuous score exists) -----------------------

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

# -- AIC-based selection (stepAIC), BSS (no continuous counterpart) --
m_log_bss_step <- stepAIC(m_log_bss, direction = "backward", trace = FALSE)
summary(m_log_bss_step)\
# result model: bss_score ~ bmi + recurrence_of_disease + renal_failure
exp(cbind(OR = coef(m_log_bss_step), confint(m_log_bss_step)))


# ---- AIS -----------------------------------------------------------------

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

# -- AIC-based selection (stepAIC), AIS --
m_log_ais_step <- stepAIC(m_log_ais, direction = "backward", trace = FALSE)
summary(m_log_ais_step)
# result model: is_binary ~ age + recurrence_of_disease + depression + corticoid
exp(cbind(OR = coef(m_log_ais_step), confint(m_log_ais_step)))

# -- Continuous version: AIS score as a continuous outcome, same predictors --
ais_formula_cont <- reformulate(
  termlabels = c(demographic_variables, clinical_variables),
  response = "ais_score"
)
df_ais_cont_cc <- complete_case_data(ais_formula_cont, df_cleaned)
m_lin_ais <- lm(ais_formula_cont, data = df_ais_cont_cc)
summary(m_lin_ais)
vif(m_lin_ais)
par(mfrow = c(2, 2)); plot(m_lin_ais); par(mfrow = c(1, 1))

m_lin_ais_step <- stepAIC(m_lin_ais, direction = "backward", trace = FALSE)
summary(m_lin_ais_step)
# result model: ais_score ~ age + bmi + recurrence_of_disease + corticoid

# -- Compare: continuous vs. binary predictor sets, AIS --
names(coef(m_log_ais_step))
names(coef(m_lin_ais_step))

# the shared predictors are age, reccurence of disease, and corticoid

# Which predictors survived AIC-based selection vs. the full model, per
# instrument? Divergence here is exactly what the assignment prompt
# anticipates ("some predictors better associated with some measures").
AIC(m_log_psqi, m_log_psqi_step)
AIC(m_log_ess,  m_log_ess_step)
AIC(m_log_bss,  m_log_bss_step)
AIC(m_log_ais,  m_log_ais_step)

# all the stepwise selected models have less predictors and smaller AIC values

# RQ2: Relationship of sleep disturbance with quality of life -----------------

# Multivariable linear regression ----------------------------------------------

# --- Crude models: sleep scores only, no adjustment for other patient traits --
# These answer "on their own, do sleep scores line up with QoL scores?" but
# cannot separate a real sleep effect from confounding (e.g. depression could
# independently affect both sleep and QoL).

m_lin_pcs <- lm(sf36_pcs ~ psqi_score + ess_score + bss_score + ais_score,
                data = df_cleaned
)
summary(m_lin_pcs)

m_lin_mcs <- lm(sf36_mcs ~ psqi_score + ess_score + bss_score + ais_score,
                data = df_cleaned
)
summary(m_lin_mcs)

# --- Adjusted models: sleep scores + demographic/clinical covariates ----------
# These ask "does sleep still matter once we account for age, gender, BMI,
# depression, etc.?" -- this is the more defensible version for drawing
# conclusions about sleep's relationship with QoL, since it controls for
# other patient characteristics that could confound the relationship.
#
# Clinical justification for including each covariate: each one is plausibly a 
# common cause of both sleep disturbance and poor QoL.
#   - depression: well-documented to independently worsen both sleep
#     (a core diagnostic symptom of depression) and QoL (especially
#     SF36 MCS). Without adjusting for it, some of what looks like a
#     "sleep effect" on QoL could really be a depression effect acting
#     on both.
#   - corticoid (corticosteroid usage): a drug of anti-rejection
#     regimens in transplant patients, and a well-known pharmacological
#     cause of insomnia/sleep disruption as a side effect, while also
#     independently affecting QoL through other side effects (mood
#     changes, weight gain, edema). A double-acting confounder specific
#     to this population.
#   - renal_failure, any_fibrosis, rejection_graft_dysfunction,
#     recurrence_of_disease: all markers of graft/disease severity.
#     Sicker patients could have both worse sleep (symptom burden,
#     more medications, more hospital visits) and worse physical QoL
#     (SF36 PCS) independent of any true causal sleep-QoL link.
#   - time_from_transplant: a proxy for recovery trajectory. Patients
#     early post-transplant may have both disrupted sleep (acute
#     recovery/hospitalization effects) and lower QoL that naturally
#     improves with time, so omitting it risks confounding a
#     cross-sectional sleep-QoL association with a "time since surgery"
#     effect.
#   - age, gender, bmi: standard demographic confounders. These act as baseline
#     QoL norms and sleep patterns both vary by these independent of
#     any transplant-specific factors.
#   - liver_diagnosis: underlying disease etiology may correlate with
#     different comorbidity/symptom profiles that could affect both
#     sleep and QoL differently.

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


# -- AIC-based selection (stepAIC), adjusted PCS model --
m_lin_pcs_adj_step <- stepAIC(m_lin_pcs_adj, direction = "backward", trace = FALSE)
summary(m_lin_pcs_adj_step)

# result model: sf36_pcs ~ ess_score + ais_score + age + bmi + recurrence_of_disease

# -- Residual diagnostics, adjusted PCS model --
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


# -- AIC-based selection (stepAIC), adjusted MCS model --
m_lin_mcs_adj_step <- stepAIC(m_lin_mcs_adj, direction = "backward", trace = FALSE)
summary(m_lin_mcs_adj_step)

# result model: sf36_mcs ~ psqi_score + ais_score + age + liver_diagnosis + 
# rejection_graft_dysfunction + any_fibrosis + depression + time_from_transplant

# -- Residual diagnostics, adjusted MCS model --
par(mfrow = c(2, 2)); plot(m_lin_mcs_adj); par(mfrow = c(1, 1))

AIC(m_lin_pcs_adj, m_lin_pcs_adj_step)
AIC(m_lin_mcs_adj, m_lin_mcs_adj_step)

# the backward stepwise AIC models are better, only ais_score survive in both models

# Do the sleep-score terms (psqi_score, ess_score, bss_score, ais_score)
# survive together in the stepAIC-reduced QoL models? If AIC drops some
# of them, that itself is informative -- it suggests those instruments
# don't add explanatory power once the others are accounted for. Worth
# checking m_pcs_adj_vif / m_mcs_adj_vif (computed above) specifically
# for the four sleep-score rows -- if any are high, the four instruments
# may be too collinear with each other to cleanly separate in one joint
# model, which would be an argument for also running each sleep
# instrument in its own separate QoL model as a complement to this one.


# Compare crude vs. adjusted fit -- if R^2 improves a lot and/or sleep-score
# coefficients change substantially once covariates are added, that's a sign
# the crude (sleep-only) models were confounded and the adjusted ones should
# be the ones you report as primary.
summary(m_lin_pcs)$r.squared #0.207
summary(m_lin_pcs_adj)$r.squared #0.292
summary(m_lin_pcs_adj_step)$r.squared #0.245
summary(m_lin_mcs)$r.squared #0.319
summary(m_lin_mcs_adj)$r.squared #0.474
summary(m_lin_mcs_adj_step)$r.squared #0.470


# ================================================================
# Sensitivity analysis: multiple imputation vs. complete-case, for PSQI
# ================================================================
# PSQI is missing in 85/268 patients (32%) -- by far the most-missing
# variable used anywhere in this analysis (compare: age 2 missing, bmi 23
# missing, ess 17 missing, ais 6 missing). Every PSQI model above
# (m_log_psqi, m_lin_psqi, and psqi_score's role inside m_lin_pcs_adj /
# m_lin_mcs_adj) is fit on a complete-case subset that silently drops
# those 85 patients. This section checks whether that matters: does using
# multiple imputation instead of complete-case meaningfully change the conclusions
# involving PSQI?
#
# We impute psqi_score together with the other variables it's analyzed
# alongside (predictors + other sleep scores + QoL outcomes), refit the
# same formulas already used above on each imputed dataset, and pool the
# results with Rubin's rules (mice::pool()).
#
# NOTE: this covers the continuous-score models only (psqi_score as an
# outcome, and as a predictor of QoL). Extending this to the binary
# m_log_psqi model would require re-deriving psqi_binary from each
# imputed psqi_score (passive imputation) rather than imputing the binary
# indicator directly, since imputing a threshold variable on its own can
# produce a value inconsistent with the underlying score. Flagged here as
# a possible extension rather than built out, to keep this section
# focused.

library(mice)

imp_vars <- c(
  demographic_variables, clinical_variables,
  "psqi_score", "ess_score", "ais_score", "bss_score",
  "sf36_pcs", "sf36_mcs"
)

df_for_imputation <- df_cleaned |> dplyr::select(all_of(imp_vars))

# With 32% missingness on psqi_score we're imputing a much larger share
# of the data than that toy example, so more imputations are warranted
# for stable pooled estimates. method = "pmm" (predictive mean matching,
# mice's default for numeric variables) draws imputed values from
# actually-observed donors, so imputed PSQI scores stay within the valid
# 0-21 range -- unlike norm.nob/norm.predict (shown in Tutorial 10 on a
# toy example), which don't respect that boundary.
imp <- mice(df_for_imputation, m = 20, method = "pmm", seed = 431859, print = FALSE)

## ---- PSQI as outcome (RQ1b, continuous) ----
fits_psqi_mi <- lapply(seq_len(imp$m), function(i) {
  lm(psqi_formula_cont, data = complete(imp, i))
})
pooled_psqi <- pool(as.mira(fits_psqi_mi))
summary(pooled_psqi)

# Compare against the complete-case version:
summary(m_lin_psqi)$coefficients
summary(pooled_psqi)

## ---- PSQI as a predictor of QoL (RQ2, adjusted models) ----
fits_pcs_mi <- lapply(seq_len(imp$m), function(i) {
  lm(pcs_formula_adj, data = complete(imp, i))
})
pooled_pcs <- pool(as.mira(fits_pcs_mi))
summary(pooled_pcs)

fits_mcs_mi <- lapply(seq_len(imp$m), function(i) {
  lm(mcs_formula_adj, data = complete(imp, i))
})
pooled_mcs <- pool(as.mira(fits_mcs_mi))
summary(pooled_mcs)

# Compare against the complete-case versions:
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

library(dplyr)
library(gtsummary)
library(gt)
library(broom)

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

tbl_psqi_log
tbl_psqi_lin
tbl_ess_log
tbl_ess_lin
tbl_bss_log
tbl_ais_log
tbl_ais_lin

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

tbl_pcs
tbl_mcs

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

tbl_sens_psqi
tbl_sens_pcs
tbl_sens_mcs
