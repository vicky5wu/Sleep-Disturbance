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
sum(is.na(df))
colSums(is.na(df))

# Subset data ------------------------------------------------------------------

# Keeping important variables (added Age it's required by the assignment as a 
# demographic variable, but wasn't in the original select())

df_updated <- df |>
  select(
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
  "corticoid"
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

df_cleaned |>
  summarise(across(
    all_of(sleep_measure_binary),
    ~ sum(is.na(.))
  ))

df_cleaned |> tbl_summary(
  include = sleep_measure_scores,
  statistic = list(all_continuous() ~ "{mean} ({sd})")
)

# Multivariable logistic regression --------------------------------------------

# Create logistic regression model for PSQI
psqi_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "psqi_binary"
)

m_log_psqi <- glm(psqi_formula,
  df_cleaned,
  family = "binomial"
)

m_psqi_vif <- vif(m_log_psqi)

res_psqi_p <- residuals(m_log_psqi, type = "pearson")
plot(fitted(m_log_psqi), res_psqi_p)
res_psqi_d <- residuals(m_log_psqi, type = "deviance")
plot(fitted(m_log_psqi), res_psqi_d)

# Create logistic regression model for ESS
ess_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "ess_binary"
)

m_log_ess <- glm(ess_formula,
  df_cleaned,
  family = "binomial"
)

m_ess_vif <- vif(m_log_ess)

# Create logistic regression model for BSS
bss_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "bss_score"
)

m_log_bss <- glm(bss_formula,
  df_cleaned,
  family = "binomial"
)

m_bss_vif <- vif(m_log_bss)

# Create logistic regression model for AIS
ais_formula <- reformulate(
  termlabels = c(
    demographic_variables,
    clinical_variables
  ),
  response = "ais_binary"
)

m_log_ais <- glm(ais_formula,
  df_cleaned,
  family = "binomial"
)

m_ais_vif <- vif(m_log_ais)

# Multivariable linear regression ----------------------------------------------

m_lin_pcs <- lm(sf36_pcs ~ psqi_score + ess_score + bss_score + ais_score,
  data = df_cleaned
)

m_lin_mcs <- lm(sf36_mcs ~ psqi_score + ess_score + bss_score + ais_score,
  data = df_cleaned
)

