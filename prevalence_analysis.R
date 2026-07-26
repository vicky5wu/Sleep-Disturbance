library(dplyr)
library(tidyr)
library(janitor)
library(binom) # for Wilson 95% CIs

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
  "ais_binary"
)

clinical_variables <- c(
  "liver_diagnosis",
  "reccurence_of_disease",
  "rejection_graft_dysfunction",
  "any_fibrosis",
  "renal_failure",
  "depression",
  "corticoid"
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
  prevalence_ci(df_updated$ESS_binary, "ESS > 10 (daytime sleepiness)"),
  prevalence_ci(df_updated$PSQI_binary, "PSQI > 4 (poor sleep quality)"),
  prevalence_ci(df_updated$AIS_binary, "AIS > 5 (insomnia)"),
  prevalence_ci(df_updated$Berlin.Sleepiness.Scale, "Berlin (high SDB risk)")
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

m_cont <- lm(SF36.PCS ~ Pittsburgh.Sleep.Quality.Index.Score, data = df_updated)
m_bin  <- lm(SF36.PCS ~ PSQI_binary, data = df_updated)

summary(m_cont)$r.squared
summary(m_bin)$r.squared
AIC(m_cont, m_bin)   # lower AIC = better fit, same data/outcome so comparable