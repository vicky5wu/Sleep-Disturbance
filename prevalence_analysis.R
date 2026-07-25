library(dplyr)
library(tidyr)
library(binom)   # for Wilson 95% CIs

df <- read.csv("project_data.csv")

str(df)
dim(df)
sum(is.na(df))
colSums(is.na(df))


## ------------------------------------------------------------
## Keeping important variables
## (added Age -- it's required by the assignment as a demographic
##  variable, but wasn't in the original select())
## ------------------------------------------------------------
df_updated <- df %>%
  select(Subject, Age, Gender, BMI, Time.from.transplant,
         Liver.Diagnosis, Recurrence.of.disease, Rejection.graft.dysfunction,
         Any.fibrosis, Renal.Failure, Depression, Corticoid,

         Pittsburgh.Sleep.Quality.Index.Score, Epworth.Sleepiness.Scale,
         Berlin.Sleepiness.Scale, Athens.Insomnia.Scale,

         SF36.PCS, SF36.MCS)


## ------------------------------------------------------------
## Creating categorical variables for sleep disturbance instruments
## ------------------------------------------------------------
df_updated$ESS_binary  <- ifelse(df$Epworth.Sleepiness.Scale > 10, 1, 0)
df_updated$PSQI_binary <- ifelse(df$Pittsburgh.Sleep.Quality.Index.Score > 4, 1, 0)
df_updated$AIS_binary  <- ifelse(df$Athens.Insomnia.Scale > 5, 1, 0)

## ------------------------------------------------------------
## RQ1a: Prevalence with proper 95% CIs (Wilson score interval)
## ------------------------------------------------------------
prevalence_ci <- function(x_vec, label) {
  x <- sum(x_vec, na.rm = TRUE)
  n <- sum(!is.na(x_vec))
  ci <- binom.confint(x, n, methods = "wilson")
  data.frame(instrument = label, positive = x, n = n,
             prevalence_percentage = round(100 * ci$mean, 1),
             ci_lower = round(100 * ci$lower, 1),
             ci_upper = round(100 * ci$upper, 1))
}

prevalence_table <- bind_rows(
  prevalence_ci(df_updated$ESS_binary,             "ESS > 10 (daytime sleepiness)"),
  prevalence_ci(df_updated$PSQI_binary,             "PSQI > 4 (poor sleep quality)"),
  prevalence_ci(df_updated$AIS_binary,              "AIS > 5 (insomnia)"),
  prevalence_ci(df_updated$Berlin.Sleepiness.Scale, "Berlin (high SDB risk)")
)
print(prevalence_table)
write.csv(prevalence_table, "prevalence_estimates.csv", row.names = FALSE)


## ============================================================
## Binary vs. continuous: which should we use downstream?
## ============================================================
## RQ1a itself ("what is the PREVALENCE") requires a binary case
## definition by construction -- prevalence is a %, so the table above
## stands regardless. The real choice is for RQ1b (predictors) and RQ2
## (relationship with QoL), where the outcome could be modeled either
## way. There's no single hypothesis test that declares one "correct" --
## it's a modeling decision -- but here's a principled way to handle it:
##
## (A) The general statistical argument (goes in your Discussion):
##     dichotomizing a continuous score at a cutoff throws away
##     information. Two patients on opposite sides of PSQI = 4 (say 3
##     vs. 5) get treated as maximally different, while PSQI = 5 and
##     PSQI = 25 get treated as identical. This is a well-documented way
##     to lose statistical power and distort effect estimates (see e.g.
##     Royston, Altman & Sauerbrei, "Dichotomizing continuous predictors
##     in multiple regression: a bad idea," Statistics in Medicine,
##     2006 -- confirm this citation yourselves before using it).
##
## (B) An empirical check on YOUR data: for a given outcome, fit the
##     same relationship once with the continuous score and once with
##     the binarized version, and compare fit. Example below uses SF36
##     PCS as the outcome and PSQI as the predictor:

m_cont <- lm(SF36.PCS ~ Pittsburgh.Sleep.Quality.Index.Score, data = df_updated)
m_bin  <- lm(SF36.PCS ~ PSQI_binary, data = df_updated)

summary(m_cont)$r.squared
summary(m_bin)$r.squared
AIC(m_cont, m_bin)   # lower AIC = better fit, same data/outcome so comparable