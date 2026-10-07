# ============================================================
# 04. Application to real survival-time data (Subbab 3.7)
#
# Data: Veterans' Administration lung cancer trial (Kalbfleisch and
# Prentice, 2002), object `veteran` in the survival package. Only the 128
# patients who died during follow-up are used, because fitnmeweibull()
# handles complete (uncensored) data. data/real_data.csv holds these 128
# survival times in days (column `time`).
#
# Steps (Subbab 3.7):
#   1. descriptive statistics;
#   2. MLE of NME-Weibull with fitnmeweibull() and of Weibull with the same
#      method;
#   3. comparison by log-likelihood, AIC, and BIC, with AIC and BIC
#      differences;
#   4. fitted versus empirical survival functions;
#   5. VaR and TVaR at 0.90, 0.95, 0.99 and RVaR on (0.90, 0.95) and
#      (0.95, 0.99).
#
# Output: results/04_real_data/, figures/04_real_data/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "04_real_data")
fig_dir <- output_dir("figures", "04_real_data")

# ------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------
real_data <- utils::read.csv("data/real_data.csv")
x <- real_data$time

if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0)) {
  stop("data/real_data.csv must contain positive survival times in a ",
       "column named time.", call. = FALSE)
}

# Cross-check against the survival package when it is installed.
if (requireNamespace("survival", quietly = TRUE)) {
  veteran <- survival::veteran
  deaths <- veteran$time[veteran$status == 1]
  cat("Matches survival::veteran (status == 1):",
      identical(sort(as.numeric(deaths)), sort(as.numeric(x))), "\n")
  data_source <- data.frame(
    patients = nrow(veteran),
    deaths = sum(veteran$status == 1),
    censored = sum(veteran$status == 0),
    used = length(x)
  )
  save_table(data_source, file.path(result_dir, "data_source.csv"))
}

descriptive <- describe_data(x)
save_table(descriptive, file.path(result_dir, "descriptive_statistics.csv"))

# ------------------------------------------------------------
# 2-3. Estimation and model comparison
# ------------------------------------------------------------
fit_nme <- fitnmeweibull(x)
fit_weibull <- fit_weibull_mle(x)
print(fit_nme)

# The Weibull fit is checked against survreg(), which uses a different
# algorithm; with complete data both maximise the same likelihood.
if (requireNamespace("survival", quietly = TRUE)) {
  reference <- survival::survreg(survival::Surv(x, rep(1, length(x))) ~ 1,
                                 dist = "weibull")
  cat("Weibull logLik, fit_weibull_mle():", fit_weibull$logLik,
      " survreg():", as.numeric(stats::logLik(reference)), "\n")
}

comparison <- compare_models(fit_nme, fit_weibull)
save_table(comparison, file.path(result_dir, "model_comparison.csv"))

# ------------------------------------------------------------
# 4. Figures: fitted PDF, survival, and hazard
# ------------------------------------------------------------
plot_fitted_models(x, fit_nme, fit_weibull, fig_dir,
                   unit_label = "Waktu hidup (hari)")

# Survival probabilities at selected days, fitted versus empirical.
days <- c(30, 90, 180, 365)
survival_table <- data.frame(
  day = days,
  empirical = vapply(days, function(d) mean(x > d), numeric(1)),
  nme_weibull = snmeweibull(days, fit_nme$estimate["alpha"],
                            fit_nme$estimate["beta"],
                            fit_nme$estimate["lambda"]),
  weibull = sweibull_rate(days, fit_weibull$estimate["beta"],
                          fit_weibull$estimate["lambda"])
)
save_table(survival_table, file.path(result_dir, "survival_at_days.csv"))

# ------------------------------------------------------------
# 5. Risk measures (in days)
# ------------------------------------------------------------
risk <- risk_measures(x, fit_nme, fit_weibull)
save_table(risk$var_tvar, file.path(result_dir, "var_tvar.csv"))
save_table(risk$rvar, file.path(result_dir, "rvar.csv"))

cat("\nResults written to", result_dir, "and", fig_dir, "\n")
