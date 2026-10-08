# ============================================================
# 05. Actuarial application and risk measures: medical cost data
#
# Data: Medical Cost Personal Datasets, individual medical costs billed by
# health insurance for 1,338 policyholders in US dollars
# (data/medical_cost_charges.csv, column `charges`). All values are
# observed, so the data are complete.
#
# Steps:
#   1. descriptive statistics;
#   2. MLE of NME-Weibull with fitnmeweibull() and of Weibull with the same
#      method;
#   3. comparison by log-likelihood, AIC, and BIC;
#   4. fitted versus empirical density and survival (exceedance) function;
#   5. VaR and TVaR at 0.90, 0.95, 0.99 and RVaR on (0.90, 0.95) and
#      (0.95, 0.99), compared with the Weibull model and empirical values.
#
# fitnmeweibull() rescales the data internally, so the costs are fitted
# directly in dollars and all risk measures are in dollars.
#
# Output: results/05_risk_measures/, figures/05_risk_measures/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "05_risk_measures")
fig_dir <- output_dir("figures", "05_risk_measures")
unit_label <- "Biaya medis (USD)"

# ------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------
x <- utils::read.csv("data/medical_cost_charges.csv")$charges

if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0)) {
  stop("data/medical_cost_charges.csv must contain positive costs in a ",
       "column named charges.", call. = FALSE)
}

save_table(describe_data(x), file.path(result_dir, "descriptive_statistics.csv"))

# ------------------------------------------------------------
# 2-3. Estimation and model comparison
# ------------------------------------------------------------
fit_nme <- fitnmeweibull(x)
fit_weibull <- fit_weibull_mle(x)
print(fit_nme)

comparison <- compare_models(fit_nme, fit_weibull)
save_table(comparison, file.path(result_dir, "model_comparison.csv"))

# ------------------------------------------------------------
# 4. Figures
# ------------------------------------------------------------
plot_fitted_models(x, fit_nme, fit_weibull, fig_dir, unit_label)

# ------------------------------------------------------------
# 5. Risk measures (USD)
# ------------------------------------------------------------
risk <- risk_measures(x, fit_nme, fit_weibull)

# Relative difference of each model from the empirical value, in percent.
risk$var_tvar$TVaR_nme_vs_empirical_pct <-
  100 * (risk$var_tvar$TVaR_nme / risk$var_tvar$TVaR_empirical - 1)
risk$var_tvar$TVaR_weibull_vs_empirical_pct <-
  100 * (risk$var_tvar$TVaR_weibull / risk$var_tvar$TVaR_empirical - 1)
risk$rvar$RVaR_nme_vs_empirical_pct <-
  100 * (risk$rvar$RVaR_nme / risk$rvar$RVaR_empirical - 1)
risk$rvar$RVaR_weibull_vs_empirical_pct <-
  100 * (risk$rvar$RVaR_weibull / risk$rvar$RVaR_empirical - 1)

save_table(risk$var_tvar, file.path(result_dir, "var_tvar.csv"))
save_table(risk$rvar, file.path(result_dir, "rvar.csv"))

# Probability that a claim exceeds selected cost levels.
levels_usd <- c(10000, 20000, 30000, 40000, 50000)
exceedance <- data.frame(
  cost_usd = levels_usd,
  empirical = vapply(levels_usd, function(v) mean(x > v), numeric(1)),
  nme_weibull = snmeweibull(levels_usd, fit_nme$estimate["alpha"],
                            fit_nme$estimate["beta"],
                            fit_nme$estimate["lambda"]),
  weibull = sweibull_rate(levels_usd, fit_weibull$estimate["beta"],
                          fit_weibull$estimate["lambda"])
)
save_table(exceedance, file.path(result_dir, "exceedance_probability.csv"))

cat("\nResults written to", result_dir, "and", fig_dir, "\n")
