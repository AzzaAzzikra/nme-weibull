# ============================================================
# 05. Supplementary applications (not part of the design in Bab III)
#
# The same comparison as script 04 on two further data sets in data/:
#   - medical_cost_charges.csv: individual medical costs billed by health
#     insurance (USD), Medical Cost Personal Datasets (Kaggle). Complete
#     data; relevant for the actuarial risk measures.
#   - gbsg2_recurrence_free_time.csv: recurrence-free survival time (days)
#     of the German Breast Cancer Study Group 2 trial. 56% of the patients
#     are censored, so only the observed recurrences (status == 1) are
#     used, as for the veteran data.
#
# Output: results/05_additional_data/, figures/05_additional_data/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "05_additional_data")
fig_dir <- output_dir("figures", "05_additional_data")

charges <- utils::read.csv("data/medical_cost_charges.csv")$charges
gbsg2 <- utils::read.csv("data/gbsg2_recurrence_free_time.csv")

data_sets <- list(
  medical_cost = list(x = charges, unit = "Biaya medis (USD)"),
  gbsg2_recurrence = list(x = gbsg2$time[gbsg2$status == 1],
                          unit = "Waktu hingga kekambuhan (hari)")
)

summary_rows <- list()

for (name in names(data_sets)) {
  x <- data_sets[[name]]$x
  x <- x[is.finite(x) & x > 0]

  fit_nme <- fitnmeweibull(x)
  fit_weibull <- fit_weibull_mle(x)
  comparison <- compare_models(fit_nme, fit_weibull)
  risk <- risk_measures(x, fit_nme, fit_weibull)

  save_table(describe_data(x),
             file.path(result_dir, paste0(name, "_descriptive.csv")))
  save_table(comparison,
             file.path(result_dir, paste0(name, "_model_comparison.csv")))
  save_table(risk$var_tvar, file.path(result_dir, paste0(name, "_var_tvar.csv")))
  save_table(risk$rvar, file.path(result_dir, paste0(name, "_rvar.csv")))
  plot_fitted_models(x, fit_nme, fit_weibull, fig_dir,
                     unit_label = data_sets[[name]]$unit,
                     prefix = paste0(name, "_"))

  summary_rows[[name]] <- data.frame(
    data = name,
    n = length(x),
    alpha = unname(fit_nme$estimate["alpha"]),
    alpha_at_bound = fit_nme$boundary,
    AIC_nme = fit_nme$AIC,
    AIC_weibull = fit_weibull$AIC,
    BIC_nme = fit_nme$BIC,
    BIC_weibull = fit_weibull$BIC,
    AIC_weibull_minus_nme = fit_weibull$AIC - fit_nme$AIC,
    BIC_weibull_minus_nme = fit_weibull$BIC - fit_nme$BIC
  )
}

save_table(do.call(rbind, summary_rows),
           file.path(result_dir, "summary.csv"))

cat("\nResults written to", result_dir, "and", fig_dir, "\n")
