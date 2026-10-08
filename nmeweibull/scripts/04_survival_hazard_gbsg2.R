# ============================================================
# 04. Survival and hazard modelling: GBSG2 data
#
# Data: German Breast Cancer Study Group 2 trial, 686 patients with
# node-positive breast cancer. data/gbsg2_recurrence_free_time.csv holds the
# recurrence-free survival time in days (`time`) and the status (`status`):
# 1 = recurrence or death, 0 = censored.
#
# fitnmeweibull() handles complete (uncensored) data, so only the 299
# patients with an observed event are used. The results therefore describe
# the time to recurrence or death among patients who experienced the event
# during follow-up.
#
# Steps:
#   1. descriptive statistics;
#   2. MLE of NME-Weibull with fitnmeweibull() and of Weibull with the same
#      method;
#   3. comparison by log-likelihood, AIC, and BIC;
#   4. fitted versus empirical survival function, Nelson-Aalen cumulative
#      hazard, and the fitted hazard functions.
#
# Output: results/04_survival_hazard/, figures/04_survival_hazard/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "04_survival_hazard")
fig_dir <- output_dir("figures", "04_survival_hazard")
unit_label <- "Waktu hingga kekambuhan atau kematian (hari)"

# ------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------
gbsg2 <- utils::read.csv("data/gbsg2_recurrence_free_time.csv")

if (!all(c("time", "status") %in% names(gbsg2)) ||
    any(!is.finite(gbsg2$time)) || any(gbsg2$time <= 0) ||
    any(!gbsg2$status %in% c(0, 1))) {
  stop("data/gbsg2_recurrence_free_time.csv must contain positive times in ",
       "column time and 0/1 values in column status.", call. = FALSE)
}

x <- gbsg2$time[gbsg2$status == 1]

data_source <- data.frame(
  patients = nrow(gbsg2),
  events = sum(gbsg2$status == 1),
  censored = sum(gbsg2$status == 0),
  censored_percent = 100 * mean(gbsg2$status == 0),
  used = length(x)
)
save_table(data_source, file.path(result_dir, "data_source.csv"))
save_table(describe_data(x), file.path(result_dir, "descriptive_statistics.csv"))

# ------------------------------------------------------------
# 2-3. Estimation and model comparison
# ------------------------------------------------------------
fit_nme <- fitnmeweibull(x)
fit_weibull <- fit_weibull_mle(x)
print(fit_nme)

# Check of the Weibull fit against survreg(), a different algorithm that
# maximises the same likelihood for complete data.
if (requireNamespace("survival", quietly = TRUE)) {
  reference <- survival::survreg(survival::Surv(x, rep(1, length(x))) ~ 1,
                                 dist = "weibull")
  cat("Weibull logLik, fit_weibull_mle():", fit_weibull$logLik,
      " survreg():", as.numeric(stats::logLik(reference)), "\n")
}

comparison <- compare_models(fit_nme, fit_weibull)
save_table(comparison, file.path(result_dir, "model_comparison.csv"))

# ------------------------------------------------------------
# 4. Survival, cumulative hazard, and hazard
# ------------------------------------------------------------
plot_fitted_models(x, fit_nme, fit_weibull, fig_dir, unit_label,
                   hazard_breaks = seq(0, 2557.5, by = 182.5))
plot_cumulative_hazard(x, fit_nme, fit_weibull, fig_dir, unit_label)

a <- unname(fit_nme$estimate["alpha"])
b <- unname(fit_nme$estimate["beta"])
l <- unname(fit_nme$estimate["lambda"])
bw <- unname(fit_weibull$estimate["beta"])
lw <- unname(fit_weibull$estimate["lambda"])

# One to five years after surgery.
days <- c(365, 730, 1095, 1460, 1825)
na <- nelson_aalen(x)
na_at <- function(t) {
  idx <- findInterval(t, na$time)
  ifelse(idx == 0, 0, na$cumulative_hazard[pmax(idx, 1)])
}

survival_hazard <- data.frame(
  day = days,
  survival_empirical = vapply(days, function(d) mean(x > d), numeric(1)),
  survival_nme = snmeweibull(days, a, b, l),
  survival_weibull = sweibull_rate(days, bw, lw),
  cumhaz_nelson_aalen = na_at(days),
  cumhaz_nme = -snmeweibull(days, a, b, l, log = TRUE),
  cumhaz_weibull = lw * days^bw,
  hazard_nme = hnmeweibull(days, a, b, l),
  hazard_weibull = hweibull_rate(days, bw, lw)
)
save_table(survival_hazard, file.path(result_dir, "survival_hazard_at_days.csv"))

# Empirical hazard per six-month interval against both fitted models.
breaks <- seq(0, 2557.5, by = 182.5)
interval <- interval_hazard(x, breaks)
interval$hazard_nme <- hnmeweibull(interval$midpoint, a, b, l)
interval$hazard_weibull <- hweibull_rate(interval$midpoint, bw, lw)
save_table(interval, file.path(result_dir, "hazard_by_interval.csv"))

# Shape of the fitted NME-Weibull hazard over the observed range: location
# of its local maximum and minimum, if any.
grid_x <- seq(1, max(x), length.out = 5000)
h_grid <- hnmeweibull(grid_x, a, b, l)
turning <- which(diff(sign(diff(h_grid))) != 0) + 1
hazard_shape <- data.frame(
  day = grid_x[turning],
  hazard = h_grid[turning],
  type = ifelse(diff(sign(diff(h_grid)))[turning - 1] < 0, "maximum",
                "minimum")
)
if (nrow(hazard_shape) == 0) {
  hazard_shape <- data.frame(day = NA_real_, hazard = NA_real_,
                             type = "monotone on the observed range")
}
save_table(hazard_shape, file.path(result_dir, "hazard_turning_points.csv"))

cat("\nResults written to", result_dir, "and", fig_dir, "\n")
