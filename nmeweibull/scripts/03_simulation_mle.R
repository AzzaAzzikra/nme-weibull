# ============================================================
# 03. Monte Carlo simulation of the MLE in fitnmeweibull()
#     (Subbab 3.6.3, Tabel 4; ADEMP framework of Morris et al., 2019)
#
# Data-generating mechanisms:
#   S1: alpha = 0.5, beta = 1, lambda = 0.8  (hazard decreases to lambda)
#   S2: alpha = 1,   beta = 2, lambda = 0.5  (moderately increasing)
#   S3: alpha = 1.5, beta = 3, lambda = 0.3  (sharply increasing)
# Sample sizes n = 100, 500, 1000; M = 1000 replications each.
#
# Performance measures per parameter: bias and MSE (Equation (43)) with
# their Monte Carlo standard errors (Equation (44)), median estimate and
# median bias, convergence rate (code 0), and the proportion of boundary
# estimates (alpha_hat >= 0.99 * 1e4).
#
# The run takes roughly 40 minutes on 4 cores. Results of each
# scenario/sample-size cell are cached in results/03_simulation/raw/, so an
# interrupted run continues where it stopped. Set the environment variable
# NME_SIM_M to a small number (e.g. 20) for a quick trial run.
#
# Output: results/03_simulation/, figures/03_simulation/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "03_simulation")
raw_dir <- output_dir(result_dir, "raw")
fig_dir <- output_dir("figures", "03_simulation")

M <- as.integer(Sys.getenv("NME_SIM_M", "1000"))
n_cores <- max(1L, parallel::detectCores() - 1L)
alpha_bound <- 1e4

scenarios <- data.frame(
  scenario = c("S1", "S2", "S3"),
  alpha = c(0.5, 1.0, 1.5),
  beta = c(1.0, 2.0, 3.0),
  lambda = c(0.8, 0.5, 0.3)
)
sample_sizes <- c(100L, 500L, 1000L)

# One replication: the seed depends only on (scenario, n, replication), so
# results do not depend on the number of cores.
run_replication <- function(m, scenario_index, n, alpha, beta, lambda) {
  set.seed(20260000L + 100000L * scenario_index + n + m)
  x <- rnmeweibull(n, alpha, beta, lambda)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(fitnmeweibull(x), error = function(e) NULL)
  seconds <- proc.time()[["elapsed"]] - started

  if (is.null(fit)) {
    return(data.frame(replication = m, alpha_hat = NA_real_,
                      beta_hat = NA_real_, lambda_hat = NA_real_,
                      logLik = NA_real_, convergence = NA_integer_,
                      seconds = seconds))
  }

  data.frame(replication = m,
             alpha_hat = unname(fit$estimate["alpha"]),
             beta_hat = unname(fit$estimate["beta"]),
             lambda_hat = unname(fit$estimate["lambda"]),
             logLik = fit$logLik,
             convergence = fit$convergence,
             seconds = seconds)
}

cluster <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterEvalQ(cluster, {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(nmeweibull)
  }
  NULL
})

raw_all <- list()

for (i in seq_len(nrow(scenarios))) {
  for (n in sample_sizes) {
    sc <- scenarios[i, ]
    cache <- file.path(raw_dir, sprintf("%s_n%d_M%d.csv", sc$scenario, n, M))

    if (file.exists(cache)) {
      cell <- utils::read.csv(cache)
      cat("Loaded", cache, "\n")
    } else {
      cat(format(Sys.time(), "%H:%M:%S"), "Running", sc$scenario, "n =", n,
          "M =", M, "on", n_cores, "cores\n")
      cell <- do.call(rbind, parallel::parLapply(
        cluster, seq_len(M), run_replication, scenario_index = i, n = n,
        alpha = sc$alpha, beta = sc$beta, lambda = sc$lambda
      ))
      cell <- cbind(scenario = sc$scenario, n = n, alpha_true = sc$alpha,
                    beta_true = sc$beta, lambda_true = sc$lambda, cell)
      utils::write.csv(cell, cache, row.names = FALSE)
    }

    raw_all[[length(raw_all) + 1L]] <- cell
  }
}

raw <- do.call(rbind, raw_all)
raw$boundary <- raw$alpha_hat >= 0.99 * alpha_bound
utils::write.csv(raw, file.path(result_dir, "simulation_raw.csv"),
                 row.names = FALSE)

# ------------------------------------------------------------
# Performance measures
# ------------------------------------------------------------
performance <- function(estimates, true_value) {
  estimates <- estimates[is.finite(estimates)]
  M_ok <- length(estimates)
  errors <- estimates - true_value
  data.frame(
    mean_estimate = mean(estimates),
    bias = mean(errors),
    mcse_bias = stats::sd(estimates) / sqrt(M_ok),
    mse = mean(errors^2),
    mcse_mse = stats::sd(errors^2) / sqrt(M_ok),
    median_estimate = stats::median(estimates),
    median_bias = stats::median(estimates) - true_value
  )
}

cells <- unique(raw[, c("scenario", "n")])
summary_rows <- list()

for (j in seq_len(nrow(cells))) {
  cell <- raw[raw$scenario == cells$scenario[j] & raw$n == cells$n[j], ]
  interior <- cell[!cell$boundary & !is.na(cell$boundary), ]

  for (parameter in c("alpha", "beta", "lambda")) {
    true_value <- cell[[paste0(parameter, "_true")]][1]
    all_reps <- performance(cell[[paste0(parameter, "_hat")]], true_value)
    interior_reps <- performance(interior[[paste0(parameter, "_hat")]],
                                 true_value)

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      scenario = cells$scenario[j],
      n = cells$n[j],
      parameter = parameter,
      true_value = true_value,
      all_reps,
      bias_interior = interior_reps$bias,
      mse_interior = interior_reps$mse,
      replications = nrow(cell),
      convergence_rate = mean(cell$convergence == 0, na.rm = TRUE),
      boundary_rate = mean(cell$boundary, na.rm = TRUE),
      failed_fits = sum(is.na(cell$alpha_hat)),
      mean_seconds = mean(cell$seconds)
    )
  }
}

simulation_summary <- do.call(rbind, summary_rows)
simulation_summary <- simulation_summary[
  order(simulation_summary$scenario, simulation_summary$parameter,
        simulation_summary$n), ]
save_table(simulation_summary,
           file.path(result_dir, "simulation_summary.csv"), digits = 4)

# Compact table of convergence and boundary rates per cell.
stability <- unique(simulation_summary[, c("scenario", "n", "replications",
                                           "convergence_rate",
                                           "boundary_rate", "failed_fits",
                                           "mean_seconds")])
save_table(stability, file.path(result_dir, "simulation_stability.csv"),
           digits = 4)

# Check of the success criteria in Subbab 3.6.3.
criteria <- do.call(rbind, lapply(split(simulation_summary,
                                        list(simulation_summary$scenario,
                                             simulation_summary$parameter)),
                                  function(d) {
  d <- d[order(d$n), ]
  data.frame(
    scenario = d$scenario[1],
    parameter = d$parameter[1],
    abs_bias_decreasing = all(diff(abs(d$bias)) < 0),
    mse_decreasing = all(diff(d$mse) < 0),
    abs_median_bias_decreasing = all(diff(abs(d$median_bias)) <= 0),
    min_convergence_rate = min(d$convergence_rate),
    boundary_rate_non_increasing = all(diff(d$boundary_rate) <= 0)
  )
}))
save_table(criteria, file.path(result_dir, "simulation_criteria.csv"))

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------
cols <- c(S1 = "#1b6ca8", S2 = "#d1495b", S3 = "#2e8b57")
parameter_labels <- c(alpha = "alpha", beta = "beta", lambda = "lambda")

plot_measure <- function(measure, ylab, file, log_y = FALSE,
                         relative = FALSE) {
  open_png(file.path(fig_dir, file), width = 2400, height = 900)
  graphics::par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
  for (parameter in names(parameter_labels)) {
    d <- simulation_summary[simulation_summary$parameter == parameter, ]
    values <- d[[measure]]
    if (relative) values <- values / d$true_value
    if (log_y) values <- abs(values)
    graphics::plot(range(sample_sizes), range(values, finite = TRUE),
                   type = "n", log = if (log_y) "xy" else "x",
                   xlab = "Ukuran sampel (n)", ylab = ylab,
                   main = parameter_labels[parameter], xaxt = "n")
    graphics::axis(1, at = sample_sizes)
    if (!log_y) graphics::abline(h = 0, col = "grey70", lty = 3)
    for (sc in names(cols)) {
      ds <- d[d$scenario == sc, ]
      v <- if (relative) ds[[measure]] / ds$true_value else ds[[measure]]
      if (log_y) v <- abs(v)
      graphics::lines(ds$n, v, type = "b", pch = 19, col = cols[sc], lwd = 2)
    }
    graphics::legend("topright", names(cols), col = cols, lwd = 2, pch = 19,
                     bty = "n")
  }
  grDevices::dev.off()
}

plot_measure("mse", "MSE (skala log)", "mse_by_n.png", log_y = TRUE)
plot_measure("bias", "|Bias| (skala log)", "abs_bias_by_n.png", log_y = TRUE)
plot_measure("median_bias", "Bias median relatif (bias median / nilai aktual)",
             "relative_median_bias_by_n.png", relative = TRUE)

open_png(file.path(fig_dir, "boundary_rate_by_n.png"))
graphics::plot(range(sample_sizes), c(0, max(stability$boundary_rate, 0.01)),
               type = "n", log = "x", xaxt = "n", xlab = "Ukuran sampel (n)",
               ylab = "Proporsi estimasi batas alpha")
graphics::axis(1, at = sample_sizes)
for (sc in names(cols)) {
  ds <- stability[stability$scenario == sc, ]
  graphics::lines(ds$n, ds$boundary_rate, type = "b", pch = 19,
                  col = cols[sc], lwd = 2)
}
graphics::legend("topright", names(cols), col = cols, lwd = 2, pch = 19,
                 bty = "n")
grDevices::dev.off()

# Boxplots of the estimates of beta and lambda, and of log10(alpha_hat).
open_png(file.path(fig_dir, "estimate_boxplots.png"), width = 2400,
         height = 900)
graphics::par(mfrow = c(1, 3), mar = c(6, 4.5, 2.5, 1))
for (parameter in names(parameter_labels)) {
  values <- raw[[paste0(parameter, "_hat")]]
  ylab <- parameter_labels[parameter]
  if (parameter == "alpha") {
    values <- log10(values)
    ylab <- "log10(alpha)"
  }
  groups <- factor(paste(raw$scenario, raw$n),
                   levels = paste(rep(scenarios$scenario, each = 3),
                                  sample_sizes))
  graphics::boxplot(values ~ groups, las = 2, xlab = "", ylab = ylab,
                    col = rep(cols, each = 3), outline = FALSE,
                    main = parameter_labels[parameter])
  true_values <- scenarios[[parameter]]
  if (parameter == "alpha") true_values <- log10(true_values)
  graphics::segments(seq(0.6, by = 1, length.out = 9),
                     rep(true_values, each = 3),
                     seq(1.4, by = 1, length.out = 9),
                     rep(true_values, each = 3), lwd = 2, lty = 2)
}
grDevices::dev.off()

cat("\nResults written to", result_dir, "and", fig_dir, "\n")
