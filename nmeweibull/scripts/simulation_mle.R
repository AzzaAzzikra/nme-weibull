# ============================================================
# Simulation Study for NME-Weibull MLE
# ============================================================

library(nmeweibull)

set.seed(123)

# Parameter scenarios
scenarios <- data.frame(
  scenario = c("S1", "S2", "S3"),
  alpha = c(0.5, 1.0, 1.5),
  beta = c(1.0, 2.0, 3.0),
  lambda = c(0.8, 0.5, 0.3)
)

# Sample sizes and number of replications
sample_sizes <- c(100, 500, 1000)
B <- 30

simulation_results <- list()
index <- 1

for (s in seq_len(nrow(scenarios))) {
  for (n in sample_sizes) {
    for (b in seq_len(B)) {

      alpha_true <- scenarios$alpha[s]
      beta_true <- scenarios$beta[s]
      lambda_true <- scenarios$lambda[s]

      x <- rnmeweibull(
        n = n,
        alpha = alpha_true,
        beta = beta_true,
        lambda = lambda_true
      )

      fit <- try(
        fitnmeweibull(x),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        simulation_results[[index]] <- data.frame(
          scenario = scenarios$scenario[s],
          n = n,
          replication = b,
          alpha_true = alpha_true,
          beta_true = beta_true,
          lambda_true = lambda_true,
          alpha_hat = NA_real_,
          beta_hat = NA_real_,
          lambda_hat = NA_real_,
          logLik = NA_real_,
          AIC = NA_real_,
          BIC = NA_real_,
          convergence = NA_integer_
        )
      } else {
        simulation_results[[index]] <- data.frame(
          scenario = scenarios$scenario[s],
          n = n,
          replication = b,
          alpha_true = alpha_true,
          beta_true = beta_true,
          lambda_true = lambda_true,
          alpha_hat = unname(fit$estimate["alpha"]),
          beta_hat = unname(fit$estimate["beta"]),
          lambda_hat = unname(fit$estimate["lambda"]),
          logLik = fit$logLik,
          AIC = fit$AIC,
          BIC = fit$BIC,
          convergence = fit$convergence
        )
      }

      index <- index + 1
    }
  }
}

simulation_results <- do.call(rbind, simulation_results)

# Save raw simulation results
dir.create("results", showWarnings = FALSE)
write.csv(
  simulation_results,
  "results/simulation_mle_raw.csv",
  row.names = FALSE
)

# Function to calculate summary statistics
summarise_parameter <- function(estimate, true_value) {
  bias <- mean(estimate - true_value, na.rm = TRUE)
  abs_bias <- mean(abs(estimate - true_value), na.rm = TRUE)
  mse <- mean((estimate - true_value)^2, na.rm = TRUE)

  c(
    mean_estimate = mean(estimate, na.rm = TRUE),
    bias = bias,
    abs_bias = abs_bias,
    mse = mse
  )
}

# Build summary table
summary_results <- do.call(
  rbind,
  lapply(split(simulation_results, list(simulation_results$scenario, simulation_results$n)), function(df) {

    alpha_summary <- summarise_parameter(df$alpha_hat, df$alpha_true[1])
    beta_summary <- summarise_parameter(df$beta_hat, df$beta_true[1])
    lambda_summary <- summarise_parameter(df$lambda_hat, df$lambda_true[1])

    data.frame(
      scenario = df$scenario[1],
      n = df$n[1],
      parameter = c("alpha", "beta", "lambda"),
      true_value = c(df$alpha_true[1], df$beta_true[1], df$lambda_true[1]),
      mean_estimate = c(
        alpha_summary["mean_estimate"],
        beta_summary["mean_estimate"],
        lambda_summary["mean_estimate"]
      ),
      bias = c(
        alpha_summary["bias"],
        beta_summary["bias"],
        lambda_summary["bias"]
      ),
      abs_bias = c(
        alpha_summary["abs_bias"],
        beta_summary["abs_bias"],
        lambda_summary["abs_bias"]
      ),
      mse = c(
        alpha_summary["mse"],
        beta_summary["mse"],
        lambda_summary["mse"]
      ),
      convergence_rate = mean(df$convergence == 0, na.rm = TRUE),
      mean_logLik = mean(df$logLik, na.rm = TRUE),
      mean_AIC = mean(df$AIC, na.rm = TRUE),
      mean_BIC = mean(df$BIC, na.rm = TRUE)
    )
  })
)

rownames(summary_results) <- NULL

write.csv(
  summary_results,
  "results/simulation_mle_summary.csv",
  row.names = FALSE
)

print(summary_results)
