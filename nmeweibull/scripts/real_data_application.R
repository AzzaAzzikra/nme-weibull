# ============================================================
# Real Data Application: NME-Weibull vs Weibull
# ============================================================

devtools::load_all()

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read real data
# ------------------------------------------------------------

real_data <- read.csv("data/real_data.csv")

# Data column must be named "time"
x <- real_data$time
x <- x[is.finite(x) & x > 0]

if (length(x) == 0) {
  stop("Data must contain at least one positive finite value.")
}

n <- length(x)

# ------------------------------------------------------------
# 2. Descriptive statistics
# ------------------------------------------------------------

skewness_base <- function(z) {
  m <- mean(z)
  s <- stats::sd(z)
  mean(((z - m) / s)^3)
}

kurtosis_base <- function(z) {
  m <- mean(z)
  s <- stats::sd(z)
  mean(((z - m) / s)^4)
}

descriptive_stats <- data.frame(
  n = n,
  minimum = min(x),
  q1 = unname(quantile(x, 0.25)),
  median = median(x),
  mean = mean(x),
  q3 = unname(quantile(x, 0.75)),
  maximum = max(x),
  variance = stats::var(x),
  skewness = skewness_base(x),
  kurtosis = kurtosis_base(x)
)

write.csv(
  descriptive_stats,
  "results/real_data_descriptive.csv",
  row.names = FALSE
)

print(descriptive_stats)

# ------------------------------------------------------------
# 3. Fit Weibull baseline model
# Parameterization:
# F(x) = 1 - exp(-lambda * x^beta)
# f(x) = beta * lambda * x^(beta - 1) * exp(-lambda * x^beta)
# ------------------------------------------------------------

fitweibull_baseline <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }

  x <- x[is.finite(x)]

  if (any(x <= 0)) {
    stop("All observations must be positive.")
  }

  neg_loglik <- function(log_par) {
    par <- exp(log_par)

    beta <- par[1]
    lambda <- par[2]

    log_density <-
      log(beta) +
      log(lambda) +
      (beta - 1) * log(x) -
      lambda * x^beta

    if (any(!is.finite(log_density))) {
      return(1e100)
    }

    -sum(log_density)
  }

  beta_grid <- c(0.5, 1, 1.5, 2, 3, 5)
  lambda_grid <- c(
    0.25 / mean(x),
    0.5 / mean(x),
    1 / mean(x),
    2 / mean(x),
    4 / mean(x)
  )

  starts <- expand.grid(beta = beta_grid, lambda = lambda_grid)

  best_fit <- NULL

  for (i in seq_len(nrow(starts))) {
    fit <- try(
      stats::optim(
        par = log(as.numeric(starts[i, ])),
        fn = neg_loglik,
        method = "BFGS",
        control = list(maxit = 5000, reltol = 1e-10)
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || !is.finite(fit$value)) {
      next
    }

    if (is.null(best_fit) || fit$value < best_fit$value) {
      best_fit <- fit
    }
  }

  if (is.null(best_fit)) {
    stop("Weibull optimization failed.")
  }

  estimate <- exp(best_fit$par)
  names(estimate) <- c("beta", "lambda")

  loglik <- -best_fit$value
  k <- length(estimate)
  n <- length(x)

  list(
    estimate = estimate,
    logLik = loglik,
    AIC = 2 * k - 2 * loglik,
    BIC = log(n) * k - 2 * loglik,
    convergence = best_fit$convergence,
    optim = best_fit
  )
}

dweibull_baseline <- function(x, beta, lambda) {
  density <- rep(NA_real_, length(x))
  density[x < 0 & !is.na(x)] <- 0

  valid <- x >= 0 & is.finite(x) & !is.na(x)

  density[valid] <-
    beta * lambda * x[valid]^(beta - 1) *
    exp(-lambda * x[valid]^beta)

  density
}

pweibull_baseline <- function(q, beta, lambda) {
  cdf <- rep(NA_real_, length(q))
  cdf[q < 0 & !is.na(q)] <- 0
  cdf[is.infinite(q) & q > 0] <- 1

  valid <- q >= 0 & is.finite(q) & !is.na(q)

  cdf[valid] <- 1 - exp(-lambda * q[valid]^beta)

  cdf
}

sweibull_baseline <- function(x, beta, lambda) {
  1 - pweibull_baseline(x, beta, lambda)
}

# ------------------------------------------------------------
# 4. Fit NME-Weibull and Weibull
# ------------------------------------------------------------

fit_nme <- fitnmeweibull(x)
fit_weibull <- fitweibull_baseline(x)

model_comparison <- data.frame(
  model = c("NME-Weibull", "Weibull"),
  alpha = c(unname(fit_nme$estimate["alpha"]), NA),
  beta = c(unname(fit_nme$estimate["beta"]), unname(fit_weibull$estimate["beta"])),
  lambda = c(unname(fit_nme$estimate["lambda"]), unname(fit_weibull$estimate["lambda"])),
  logLik = c(fit_nme$logLik, fit_weibull$logLik),
  AIC = c(fit_nme$AIC, fit_weibull$AIC),
  BIC = c(fit_nme$BIC, fit_weibull$BIC),
  convergence = c(fit_nme$convergence, fit_weibull$convergence)
)

write.csv(
  model_comparison,
  "results/real_data_model_comparison.csv",
  row.names = FALSE
)

print(model_comparison)

# ------------------------------------------------------------
# 5. Risk measures using fitted NME-Weibull parameters
# ------------------------------------------------------------

alpha_hat <- unname(fit_nme$estimate["alpha"])
beta_hat <- unname(fit_nme$estimate["beta"])
lambda_hat <- unname(fit_nme$estimate["lambda"])

risk_levels <- c(0.90, 0.95, 0.99)

risk_table <- data.frame(
  p = risk_levels,
  VaR = varnmeweibull(
    p = risk_levels,
    alpha = alpha_hat,
    beta = beta_hat,
    lambda = lambda_hat
  ),
  TVaR = tvarnmeweibull(
    p = risk_levels,
    alpha = alpha_hat,
    beta = beta_hat,
    lambda = lambda_hat
  )
)

rvar_table <- data.frame(
  p_lower = c(0.90, 0.95),
  p_upper = c(0.95, 0.99)
)

rvar_table$RVaR <- rvarnmeweibull(
  p_lower = rvar_table$p_lower,
  p_upper = rvar_table$p_upper,
  alpha = alpha_hat,
  beta = beta_hat,
  lambda = lambda_hat
)

write.csv(
  risk_table,
  "results/real_data_var_tvar.csv",
  row.names = FALSE
)

write.csv(
  rvar_table,
  "results/real_data_rvar.csv",
  row.names = FALSE
)

print(risk_table)
print(rvar_table)

# ------------------------------------------------------------
# 6. Plots: fitted PDF, CDF, and survival function
# ------------------------------------------------------------

grid_x <- seq(min(x), max(x), length.out = 300)

nme_pdf <- dnmeweibull(
  x = grid_x,
  alpha = alpha_hat,
  beta = beta_hat,
  lambda = lambda_hat
)

weibull_pdf <- dweibull_baseline(
  x = grid_x,
  beta = unname(fit_weibull$estimate["beta"]),
  lambda = unname(fit_weibull$estimate["lambda"])
)

png(
  filename = "figures/real_data_fitted_pdf.png",
  width = 1000,
  height = 700,
  res = 150
)

hist(
  x,
  probability = TRUE,
  breaks = "FD",
  xlab = "Data",
  main = "Fitted PDF pada Data Riil"
)

lines(grid_x, nme_pdf, lwd = 2)
lines(grid_x, weibull_pdf, lwd = 2, lty = 2)

legend(
  "topright",
  legend = c("NME-Weibull", "Weibull"),
  lty = c(1, 2),
  lwd = 2
)

dev.off()

png(
  filename = "figures/real_data_fitted_cdf.png",
  width = 1000,
  height = 700,
  res = 150
)

plot(
  stats::ecdf(x),
  main = "Fitted CDF pada Data Riil",
  xlab = "Data",
  ylab = "F(x)"
)

lines(
  grid_x,
  pnmeweibull(grid_x, alpha = alpha_hat, beta = beta_hat, lambda = lambda_hat),
  lwd = 2
)

lines(
  grid_x,
  pweibull_baseline(
    grid_x,
    beta = unname(fit_weibull$estimate["beta"]),
    lambda = unname(fit_weibull$estimate["lambda"])
  ),
  lwd = 2,
  lty = 2
)

legend(
  "bottomright",
  legend = c("Empirical CDF", "NME-Weibull", "Weibull"),
  lty = c(1, 1, 2),
  lwd = c(1, 2, 2)
)

dev.off()

png(
  filename = "figures/real_data_fitted_survival.png",
  width = 1000,
  height = 700,
  res = 150
)

plot(
  sort(x),
  1 - seq_along(sort(x)) / length(x),
  type = "s",
  xlab = "Data",
  ylab = "S(x)",
  main = "Fitted Survival Function pada Data Riil"
)

lines(
  grid_x,
  snmeweibull(grid_x, alpha = alpha_hat, beta = beta_hat, lambda = lambda_hat),
  lwd = 2
)

lines(
  grid_x,
  sweibull_baseline(
    grid_x,
    beta = unname(fit_weibull$estimate["beta"]),
    lambda = unname(fit_weibull$estimate["lambda"])
  ),
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c("Empirical Survival", "NME-Weibull", "Weibull"),
  lty = c(1, 1, 2),
  lwd = c(1, 2, 2)
)

dev.off()
