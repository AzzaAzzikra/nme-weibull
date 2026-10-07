# ============================================================
# RVaR Optimization on Real Data
# ============================================================

devtools::load_all()

dir.create("results", showWarnings = FALSE)

model_comparison <- read.csv("results/real_data_model_comparison.csv")

nme_row <- model_comparison[model_comparison$model == "NME-Weibull", ]

alpha_hat <- nme_row$alpha
beta_hat <- nme_row$beta
lambda_hat <- nme_row$lambda

# ------------------------------------------------------------
# Optimize upper bound of RVaR interval
# ------------------------------------------------------------

optimize_rvarnmeweibull <- function(
    p_lower,
    p_upper_range,
    alpha,
    beta,
    lambda,
    rel.tol = .Machine$double.eps^0.25,
    subdivisions = 100L
) {
  if (length(p_upper_range) != 2) {
    stop("p_upper_range must contain two values: lower and upper bounds.")
  }

  p_upper_min <- p_upper_range[1]
  p_upper_max <- p_upper_range[2]

  if (p_lower <= 0 || p_lower >= 1) {
    stop("p_lower must be between 0 and 1.")
  }

  if (p_upper_min <= p_lower || p_upper_max <= p_lower) {
    stop("p_upper values must be greater than p_lower.")
  }

  if (p_upper_min >= p_upper_max) {
    stop("p_upper_range must be increasing.")
  }

  objective <- function(p_upper) {
    rvarnmeweibull(
      p_lower = p_lower,
      p_upper = p_upper,
      alpha = alpha,
      beta = beta,
      lambda = lambda,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    )
  }

  opt <- stats::optimize(
    f = objective,
    interval = c(p_upper_min, p_upper_max),
    maximum = FALSE
  )

  # Because RVaR generally increases with p_upper,
  # check the boundary values explicitly.
  candidates <- c(p_upper_min, opt$minimum, p_upper_max)

  values <- sapply(candidates, objective)

  best_index <- which.min(values)

  p_upper_opt <- candidates[best_index]
  rvar_opt <- values[best_index]

  var_lower <- varnmeweibull(
    p = p_lower,
    alpha = alpha,
    beta = beta,
    lambda = lambda
  )

  var_upper <- varnmeweibull(
    p = p_upper_opt,
    alpha = alpha,
    beta = beta,
    lambda = lambda
  )

  tvar_lower <- tvarnmeweibull(
    p = p_lower,
    alpha = alpha,
    beta = beta,
    lambda = lambda,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )

  data.frame(
    p_lower = p_lower,
    p_upper_min = p_upper_min,
    p_upper_max = p_upper_max,
    p_upper_opt = p_upper_opt,
    VaR_lower = var_lower,
    RVaR_opt = rvar_opt,
    VaR_upper_opt = var_upper,
    TVaR_lower = tvar_lower,
    saving_vs_TVaR = tvar_lower - rvar_opt,
    saving_pct_vs_TVaR = 100 * (tvar_lower - rvar_opt) / tvar_lower
  )
}

# ------------------------------------------------------------
# Optimization scenarios
# ------------------------------------------------------------

rvar_optimization <- do.call(
  rbind,
  list(
    optimize_rvarnmeweibull(
      p_lower = 0.90,
      p_upper_range = c(0.95, 0.99),
      alpha = alpha_hat,
      beta = beta_hat,
      lambda = lambda_hat
    ),
    optimize_rvarnmeweibull(
      p_lower = 0.95,
      p_upper_range = c(0.96, 0.99),
      alpha = alpha_hat,
      beta = beta_hat,
      lambda = lambda_hat
    )
  )
)

write.csv(
  rvar_optimization,
  "results/real_data_rvar_optimization.csv",
  row.names = FALSE
)

print(rvar_optimization)