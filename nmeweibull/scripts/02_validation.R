# ============================================================
# 02. Statistical validation rules V1-V9 (Subbab 3.6, Tabel 3)
#
# Each rule is applied to K = 48 parameter combinations:
#   alpha in {0.1, 0.5, 1, 5}, beta in {0.5, 1, 2, 3},
#   lambda in {0.5, 1, 2}.
# For each combination the evaluation points are 100 equally spaced
# values from Q(0.001) to Q(0.999). A rule passes when the maximum error
# e_max does not exceed its tolerance:
#   tau1 = sqrt(eps)  ~ 1.49e-8  (closed-form formulas)
#   tau2 = eps^(1/4)  ~ 1.22e-4  (numerical integration)
#   tau3 = 1e-6                  (numerical derivative and limit case)
#
# Output: results/02_validation/
# ============================================================

source("scripts/helpers.R")

result_dir <- output_dir("results", "02_validation")

eps <- .Machine$double.eps
tau1 <- sqrt(eps)
tau2 <- eps^(1 / 4)
tau3 <- 1e-6

grid <- expand.grid(lambda = c(0.5, 1, 2), beta = c(0.5, 1, 2, 3),
                    alpha = c(0.1, 0.5, 1, 5))
grid <- grid[, c("alpha", "beta", "lambda")]
K <- nrow(grid)

p_grid <- seq(0.01, 0.99, by = 0.01)
risk_p <- c(0.90, 0.95, 0.99)
risk_lower <- c(0.90, 0.95)
risk_upper <- c(0.95, 0.99)
ks_level <- 0.05 / K
ks_n <- 10000

relative_error <- function(a, b) max(abs(a - b) / abs(b))

validate_one <- function(alpha, beta, lambda, index) {
  x <- seq(qnmeweibull(0.001, alpha, beta, lambda),
           qnmeweibull(0.999, alpha, beta, lambda), length.out = 100)
  pdf <- dnmeweibull(x, alpha, beta, lambda)
  cdf <- pnmeweibull(x, alpha, beta, lambda)
  surv <- snmeweibull(x, alpha, beta, lambda)
  haz <- hnmeweibull(x, alpha, beta, lambda)

  # V1: f >= 0 and the density integrates to one. The integral is split at
  # the median because the density is unbounded at zero when beta < 1.
  median_x <- qnmeweibull(0.5, alpha, beta, lambda)
  total <- stats::integrate(dnmeweibull, 0, median_x, alpha = alpha,
                            beta = beta, lambda = lambda)$value +
    stats::integrate(dnmeweibull, median_x, Inf, alpha = alpha, beta = beta,
                     lambda = lambda)$value
  v1_error <- abs(total - 1)
  v1_pass <- all(pdf >= 0) && v1_error <= tau2

  # V2: F(0) = 0, F non-decreasing, and F(Q(0.999)) = 0.999.
  v2_error <- max(abs(pnmeweibull(0, alpha, beta, lambda)),
                  abs(pnmeweibull(qnmeweibull(0.999, alpha, beta, lambda),
                                  alpha, beta, lambda) - 0.999))
  v2_pass <- all(diff(cdf) >= 0) && v2_error <= tau1

  # V3: S(x) = 1 - F(x).
  v3_error <- max(abs(surv - (1 - cdf)))
  v3_pass <- v3_error <= tau1

  # V4: h = f / S (relative, tau1) and f equals the central difference of
  # F with step delta = 1e-5 * x (relative, tau3).
  v4a_error <- relative_error(haz, pdf / surv)
  delta <- 1e-5 * x
  numerical_pdf <- (pnmeweibull(x + delta, alpha, beta, lambda) -
                      pnmeweibull(x - delta, alpha, beta, lambda)) / (2 * delta)
  v4b_error <- relative_error(numerical_pdf, pdf)
  v4_pass <- v4a_error <= tau1 && v4b_error <= tau3

  # V5: F(Q(p)) = p and Q(F(x)) = x (relative, tau1).
  v5a_error <- relative_error(
    pnmeweibull(qnmeweibull(p_grid, alpha, beta, lambda), alpha, beta,
                lambda), p_grid)
  v5b_error <- relative_error(qnmeweibull(cdf, alpha, beta, lambda), x)
  v5_pass <- max(v5a_error, v5b_error) <= tau1

  # V7: Kolmogorov-Smirnov test on n = 10,000 random values, Bonferroni
  # corrected level 0.05 / 48. The seed is fixed per combination.
  set.seed(20260000 + index)
  ks <- stats::ks.test(rnmeweibull(ks_n, alpha, beta, lambda), pnmeweibull,
                       alpha = alpha, beta = beta, lambda = lambda)
  v7_pass <- ks$p.value > ks_level

  # V8: VaR equals Q(p), increases with p, and
  # VaR(p1) <= RVaR(p1, p2) <= VaR(p2).
  var_p <- varnmeweibull(risk_p, alpha, beta, lambda)
  v8a_error <- max(abs(var_p - qnmeweibull(risk_p, alpha, beta, lambda)))
  rvar <- rvarnmeweibull(risk_lower, risk_upper, alpha, beta, lambda)
  var_lower <- varnmeweibull(risk_lower, alpha, beta, lambda)
  var_upper <- varnmeweibull(risk_upper, alpha, beta, lambda)
  v8_pass <- v8a_error <= tau1 && all(diff(var_p) > 0) &&
    all(var_lower <= rvar) && all(rvar <= var_upper)

  # V9: TVaR from Equation (31) equals TVaR from Equation (32)
  # (relative, tau2), and TVaR >= VaR.
  tvar_quantile <- tvarnmeweibull(risk_p, alpha, beta, lambda)
  tvar_survival <- var_p + vapply(seq_along(risk_p), function(i) {
    stats::integrate(snmeweibull, var_p[i], Inf, alpha = alpha, beta = beta,
                     lambda = lambda)$value / (1 - risk_p[i])
  }, numeric(1))
  v9_error <- relative_error(tvar_quantile, tvar_survival)
  v9_pass <- v9_error <= tau2 && all(tvar_quantile >= var_p)

  data.frame(
    alpha = alpha, beta = beta, lambda = lambda,
    V1_error = v1_error, V1_pass = v1_pass,
    V2_error = v2_error, V2_pass = v2_pass,
    V3_error = v3_error, V3_pass = v3_pass,
    V4_ratio_error = v4a_error, V4_derivative_error = v4b_error,
    V4_pass = v4_pass,
    V5_FQ_error = v5a_error, V5_QF_error = v5b_error, V5_pass = v5_pass,
    V7_ks_statistic = unname(ks$statistic), V7_p_value = ks$p.value,
    V7_pass = v7_pass,
    V8_error = v8a_error, V8_pass = v8_pass,
    V9_error = v9_error, V9_pass = v9_pass
  )
}

details <- do.call(rbind, lapply(seq_len(K), function(i) {
  validate_one(grid$alpha[i], grid$beta[i], grid$lambda[i], i)
}))

# V6: limit case. With alpha = 1e8 the NME-Weibull CDF must match pweibull()
# with scale = lambda^(-1/beta). Alpha does not enter, so the 12 (beta,
# lambda) pairs are checked.
limit_grid <- unique(grid[, c("beta", "lambda")])
limit_alpha <- 1e8
v6 <- do.call(rbind, lapply(seq_len(nrow(limit_grid)), function(i) {
  b <- limit_grid$beta[i]
  l <- limit_grid$lambda[i]
  x <- seq(qnmeweibull(0.001, limit_alpha, b, l),
           qnmeweibull(0.999, limit_alpha, b, l), length.out = 100)
  error <- max(abs(pnmeweibull(x, limit_alpha, b, l) -
                     stats::pweibull(x, shape = b, scale = l^(-1 / b))))
  data.frame(beta = b, lambda = l, V6_error = error, V6_pass = error <= tau3)
}))

utils::write.csv(details, file.path(result_dir, "validation_details.csv"),
                 row.names = FALSE)
utils::write.csv(v6, file.path(result_dir, "validation_limit_case.csv"),
                 row.names = FALSE)

summary_row <- function(code, aspect, error, pass, tolerance, note = "") {
  data.frame(
    rule = code,
    aspect = aspect,
    combinations = length(pass),
    passed = sum(pass),
    max_error = if (is.null(error)) NA_real_ else max(error),
    tolerance = tolerance,
    note = note,
    status = if (all(pass)) "Terpenuhi" else "Tidak terpenuhi"
  )
}

validation_summary <- rbind(
  summary_row("V1", "Sifat PDF", details$V1_error, details$V1_pass, tau2,
              "|integral f - 1|"),
  summary_row("V2", "Sifat CDF", details$V2_error, details$V2_pass, tau1),
  summary_row("V3", "Fungsi survival", details$V3_error, details$V3_pass,
              tau1),
  summary_row("V4", "Hazard = f/S", details$V4_ratio_error, details$V4_pass,
              tau1, "galat relatif"),
  summary_row("V4", "PDF = turunan numerik F", details$V4_derivative_error,
              details$V4_pass, tau3, "galat relatif"),
  summary_row("V5", "F(Q(p)) = p", details$V5_FQ_error, details$V5_pass,
              tau1, "galat relatif"),
  summary_row("V5", "Q(F(x)) = x", details$V5_QF_error, details$V5_pass,
              tau1, "galat relatif"),
  summary_row("V6", "Kasus limit Weibull", v6$V6_error, v6$V6_pass, tau3,
              "alpha = 1e8, 12 pasangan (beta, lambda)"),
  summary_row("V7", "Data acak (uji KS)", NULL, details$V7_pass, ks_level,
              sprintf("nilai-p minimum = %.4f", min(details$V7_p_value))),
  summary_row("V8", "VaR dan RVaR", details$V8_error, details$V8_pass, tau1),
  summary_row("V9", "TVaR persamaan (31) = (32)", details$V9_error,
              details$V9_pass, tau2, "galat relatif")
)

save_table(validation_summary,
           file.path(result_dir, "validation_summary.csv"), digits = 3)

cat("\nAll rules satisfied:", all(validation_summary$status == "Terpenuhi"),
    "\n")
