# ============================================================
# Shared helpers for the analysis scripts in scripts/
#
# Every script starts with
#   source("scripts/helpers.R")
# and must be run from the package root (the folder that contains
# DESCRIPTION, R/, data/, and scripts/), for example by opening
# nmeweibull.Rproj in RStudio.
# ============================================================

if (!file.exists("DESCRIPTION") || !dir.exists("scripts")) {
  stop("Run the scripts from the package root: open nmeweibull.Rproj first.",
       call. = FALSE)
}

# Load the development version of the package when pkgload/devtools is
# available; otherwise use the installed package.
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(nmeweibull)
}

# Create (if needed) and return an output folder.
output_dir <- function(...) {
  path <- file.path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Write a data frame to CSV without row names and print it.
save_table <- function(x, path, digits = 6) {
  utils::write.csv(x, path, row.names = FALSE)
  print(format(x, digits = digits), row.names = FALSE)
  invisible(x)
}

# Open a PNG device with the settings used for all thesis figures.
open_png <- function(path, width = 1600, height = 1000) {
  grDevices::png(path, width = width, height = height, res = 200)
}

# Moment-based sample skewness g1 = m3 / m2^(3/2).
sample_skewness <- function(x) {
  centred <- x - mean(x)
  mean(centred^3) / mean(centred^2)^(3 / 2)
}

# Descriptive statistics used in the real-data application (Subbab 3.7).
describe_data <- function(x) {
  data.frame(
    n = length(x),
    minimum = min(x),
    q1 = unname(stats::quantile(x, 0.25)),
    median = stats::median(x),
    mean = mean(x),
    q3 = unname(stats::quantile(x, 0.75)),
    maximum = max(x),
    sd = stats::sd(x),
    skewness = sample_skewness(x)
  )
}

# ------------------------------------------------------------
# Weibull maximum likelihood with the same method as fitnmeweibull():
# parameters on the log scale, data divided by their geometric mean,
# multi-start Nelder-Mead followed by BFGS. Parameterisation:
#   G(x) = 1 - exp(-lambda * x^beta).
# ------------------------------------------------------------
fit_weibull_mle <- function(x) {
  n <- length(x)
  log_m <- mean(log(x))
  y <- x / exp(log_m)
  log_y <- log(y)

  neg_loglik <- function(par) {
    beta <- exp(par[1])
    lambda_star <- exp(par[2])
    value <- -sum(log(beta) + log(lambda_star) + (beta - 1) * log_y -
                    lambda_star * exp(beta * log_y))
    if (!is.finite(value)) 1e100 else value
  }

  starts <- expand.grid(beta = c(0.5, 1, 1.5, 2, 3),
                        multiplier = c(0.25, 0.5, 1, 2, 4))
  best <- NULL

  for (i in seq_len(nrow(starts))) {
    beta0 <- starts$beta[i]
    # lambda0 = c / mean(x)^beta0, written on the internal scale.
    par0 <- c(log(beta0),
              log(starts$multiplier[i]) - beta0 * (log(mean(x)) - log_m))
    fit_nm <- stats::optim(par0, neg_loglik, method = "Nelder-Mead",
                           control = list(maxit = 5000, reltol = 1e-10))
    fit_bfgs <- stats::optim(fit_nm$par, neg_loglik, method = "BFGS",
                             control = list(maxit = 5000, reltol = 1e-10))
    candidate <- if (fit_bfgs$value <= fit_nm$value) fit_bfgs else fit_nm
    if (is.null(best) || candidate$value < best$value) best <- candidate
  }

  beta <- exp(best$par[1])
  lambda <- exp(best$par[2] - beta * log_m)
  loglik <- -best$value - n * log_m

  list(
    estimate = c(beta = beta, lambda = lambda),
    logLik = loglik,
    AIC = 2 * 2 - 2 * loglik,
    BIC = log(n) * 2 - 2 * loglik,
    convergence = best$convergence,
    n = n
  )
}

# Weibull functions in the (beta, lambda) parameterisation, via stats.
weibull_scale <- function(beta, lambda) lambda^(-1 / beta)

sweibull_rate <- function(x, beta, lambda) {
  stats::pweibull(x, beta, weibull_scale(beta, lambda), lower.tail = FALSE)
}

dweibull_rate <- function(x, beta, lambda) {
  stats::dweibull(x, beta, weibull_scale(beta, lambda))
}

hweibull_rate <- function(x, beta, lambda) beta * lambda * x^(beta - 1)

qweibull_rate <- function(p, beta, lambda) {
  stats::qweibull(p, beta, weibull_scale(beta, lambda))
}

# Weibull TVaR and RVaR by integrating the quantile function, with the
# same substitution as tvarnmeweibull() for the tail.
tvar_weibull_rate <- function(p, beta, lambda) {
  vapply(p, function(level) {
    integrand <- function(t) {
      s <- (1 - level) * exp(-t)
      out <- numeric(length(t))
      ok <- s > 0
      out[ok] <- stats::qweibull(s[ok], beta, weibull_scale(beta, lambda),
                                 lower.tail = FALSE) * exp(-t[ok])
      out
    }
    stats::integrate(integrand, 0, Inf)$value
  }, numeric(1))
}

rvar_weibull_rate <- function(p_lower, p_upper, beta, lambda) {
  vapply(seq_along(p_lower), function(i) {
    stats::integrate(qweibull_rate, p_lower[i], p_upper[i], beta = beta,
                     lambda = lambda)$value / (p_upper[i] - p_lower[i])
  }, numeric(1))
}

# Model comparison table with differences to the best model.
compare_models <- function(fit_nme, fit_weibull) {
  out <- data.frame(
    model = c("NME-Weibull", "Weibull"),
    k = c(3L, 2L),
    alpha = c(unname(fit_nme$estimate["alpha"]), NA_real_),
    beta = c(unname(fit_nme$estimate["beta"]),
             unname(fit_weibull$estimate["beta"])),
    lambda = c(unname(fit_nme$estimate["lambda"]),
               unname(fit_weibull$estimate["lambda"])),
    logLik = c(fit_nme$logLik, fit_weibull$logLik),
    AIC = c(fit_nme$AIC, fit_weibull$AIC),
    BIC = c(fit_nme$BIC, fit_weibull$BIC),
    convergence = c(fit_nme$convergence, fit_weibull$convergence),
    alpha_at_bound = c(isTRUE(fit_nme$boundary), NA)
  )
  out$delta_AIC <- out$AIC - min(out$AIC)
  out$delta_BIC <- out$BIC - min(out$BIC)
  out
}

# Empirical TVaR and RVaR from the sorted sample (for comparison only).
empirical_tvar <- function(x, p) {
  vapply(p, function(level) mean(x[x > stats::quantile(x, level)]),
         numeric(1))
}

empirical_rvar <- function(x, p_lower, p_upper) {
  vapply(seq_along(p_lower), function(i) {
    lo <- stats::quantile(x, p_lower[i])
    hi <- stats::quantile(x, p_upper[i])
    mean(x[x > lo & x <= hi])
  }, numeric(1))
}

# Risk measures of both fitted models at the levels of Subbab 3.7.
risk_measures <- function(x, fit_nme, fit_weibull,
                          p = c(0.90, 0.95, 0.99),
                          p_lower = c(0.90, 0.95),
                          p_upper = c(0.95, 0.99)) {
  a <- unname(fit_nme$estimate["alpha"])
  b <- unname(fit_nme$estimate["beta"])
  l <- unname(fit_nme$estimate["lambda"])
  bw <- unname(fit_weibull$estimate["beta"])
  lw <- unname(fit_weibull$estimate["lambda"])

  var_tvar <- data.frame(
    p = p,
    VaR_nme = varnmeweibull(p, a, b, l),
    TVaR_nme = tvarnmeweibull(p, a, b, l),
    VaR_weibull = qweibull_rate(p, bw, lw),
    TVaR_weibull = tvar_weibull_rate(p, bw, lw),
    VaR_empirical = unname(stats::quantile(x, p)),
    TVaR_empirical = empirical_tvar(x, p)
  )

  rvar <- data.frame(
    p_lower = p_lower,
    p_upper = p_upper,
    RVaR_nme = rvarnmeweibull(p_lower, p_upper, a, b, l),
    RVaR_weibull = rvar_weibull_rate(p_lower, p_upper, bw, lw),
    RVaR_empirical = empirical_rvar(x, p_lower, p_upper)
  )

  list(var_tvar = var_tvar, rvar = rvar)
}

# Figures comparing the two fitted models with the data.
plot_fitted_models <- function(x, fit_nme, fit_weibull, fig_dir, unit_label,
                               prefix = "") {
  a <- unname(fit_nme$estimate["alpha"])
  b <- unname(fit_nme$estimate["beta"])
  l <- unname(fit_nme$estimate["lambda"])
  bw <- unname(fit_weibull$estimate["beta"])
  lw <- unname(fit_weibull$estimate["lambda"])
  grid_x <- seq(min(x) / 2, max(x), length.out = 400)
  cols <- c(nme = "#1b6ca8", weibull = "#d1495b")

  open_png(file.path(fig_dir, paste0(prefix, "fitted_pdf.png")))
  graphics::hist(x, breaks = "FD", freq = FALSE, col = "grey90",
                 border = "grey60", main = "", xlab = unit_label,
                 ylab = "Kepadatan",
                 ylim = c(0, 1.1 * max(dnmeweibull(grid_x, a, b, l),
                                       dweibull_rate(grid_x, bw, lw),
                                       graphics::hist(x, breaks = "FD",
                                                      plot = FALSE)$density)))
  graphics::lines(grid_x, dnmeweibull(grid_x, a, b, l), col = cols["nme"],
                  lwd = 2)
  graphics::lines(grid_x, dweibull_rate(grid_x, bw, lw),
                  col = cols["weibull"], lwd = 2, lty = 2)
  graphics::legend("topright", c("NME-Weibull", "Weibull"),
                   col = cols, lwd = 2, lty = c(1, 2), bty = "n")
  grDevices::dev.off()

  sorted <- sort(x)
  empirical_survival <- 1 - seq_along(sorted) / length(sorted)
  open_png(file.path(fig_dir, paste0(prefix, "fitted_survival.png")))
  graphics::plot(c(0, sorted), c(1, empirical_survival), type = "s",
                 xlab = unit_label, ylab = "S(x)", main = "",
                 col = "grey30")
  graphics::lines(grid_x, snmeweibull(grid_x, a, b, l), col = cols["nme"],
                  lwd = 2)
  graphics::lines(grid_x, sweibull_rate(grid_x, bw, lw),
                  col = cols["weibull"], lwd = 2, lty = 2)
  graphics::legend("topright", c("Empiris", "NME-Weibull", "Weibull"),
                   col = c("grey30", cols), lwd = c(1, 2, 2),
                   lty = c(1, 1, 2), bty = "n")
  grDevices::dev.off()

  open_png(file.path(fig_dir, paste0(prefix, "fitted_hazard.png")))
  h_nme <- hnmeweibull(grid_x, a, b, l)
  h_wei <- hweibull_rate(grid_x, bw, lw)
  graphics::plot(grid_x, h_nme, type = "l", col = cols["nme"], lwd = 2,
                 ylim = range(c(h_nme, h_wei), finite = TRUE),
                 xlab = unit_label, ylab = "h(x)", main = "")
  graphics::lines(grid_x, h_wei, col = cols["weibull"], lwd = 2, lty = 2)
  graphics::legend("topright", c("NME-Weibull", "Weibull"), col = cols,
                   lwd = 2, lty = c(1, 2), bty = "n")
  grDevices::dev.off()

  invisible(NULL)
}
