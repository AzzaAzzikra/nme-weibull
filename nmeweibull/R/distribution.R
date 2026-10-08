#' Check NME-Weibull Parameters
#'
#' Internal helper that stops with an informative error unless `alpha`,
#' `beta`, and `lambda` are single, finite, positive numbers.
#'
#' @param alpha,beta,lambda Parameters of the NME-Weibull distribution.
#'
#' @return `NULL`, invisibly. Called for its side effect.
#' @noRd
.check_nmeweibull_params <- function(alpha, beta, lambda) {
  params <- list(alpha = alpha, beta = beta, lambda = lambda)

  for (name in names(params)) {
    value <- params[[name]]

    if (!is.numeric(value) || length(value) != 1L) {
      stop(name, " must be a single numeric value.", call. = FALSE)
    }

    if (is.na(value) || !is.finite(value)) {
      stop(name, " must be finite.", call. = FALSE)
    }

    if (value <= 0) {
      stop(name, " must be positive.", call. = FALSE)
    }
  }

  invisible(NULL)
}


#' Check a Numeric Argument
#'
#' @param x Object to check.
#' @param name Argument name used in the error message.
#'
#' @return `NULL`, invisibly.
#' @noRd
.check_numeric_arg <- function(x, name) {
  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }

  invisible(NULL)
}


#' Core Log-Scale Quantities of the NME-Weibull Distribution
#'
#' Internal helper that evaluates the log-survival, log-density, and
#' log-hazard of the NME-Weibull distribution at strictly positive, finite
#' points. All quantities are computed on the log scale for numerical
#' stability, using
#' \deqn{\eta = \lambda x^\beta, \quad \bar{G} = e^{-\eta}, \quad
#'   G = 1 - e^{-\eta},}
#' so that \eqn{\alpha + 2 - G = \alpha + 1 + \bar{G}}.
#'
#' @param x Numeric vector of strictly positive, finite values.
#' @param alpha,beta,lambda Parameters of the NME-Weibull distribution.
#'
#' @return A list with elements `log_surv`, `log_dens`, and `log_haz`.
#' @noRd
.nmeweibull_core <- function(x, alpha, beta, lambda) {
  log_x <- log(x)
  eta <- lambda * exp(beta * log_x)
  surv_base <- exp(-eta)
  cdf_base <- -expm1(-eta)

  # (beta - 1) * log(x) written so that beta = 1 never gives 0 * -Inf.
  power_term <- if (beta == 1) 0 * log_x else (beta - 1) * log_x

  log_alpha_g <- log(alpha + cdf_base)
  log_tail_factor <- log(alpha + 1 + surv_base)

  list(
    log_surv = 2 * log(alpha) - eta - 2 * log_alpha_g,
    log_dens = 2 * log(alpha) + log(beta) + log(lambda) + power_term -
      eta + log_tail_factor - 3 * log_alpha_g,
    log_haz = log(beta) + log(lambda) + power_term +
      log_tail_factor - log_alpha_g
  )
}


#' Log-Hazard of the NME-Weibull Distribution at Zero and Infinity
#'
#' Internal helper returning the limits in Equation (23) of the thesis:
#' at \eqn{x = 0} the hazard is infinite, \eqn{\lambda(\alpha + 2)/\alpha},
#' or zero for \eqn{\beta < 1}, \eqn{\beta = 1}, or \eqn{\beta > 1}; as
#' \eqn{x \to \infty} it is zero, \eqn{\lambda}, or infinite.
#'
#' @param at Either `"zero"` or `"infinity"`.
#' @param alpha,beta,lambda Parameters of the NME-Weibull distribution.
#'
#' @return A single log-hazard value.
#' @noRd
.nmeweibull_log_hazard_limit <- function(at, alpha, beta, lambda) {
  if (at == "zero") {
    if (beta < 1) return(Inf)
    if (beta > 1) return(-Inf)
    return(log(lambda) + log(alpha + 2) - log(alpha))
  }

  if (beta < 1) return(-Inf)
  if (beta > 1) return(Inf)
  log(lambda)
}


#' Density of the NME-Weibull Distribution
#'
#' Computes the probability density function (PDF) of the New Modified
#' Exponential-Weibull (NME-Weibull) distribution,
#' \deqn{f(x) = \frac{\alpha^2 \beta \lambda x^{\beta - 1} e^{-\lambda
#'   x^\beta} (\alpha + 1 + e^{-\lambda x^\beta})}{(\alpha + 1 -
#'   e^{-\lambda x^\beta})^3}, \quad x \ge 0,}
#' with \eqn{\alpha > 0}, shape \eqn{\beta > 0}, and rate
#' \eqn{\lambda > 0} (Alshanbari et al., 2023).
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Additional parameter of the NME-X family (\eqn{\alpha > 0}).
#'   The distribution tends to the Weibull distribution as
#'   \eqn{\alpha \to \infty}.
#' @param beta Shape parameter of the baseline Weibull distribution
#'   (\eqn{\beta > 0}).
#' @param lambda Rate parameter of the baseline Weibull distribution
#'   (\eqn{\lambda > 0}), so that the baseline CDF is
#'   \eqn{G(x) = 1 - e^{-\lambda x^\beta}}.
#' @param log Logical; if `TRUE`, the log-density is returned.
#'
#' @return A numeric vector of the same length as `x`. The density is zero
#'   for `x < 0`. `NA` values in `x` are propagated.
#'
#' @details The parameterisation of the baseline Weibull distribution is
#'   related to [stats::dweibull()] through
#'   \eqn{\lambda = \sigma^{-\beta}}, where \eqn{\sigma} is the `scale`
#'   argument of [stats::dweibull()].
#'
#' @references Alshanbari, H. M., Odhah, O. H., Ahmad, Z., Khan, F., and
#'   El-Bagoury, A. A.-A. H. (2023). A new probability distribution: Model,
#'   theory and analyzing the recovery time data. *Axioms*, 12(5), 477.
#'   \doi{10.3390/axioms12050477}
#'
#' @seealso [pnmeweibull()], [snmeweibull()], [hnmeweibull()],
#'   [qnmeweibull()], [rnmeweibull()].
#' @export
#'
#' @examples
#' dnmeweibull(1, alpha = 1, beta = 2, lambda = 0.5)
#' dnmeweibull(c(0.5, 1, 2), alpha = 0.5, beta = 1, lambda = 0.8, log = TRUE)
#'
#' # The density integrates to one
#' integrate(dnmeweibull, 0, Inf, alpha = 1, beta = 2, lambda = 0.5)
dnmeweibull <- function(x, alpha, beta, lambda, log = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_numeric_arg(x, "x")

  out <- rep(NA_real_, length(x))
  out[!is.na(x) & (x < 0 | x == Inf)] <- -Inf

  at_zero <- !is.na(x) & x == 0
  if (any(at_zero)) {
    out[at_zero] <- if (beta < 1) {
      Inf
    } else if (beta > 1) {
      -Inf
    } else {
      log(lambda) + log(alpha + 2) - log(alpha)
    }
  }

  inside <- !is.na(x) & x > 0 & is.finite(x)
  if (any(inside)) {
    out[inside] <- .nmeweibull_core(x[inside], alpha, beta, lambda)$log_dens
  }

  if (log) out else exp(out)
}


#' Distribution Function of the NME-Weibull Distribution
#'
#' Computes the cumulative distribution function (CDF) of the NME-Weibull
#' distribution,
#' \deqn{F(x) = 1 - \frac{\alpha^2 e^{-\lambda x^\beta}}{(\alpha + 1 -
#'   e^{-\lambda x^\beta})^2}, \quad x \ge 0.}
#'
#' @inheritParams dnmeweibull
#' @param q Numeric vector of quantiles.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise \eqn{P(X > q)}.
#' @param log.p Logical; if `TRUE`, probabilities are returned on the log
#'   scale.
#'
#' @return A numeric vector of probabilities of the same length as `q`.
#'
#' @details The upper tail (`lower.tail = FALSE`) is computed directly from
#'   the survival function rather than as `1 - F(q)`, so small tail
#'   probabilities keep full relative accuracy.
#'
#' @seealso [dnmeweibull()], [snmeweibull()], [qnmeweibull()].
#' @export
#'
#' @examples
#' pnmeweibull(c(0, 0.5, 1, 2), alpha = 1, beta = 2, lambda = 0.5)
#' pnmeweibull(5, alpha = 1, beta = 2, lambda = 0.5, lower.tail = FALSE)
pnmeweibull <- function(q, alpha, beta, lambda,
                        lower.tail = TRUE, log.p = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_numeric_arg(q, "q")

  log_surv <- .nmeweibull_log_survival(q, alpha, beta, lambda)

  if (!lower.tail) {
    return(if (log.p) log_surv else exp(log_surv))
  }

  cdf <- -expm1(log_surv)
  if (log.p) log(cdf) else cdf
}


#' Log-Survival at Arbitrary Points
#'
#' @inheritParams pnmeweibull
#' @return Log-survival values with the conventions of [snmeweibull()].
#' @noRd
.nmeweibull_log_survival <- function(q, alpha, beta, lambda) {
  out <- rep(NA_real_, length(q))
  out[!is.na(q) & q <= 0] <- 0
  out[!is.na(q) & q == Inf] <- -Inf

  inside <- !is.na(q) & q > 0 & is.finite(q)
  if (any(inside)) {
    out[inside] <- .nmeweibull_core(q[inside], alpha, beta, lambda)$log_surv
  }

  out
}


#' Survival Function of the NME-Weibull Distribution
#'
#' Computes the survival function of the NME-Weibull distribution,
#' \deqn{S(x) = P(X > x) = \frac{\alpha^2 e^{-\lambda x^\beta}}{(\alpha + 1 -
#'   e^{-\lambda x^\beta})^2}, \quad x \ge 0.}
#'
#' @inheritParams dnmeweibull
#' @param log Logical; if `TRUE`, the log-survival is returned.
#'
#' @return A numeric vector of survival probabilities of the same length as
#'   `x`. The survival function equals one for `x <= 0`.
#'
#' @seealso [pnmeweibull()], [hnmeweibull()].
#' @export
#'
#' @examples
#' snmeweibull(c(0, 0.5, 1, 2), alpha = 1, beta = 2, lambda = 0.5)
#'
#' # S(x) = 1 - F(x)
#' x <- c(0.5, 1, 2)
#' all.equal(snmeweibull(x, 1, 2, 0.5), 1 - pnmeweibull(x, 1, 2, 0.5))
snmeweibull <- function(x, alpha, beta, lambda, log = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_numeric_arg(x, "x")

  log_surv <- .nmeweibull_log_survival(x, alpha, beta, lambda)
  if (log) log_surv else exp(log_surv)
}


#' Hazard Function of the NME-Weibull Distribution
#'
#' Computes the hazard function of the NME-Weibull distribution,
#' \deqn{h(x) = \frac{f(x)}{S(x)} = \beta \lambda x^{\beta - 1}
#'   \frac{\alpha + 1 + e^{-\lambda x^\beta}}{\alpha + 1 -
#'   e^{-\lambda x^\beta}}, \quad x \ge 0.}
#'
#' @inheritParams dnmeweibull
#' @param log Logical; if `TRUE`, the log-hazard is returned.
#'
#' @return A numeric vector of hazard values of the same length as `x`. The
#'   hazard is zero for `x < 0`. At `x = 0` and `x = Inf` the limits of the
#'   hazard function are returned: for `x = 0` these are `Inf`,
#'   \eqn{\lambda(\alpha + 2)/\alpha}, or `0`, and for `x = Inf` they are
#'   `0`, \eqn{\lambda}, or `Inf`, for \eqn{\beta < 1}, \eqn{\beta = 1}, or
#'   \eqn{\beta > 1}, respectively.
#'
#' @details The hazard is the Weibull hazard \eqn{\beta \lambda x^{\beta -
#'   1}} multiplied by a factor that decreases from \eqn{(\alpha +
#'   2)/\alpha} to one. Depending on the parameters, the hazard can be
#'   decreasing, increasing, unimodal, or modified unimodal.
#'
#' @seealso [dnmeweibull()], [snmeweibull()].
#' @export
#'
#' @examples
#' x <- seq(0.1, 3, by = 0.1)
#' plot(x, hnmeweibull(x, alpha = 0.05, beta = 2, lambda = 1), type = "l",
#'      ylab = "h(x)")
#'
#' # h(x) = f(x) / S(x)
#' all.equal(hnmeweibull(1, 1, 2, 0.5),
#'           dnmeweibull(1, 1, 2, 0.5) / snmeweibull(1, 1, 2, 0.5))
hnmeweibull <- function(x, alpha, beta, lambda, log = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_numeric_arg(x, "x")

  out <- rep(NA_real_, length(x))
  out[!is.na(x) & x < 0] <- -Inf

  at_zero <- !is.na(x) & x == 0
  if (any(at_zero)) {
    out[at_zero] <- .nmeweibull_log_hazard_limit("zero", alpha, beta, lambda)
  }

  at_inf <- !is.na(x) & x == Inf
  if (any(at_inf)) {
    out[at_inf] <- .nmeweibull_log_hazard_limit("infinity", alpha, beta,
                                                lambda)
  }

  inside <- !is.na(x) & x > 0 & is.finite(x)
  if (any(inside)) {
    out[inside] <- .nmeweibull_core(x[inside], alpha, beta, lambda)$log_haz
  }

  if (log) out else exp(out)
}


#' Quantile Function of the NME-Weibull Distribution
#'
#' Computes the quantile function of the NME-Weibull distribution in closed
#' form. With \eqn{s = 1 - p},
#' \deqn{w = \frac{2\sqrt{s}(\alpha + 1)}{\alpha + \sqrt{\alpha^2 +
#'   4s(\alpha + 1)}}, \qquad Q(p) = \left[\frac{-2 \ln w}{\lambda}
#'   \right]^{1/\beta}.}
#'
#' @inheritParams dnmeweibull
#' @param p Numeric vector of probabilities in \eqn{[0, 1]}.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#' @param log.p Logical; if `TRUE`, probabilities `p` are given as
#'   `log(p)`.
#'
#' @return A numeric vector of quantiles of the same length as `p`. The
#'   quantile is `0` for a lower-tail probability of `0` and `Inf` for a
#'   lower-tail probability of `1`.
#'
#' @details The closed form follows from solving \eqn{S(x) = s} for
#'   \eqn{w = \sqrt{e^{-\lambda x^\beta}}}, which gives the quadratic
#'   equation \eqn{\sqrt{s}\,w^2 + \alpha w - \sqrt{s}(\alpha + 1) = 0}.
#'   The positive root is written in rationalised form to avoid
#'   cancellation when \eqn{s} is small. No iteration is needed.
#'
#' @seealso [pnmeweibull()], [rnmeweibull()], [varnmeweibull()].
#' @export
#'
#' @examples
#' qnmeweibull(c(0.25, 0.5, 0.75), alpha = 1, beta = 2, lambda = 0.5)
#'
#' # Q(F(x)) = x
#' qnmeweibull(pnmeweibull(1.3, 1, 2, 0.5), 1, 2, 0.5)
qnmeweibull <- function(p, alpha, beta, lambda,
                        lower.tail = TRUE, log.p = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_numeric_arg(p, "p")

  if (log.p) {
    p <- exp(p)
  }

  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("p must be between 0 and 1.", call. = FALSE)
  }

  # s is the upper-tail (survival) probability.
  s <- if (lower.tail) 1 - p else p

  out <- rep(NA_real_, length(p))
  out[!is.na(s) & s == 1] <- 0
  out[!is.na(s) & s == 0] <- Inf

  inside <- !is.na(s) & s > 0 & s < 1
  if (any(inside)) {
    s_in <- s[inside]
    w <- 2 * sqrt(s_in) * (alpha + 1) /
      (alpha + sqrt(alpha^2 + 4 * s_in * (alpha + 1)))
    out[inside] <- (-2 * log(w) / lambda)^(1 / beta)
  }

  out
}


#' Random Generation from the NME-Weibull Distribution
#'
#' Generates random observations from the NME-Weibull distribution with the
#' inverse transform method, \eqn{X = Q(U)} with
#' \eqn{U \sim \mathrm{Uniform}(0, 1)}.
#'
#' @inheritParams dnmeweibull
#' @param n Number of observations. If `length(n) > 1`, the length is taken
#'   to be the number required, as in [stats::rnorm()].
#'
#' @return A numeric vector of `n` random observations.
#'
#' @details Results are reproducible for a fixed seed set with
#'   [set.seed()].
#'
#' @seealso [qnmeweibull()].
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnmeweibull(1000, alpha = 1, beta = 2, lambda = 0.5)
#' summary(x)
rnmeweibull <- function(n, alpha, beta, lambda) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (length(n) > 1L) {
    n <- length(n)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 0) {
    stop("n must be a single non-negative finite number.", call. = FALSE)
  }

  n <- as.integer(n)
  if (n == 0L) {
    return(numeric(0))
  }

  qnmeweibull(stats::runif(n), alpha = alpha, beta = beta, lambda = lambda)
}
