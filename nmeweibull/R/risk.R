#' Check Probability Levels of a Risk Measure
#'
#' @param p Numeric vector of probability levels.
#' @param name Argument name used in error messages.
#'
#' @return `NULL`, invisibly.
#' @noRd
.check_risk_prob <- function(p, name) {
  if (!is.numeric(p)) {
    stop(name, " must be numeric.", call. = FALSE)
  }

  if (any(p <= 0 | p >= 1, na.rm = TRUE)) {
    stop(name, " must be strictly between 0 and 1.", call. = FALSE)
  }

  invisible(NULL)
}


#' Average of the NME-Weibull Quantile Function over an Interval
#'
#' Internal helper computing
#' \eqn{\frac{1}{b - a}\int_a^b Q(u)\,du} with [stats::integrate()]. The
#' integral is evaluated in the upper-tail probability \eqn{s = 1 - u}, which
#' keeps full accuracy close to \eqn{u = 1}.
#'
#' @param a,b Lower and upper probability levels, \eqn{0 < a < b \le 1}.
#' @param alpha,beta,lambda Parameters of the NME-Weibull distribution.
#' @param rel.tol,subdivisions Passed to [stats::integrate()].
#'
#' @return A single numeric value. A warning is issued if
#'   [stats::integrate()] does not report success.
#' @noRd
.mean_quantile <- function(a, b, alpha, beta, lambda, rel.tol, subdivisions) {
  integrand <- function(s) {
    qnmeweibull(s, alpha, beta, lambda, lower.tail = FALSE)
  }

  result <- stats::integrate(
    integrand,
    lower = 1 - b,
    upper = 1 - a,
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )

  if (!identical(result$message, "OK")) {
    warning("Numerical integration did not report success for the interval (",
            a, ", ", b, "): ", result$message, call. = FALSE)
  }

  result$value / (b - a)
}


#' Average of the NME-Weibull Quantile Function above a Level
#'
#' Internal helper computing \eqn{\frac{1}{1 - p}\int_p^1 Q(u)\,du}. The
#' quantile function is unbounded as \eqn{u \to 1}, so the substitution
#' \eqn{1 - u = (1 - p)e^{-t}} is used, which gives
#' \deqn{\frac{1}{1 - p}\int_p^1 Q(u)\,du = \int_0^\infty
#'   Q_{\mathrm{upper}}\{(1 - p)e^{-t}\}\, e^{-t}\,dt,}
#' a smooth integrand that decays exponentially, where
#' \eqn{Q_{\mathrm{upper}}(s)} is the quantile at upper-tail probability
#' \eqn{s}.
#'
#' @param p Probability level, \eqn{0 < p < 1}.
#' @inheritParams .mean_quantile
#'
#' @return A single numeric value. A warning is issued if
#'   [stats::integrate()] does not report success.
#' @noRd
.tail_mean_quantile <- function(p, alpha, beta, lambda, rel.tol,
                                subdivisions) {
  integrand <- function(t) {
    s <- (1 - p) * exp(-t)
    out <- numeric(length(t))
    positive <- s > 0
    out[positive] <- qnmeweibull(s[positive], alpha, beta, lambda,
                                 lower.tail = FALSE) * exp(-t[positive])
    out
  }

  result <- stats::integrate(
    integrand,
    lower = 0,
    upper = Inf,
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )

  if (!identical(result$message, "OK")) {
    warning("Numerical integration did not report success for p = ", p,
            ": ", result$message, call. = FALSE)
  }

  result$value
}


#' Value-at-Risk of the NME-Weibull Distribution
#'
#' Computes the Value-at-Risk of the NME-Weibull distribution at probability
#' level \eqn{p},
#' \deqn{\mathrm{VaR}_p(X) = \inf\{x : F(x) \ge p\} = Q(p).}
#' For survival-time data, \eqn{\mathrm{VaR}_p} is the time not exceeded by
#' \eqn{100p} percent of subjects; \eqn{\mathrm{VaR}_{0.5}} is the median.
#'
#' @inheritParams dnmeweibull
#' @param p Numeric vector of probability levels, strictly between 0 and 1.
#'
#' @return A numeric vector of VaR values of the same length as `p`.
#'
#' @seealso [tvarnmeweibull()], [rvarnmeweibull()], [qnmeweibull()].
#' @export
#'
#' @examples
#' varnmeweibull(c(0.90, 0.95, 0.99), alpha = 1, beta = 2, lambda = 0.5)
varnmeweibull <- function(p, alpha, beta, lambda) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_risk_prob(p, "p")

  qnmeweibull(p, alpha = alpha, beta = beta, lambda = lambda)
}


#' Tail Value-at-Risk of the NME-Weibull Distribution
#'
#' Computes the Tail Value-at-Risk of the NME-Weibull distribution at
#' probability level \eqn{p},
#' \deqn{\mathrm{TVaR}_p(X) = \frac{1}{1 - p}\int_p^1 Q(u)\,du =
#'   E[X \mid X > \mathrm{VaR}_p(X)],}
#' i.e. the mean of the longest \eqn{100(1 - p)} percent of survival times.
#'
#' @inheritParams varnmeweibull
#' @param rel.tol Relative accuracy requested from [stats::integrate()].
#' @param subdivisions Maximum number of subintervals used by
#'   [stats::integrate()].
#'
#' @return A numeric vector of TVaR values of the same length as `p`. A
#'   warning is issued if the numerical integration does not report success.
#'
#' @details The integral is finite because the NME-Weibull survival function
#'   is bounded by the Weibull survival function \eqn{e^{-\lambda x^\beta}}.
#'   Because \eqn{Q(u)} is unbounded as \eqn{u \to 1}, the integral is
#'   evaluated with [stats::integrate()] after the substitution
#'   \eqn{1 - u = (1 - p)e^{-t}}, which turns it into
#'   \eqn{\int_0^\infty Q(1 - (1 - p)e^{-t})\,e^{-t}\,dt} with a smooth,
#'   exponentially decaying integrand. TVaR is always at least
#'   \eqn{\mathrm{VaR}_p}.
#'
#' @seealso [varnmeweibull()], [rvarnmeweibull()].
#' @export
#'
#' @examples
#' tvarnmeweibull(c(0.90, 0.95, 0.99), alpha = 1, beta = 2, lambda = 0.5)
tvarnmeweibull <- function(p, alpha, beta, lambda,
                           rel.tol = .Machine$double.eps^0.25,
                           subdivisions = 100L) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_risk_prob(p, "p")

  vapply(p, function(level) {
    if (is.na(level)) {
      return(NA_real_)
    }
    .tail_mean_quantile(level, alpha, beta, lambda, rel.tol, subdivisions)
  }, numeric(1))
}


#' Range Value-at-Risk of the NME-Weibull Distribution
#'
#' Computes the Range Value-at-Risk of the NME-Weibull distribution over the
#' probability interval \eqn{(p_1, p_2)},
#' \deqn{\mathrm{RVaR}_{p_1, p_2}(X) = \frac{1}{p_2 - p_1}
#'   \int_{p_1}^{p_2} Q(u)\,du,}
#' the mean of the quantiles between the two levels. RVaR lies between
#' \eqn{\mathrm{VaR}_{p_1}} and \eqn{\mathrm{VaR}_{p_2}}, and is less
#' sensitive to extreme values than TVaR.
#'
#' @inheritParams tvarnmeweibull
#' @param p_lower Numeric vector of lower probability levels \eqn{p_1}.
#' @param p_upper Numeric vector of upper probability levels \eqn{p_2}.
#'   `p_lower` and `p_upper` must have the same length, or one of them must
#'   have length one.
#'
#' @return A numeric vector of RVaR values. A warning is issued if the
#'   numerical integration does not report success.
#'
#' @seealso [varnmeweibull()], [tvarnmeweibull()].
#' @export
#'
#' @examples
#' rvarnmeweibull(p_lower = c(0.90, 0.95), p_upper = c(0.95, 0.99),
#'                alpha = 1, beta = 2, lambda = 0.5)
rvarnmeweibull <- function(p_lower, p_upper, alpha, beta, lambda,
                           rel.tol = .Machine$double.eps^0.25,
                           subdivisions = 100L) {
  .check_nmeweibull_params(alpha, beta, lambda)
  .check_risk_prob(p_lower, "p_lower")
  .check_risk_prob(p_upper, "p_upper")

  if (length(p_lower) != length(p_upper)) {
    if (length(p_lower) == 1L) {
      p_lower <- rep(p_lower, length(p_upper))
    } else if (length(p_upper) == 1L) {
      p_upper <- rep(p_upper, length(p_lower))
    } else {
      stop("p_lower and p_upper must have the same length, or one of them ",
           "must have length 1.", call. = FALSE)
    }
  }

  if (any(p_lower >= p_upper, na.rm = TRUE)) {
    stop("p_lower must be smaller than p_upper.", call. = FALSE)
  }

  vapply(seq_along(p_lower), function(i) {
    if (is.na(p_lower[i]) || is.na(p_upper[i])) {
      return(NA_real_)
    }
    .mean_quantile(p_lower[i], p_upper[i], alpha, beta, lambda, rel.tol,
                   subdivisions)
  }, numeric(1))
}
