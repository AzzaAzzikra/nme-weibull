#' Value-at-Risk of the NME-Weibull Distribution
#'
#' Computes the Value-at-Risk (VaR) of the New Modified Exponential-Weibull
#' distribution at a given probability level.
#'
#' @param p Numeric vector of probability levels.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param lower.tail Logical; if TRUE, probabilities are P(X <= x), otherwise P(X > x).
#' @param log.p Logical; if TRUE, probabilities are given on the log scale.
#'
#' @return A numeric vector of VaR values.
#' @export
#'
#' @examples
#' varnmeweibull(p = 0.95, alpha = 1, beta = 2, lambda = 0.5)
#' varnmeweibull(p = c(0.90, 0.95, 0.99), alpha = 1, beta = 2, lambda = 0.5)
varnmeweibull <- function(p, alpha, beta, lambda,
                          lower.tail = TRUE, log.p = FALSE) {
  qnmeweibull(
    p = p,
    alpha = alpha,
    beta = beta,
    lambda = lambda,
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' Tail Value-at-Risk of the NME-Weibull Distribution
#'
#' Computes the Tail Value-at-Risk (TVaR) of the New Modified
#' Exponential-Weibull distribution at a given probability level.
#'
#' @param p Numeric vector of probability levels.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param rel.tol Relative accuracy requested in numerical integration.
#' @param subdivisions Maximum number of subintervals used in numerical integration.
#'
#' @return A numeric vector of TVaR values.
#' @export
#'
#' @examples
#' tvarnmeweibull(p = 0.95, alpha = 1, beta = 2, lambda = 0.5)
#' tvarnmeweibull(p = c(0.90, 0.95, 0.99), alpha = 1, beta = 2, lambda = 0.5)
tvarnmeweibull <- function(p, alpha, beta, lambda,
                           rel.tol = .Machine$double.eps^0.25,
                           subdivisions = 100L) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(p)) {
    stop("p must be numeric.", call. = FALSE)
  }

  if (any(p <= 0 | p >= 1, na.rm = TRUE)) {
    stop("p must be between 0 and 1, exclusive.", call. = FALSE)
  }

  result <- rep(NA_real_, length(p))

  for (i in seq_along(p)) {
    if (is.na(p[i])) {
      next
    }

    integrand <- function(u) {
      qnmeweibull(
        p = u,
        alpha = alpha,
        beta = beta,
        lambda = lambda
      )
    }

    integral_value <- stats::integrate(
      f = integrand,
      lower = p[i],
      upper = 1,
      rel.tol = rel.tol,
      subdivisions = subdivisions,
      stop.on.error = FALSE
    )

    if (!isTRUE(integral_value$message == "OK")) {
      warning(
        paste0("Numerical integration may not have converged for p = ", p[i], "."),
        call. = FALSE
      )
    }

    result[i] <- integral_value$value / (1 - p[i])
  }

  result
}

#' Range Value-at-Risk of the NME-Weibull Distribution
#'
#' Computes the Range Value-at-Risk (RVaR) of the New Modified
#' Exponential-Weibull distribution over a probability interval.
#'
#' @param p_lower Numeric vector of lower probability levels.
#' @param p_upper Numeric vector of upper probability levels.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param rel.tol Relative accuracy requested in numerical integration.
#' @param subdivisions Maximum number of subintervals used in numerical integration.
#'
#' @return A numeric vector of RVaR values.
#' @export
#'
#' @examples
#' rvarnmeweibull(p_lower = 0.90, p_upper = 0.95,
#'                alpha = 1, beta = 2, lambda = 0.5)
#' rvarnmeweibull(p_lower = c(0.90, 0.95), p_upper = c(0.95, 0.99),
#'                alpha = 1, beta = 2, lambda = 0.5)
rvarnmeweibull <- function(p_lower, p_upper, alpha, beta, lambda,
                           rel.tol = .Machine$double.eps^0.25,
                           subdivisions = 100L) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(p_lower) || !is.numeric(p_upper)) {
    stop("p_lower and p_upper must be numeric.", call. = FALSE)
  }

  if (length(p_lower) != length(p_upper)) {
    if (length(p_lower) == 1) {
      p_lower <- rep(p_lower, length(p_upper))
    } else if (length(p_upper) == 1) {
      p_upper <- rep(p_upper, length(p_lower))
    } else {
      stop("p_lower and p_upper must have the same length, or one of them must have length 1.",
           call. = FALSE)
    }
  }

  if (any(p_lower <= 0 | p_lower >= 1, na.rm = TRUE) ||
      any(p_upper <= 0 | p_upper >= 1, na.rm = TRUE)) {
    stop("p_lower and p_upper must be between 0 and 1, exclusive.",
         call. = FALSE)
  }

  if (any(p_lower >= p_upper, na.rm = TRUE)) {
    stop("p_lower must be smaller than p_upper.", call. = FALSE)
  }

  result <- rep(NA_real_, length(p_lower))

  for (i in seq_along(p_lower)) {
    if (is.na(p_lower[i]) || is.na(p_upper[i])) {
      next
    }

    integrand <- function(u) {
      qnmeweibull(
        p = u,
        alpha = alpha,
        beta = beta,
        lambda = lambda
      )
    }

    integral_value <- stats::integrate(
      f = integrand,
      lower = p_lower[i],
      upper = p_upper[i],
      rel.tol = rel.tol,
      subdivisions = subdivisions,
      stop.on.error = FALSE
    )

    if (!isTRUE(integral_value$message == "OK")) {
      warning(
        paste0(
          "Numerical integration may not have converged for interval [",
          p_lower[i], ", ", p_upper[i], "]."
        ),
        call. = FALSE
      )
    }

    result[i] <- integral_value$value / (p_upper[i] - p_lower[i])
  }

  result
}
