#' Check NME-Weibull Parameters
#'
#' Internal function to check whether all NME-Weibull parameters are positive.
#'
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#'
#' @noRd
.check_nmeweibull_params <- function(alpha, beta, lambda) {
  params <- c(alpha = alpha, beta = beta, lambda = lambda)

  if (any(!is.finite(params))) {
    stop("All parameters must be finite.", call. = FALSE)
  }

  if (any(params <= 0)) {
    stop("All parameters must be positive.", call. = FALSE)
  }
}


#' Density Function of the NME-Weibull Distribution
#'
#' Computes the probability density function of the New Modified
#' Exponential-Weibull distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param log Logical; if TRUE, probabilities are returned on the log scale.
#'
#' @return A numeric vector of density values.
#' @export
#'
#' @examples
#' dnmeweibull(x = 1, alpha = 1, beta = 1, lambda = 1)
#' dnmeweibull(x = c(0.5, 1, 2), alpha = 1, beta = 2, lambda = 0.5)
dnmeweibull <- function(x, alpha, beta, lambda, log = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }

  density <- rep(NA_real_, length(x))

  density[x < 0 & !is.na(x)] <- 0
  density[is.infinite(x) & x > 0] <- 0

  valid <- x >= 0 & is.finite(x) & !is.na(x)

  exp_term <- exp(-lambda * x[valid]^beta)
  weibull_cdf <- 1 - exp_term

  density[valid] <-
    (alpha^2 * beta * lambda * x[valid]^(beta - 1) * exp_term *
       (alpha + 2 - weibull_cdf)) /
    (alpha + weibull_cdf)^3

  if (log) {
    return(log(density))
  }

  density
}

#' Distribution Function of the NME-Weibull Distribution
#'
#' Computes the cumulative distribution function of the New Modified
#' Exponential-Weibull distribution.
#'
#' @param q Numeric vector of quantiles.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param lower.tail Logical; if TRUE, probabilities are P(X <= q), otherwise P(X > q).
#' @param log.p Logical; if TRUE, probabilities are returned on the log scale.
#'
#' @return A numeric vector of cumulative probabilities.
#' @export
#'
#' @examples
#' pnmeweibull(q = 1, alpha = 1, beta = 1, lambda = 1)
#' pnmeweibull(q = c(0, 0.5, 1, 2), alpha = 1, beta = 2, lambda = 0.5)
pnmeweibull <- function(q, alpha, beta, lambda,
                        lower.tail = TRUE, log.p = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(q)) {
    stop("q must be numeric.", call. = FALSE)
  }

  cdf <- rep(NA_real_, length(q))

  cdf[q < 0 & !is.na(q)] <- 0
  cdf[is.infinite(q) & q > 0] <- 1

  valid <- q >= 0 & is.finite(q) & !is.na(q)

  exp_term <- exp(-lambda * q[valid]^beta)
  weibull_cdf <- 1 - exp_term

  cdf[valid] <- 1 -
    (alpha^2 * exp_term) / (alpha + weibull_cdf)^2

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  if (log.p) {
    return(log(cdf))
  }

  cdf
}

#' Survival Function of the NME-Weibull Distribution
#'
#' Computes the survival function of the New Modified Exponential-Weibull
#' distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param log Logical; if TRUE, probabilities are returned on the log scale.
#'
#' @return A numeric vector of survival probabilities.
#' @export
#'
#' @examples
#' snmeweibull(x = 1, alpha = 1, beta = 1, lambda = 1)
#' snmeweibull(x = c(0, 0.5, 1, 2), alpha = 1, beta = 2, lambda = 0.5)
snmeweibull <- function(x, alpha, beta, lambda, log = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }

  survival <- rep(NA_real_, length(x))

  survival[x < 0 & !is.na(x)] <- 1
  survival[is.infinite(x) & x > 0] <- 0

  valid <- x >= 0 & is.finite(x) & !is.na(x)

  exp_term <- exp(-lambda * x[valid]^beta)
  weibull_cdf <- 1 - exp_term

  survival[valid] <-
    (alpha^2 * exp_term) / (alpha + weibull_cdf)^2

  if (log) {
    return(log(survival))
  }

  survival
}


#' Hazard Function of the NME-Weibull Distribution
#'
#' Computes the hazard function of the New Modified Exponential-Weibull
#' distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param log Logical; if TRUE, hazard values are returned on the log scale.
#'
#' @return A numeric vector of hazard values.
#' @export
#'
#' @examples
#' hnmeweibull(x = 1, alpha = 1, beta = 1, lambda = 1)
#' hnmeweibull(x = c(0.5, 1, 2), alpha = 1, beta = 2, lambda = 0.5)
hnmeweibull <- function(x, alpha, beta, lambda, log = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }

  hazard <- rep(NA_real_, length(x))

  hazard[x < 0 & !is.na(x)] <- NA_real_

  valid <- x >= 0 & is.finite(x) & !is.na(x)

  exp_term <- exp(-lambda * x[valid]^beta)
  weibull_cdf <- 1 - exp_term

  hazard[valid] <-
    beta * lambda * x[valid]^(beta - 1) *
    (alpha + 2 - weibull_cdf) /
    (alpha + weibull_cdf)

  hazard[is.infinite(x) & x > 0] <- NA_real_

  if (log) {
    return(log(hazard))
  }

  hazard
}

#' Quantile Function of the NME-Weibull Distribution
#'
#' Computes the quantile function of the New Modified Exponential-Weibull
#' distribution.
#'
#' @param p Numeric vector of probabilities.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#' @param lower.tail Logical; if TRUE, probabilities are P(X <= x), otherwise P(X > x).
#' @param log.p Logical; if TRUE, probabilities are given on the log scale.
#'
#' @return A numeric vector of quantiles.
#' @export
#'
#' @examples
#' qnmeweibull(p = 0.5, alpha = 1, beta = 1, lambda = 1)
#' qnmeweibull(p = c(0.25, 0.5, 0.75), alpha = 1, beta = 2, lambda = 0.5)
qnmeweibull <- function(p, alpha, beta, lambda,
                        lower.tail = TRUE, log.p = FALSE) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (!is.numeric(p)) {
    stop("p must be numeric.", call. = FALSE)
  }

  if (log.p) {
    p <- exp(p)
  }

  if (!lower.tail) {
    p <- 1 - p
  }

  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("p must be between 0 and 1.", call. = FALSE)
  }

  quantile <- rep(NA_real_, length(p))

  quantile[p == 0] <- 0
  quantile[p == 1] <- Inf

  valid <- p > 0 & p < 1 & !is.na(p)

  survival <- 1 - p[valid]
  a <- alpha + 1

  sqrt_z <- (-alpha + sqrt(alpha^2 + 4 * a * survival)) /
    (2 * sqrt(survival))

  z <- sqrt_z^2

  quantile[valid] <- (-log(z) / lambda)^(1 / beta)

  quantile
}

#' Random Generation from the NME-Weibull Distribution
#'
#' Generates random observations from the New Modified Exponential-Weibull
#' distribution using the inverse transform method.
#'
#' @param n Number of observations. If length(n) > 1, the length is taken to be
#'   the number required.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#'
#' @return A numeric vector of random observations.
#' @export
#'
#' @examples
#' rnmeweibull(n = 5, alpha = 1, beta = 1, lambda = 1)
#' rnmeweibull(n = 10, alpha = 1, beta = 2, lambda = 0.5)
rnmeweibull <- function(n, alpha, beta, lambda) {
  .check_nmeweibull_params(alpha, beta, lambda)

  if (length(n) > 1) {
    n <- length(n)
  }

  if (!is.numeric(n) || length(n) != 1 || is.na(n) || !is.finite(n)) {
    stop("n must be a non-negative finite number.", call. = FALSE)
  }

  n <- as.integer(n)

  if (n < 0) {
    stop("n must be non-negative.", call. = FALSE)
  }

  if (n == 0) {
    return(numeric(0))
  }

  u <- stats::runif(n)

  qnmeweibull(p = u, alpha = alpha, beta = beta, lambda = lambda)
}
