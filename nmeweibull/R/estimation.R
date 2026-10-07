#' Log-Density Used in the NME-Weibull Log-Likelihood
#'
#' Internal helper that evaluates the log-density of each observation, i.e.
#' the summands of the log-likelihood
#' \deqn{\ell(\theta) = 2n\ln\alpha + n\ln\beta + n\ln\lambda + (\beta - 1)
#'   \sum \ln x_i - \lambda \sum x_i^\beta + \sum \ln(\alpha + 2 - G_i) -
#'   3 \sum \ln(\alpha + G_i),}
#' with \eqn{G_i = 1 - e^{-\lambda x_i^\beta}} and
#' \eqn{x_i^\beta = \exp(\beta \ln x_i)}.
#'
#' @param x Numeric vector of strictly positive, finite observations.
#' @param alpha,beta,lambda Parameters of the NME-Weibull distribution.
#'
#' @return A numeric vector of log-density values.
#' @noRd
.log_dnmeweibull <- function(x, alpha, beta, lambda) {
  .nmeweibull_core(x, alpha, beta, lambda)$log_dens
}


#' Starting Values for NME-Weibull Maximum Likelihood Estimation
#'
#' Internal helper that builds the 75 starting points used by
#' [fitnmeweibull()]: every combination of
#' \eqn{\alpha_0 \in \{0.5, 1, 2\}}, \eqn{\beta_0 \in \{0.5, 1, 1.5, 2,
#' 3\}}, and \eqn{c \in \{0.25, 0.5, 1, 2, 4\}}, with
#' \eqn{\lambda_0 = c / \bar{x}^{\beta_0}} so that the starting value of
#' \eqn{\lambda_0 \bar{x}^{\beta_0}} is neither close to zero nor large.
#'
#' @param x Numeric vector of strictly positive observations.
#'
#' @return A 75-row matrix with columns `alpha`, `beta`, and `lambda`.
#' @noRd
.make_nmeweibull_starts <- function(x) {
  grid <- expand.grid(
    multiplier = c(0.25, 0.5, 1, 2, 4),
    beta = c(0.5, 1, 1.5, 2, 3),
    alpha = c(0.5, 1, 2)
  )

  cbind(
    alpha = grid$alpha,
    beta = grid$beta,
    lambda = grid$multiplier / mean(x)^grid$beta
  )
}


#' Check Bounds Supplied to fitnmeweibull()
#'
#' @param bounds Named numeric vector.
#' @param name Argument name used in error messages.
#'
#' @return `bounds` reordered as `alpha`, `beta`, `lambda`.
#' @noRd
.check_fit_bounds <- function(bounds, name) {
  required <- c("alpha", "beta", "lambda")

  if (!is.numeric(bounds) || !all(required %in% names(bounds))) {
    stop(name, " must be a named numeric vector with names alpha, beta, ",
         "and lambda.", call. = FALSE)
  }

  bounds <- bounds[required]

  if (anyNA(bounds) || any(bounds < 0)) {
    stop(name, " must be non-negative.", call. = FALSE)
  }

  bounds
}


#' Fit the NME-Weibull Distribution by Maximum Likelihood
#'
#' Estimates the parameters \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda} of
#' the NME-Weibull distribution from complete (uncensored) data by maximum
#' likelihood, and reports the maximised log-likelihood, AIC, BIC, and the
#' convergence code.
#'
#' @param x Numeric vector of positive observations. `NA`, `NaN`, and
#'   infinite values are removed before fitting.
#' @param start Optional named numeric vector with elements `alpha`, `beta`,
#'   and `lambda`. When supplied, it is used as an additional starting point
#'   on top of the 75 automatic ones.
#' @param lower,upper Named numeric vectors with elements `alpha`, `beta`,
#'   and `lambda` giving the parameter space searched. The default upper
#'   bound of `alpha` is \eqn{10^4}; because the Weibull distribution is the
#'   limit of the NME-Weibull distribution as \eqn{\alpha \to \infty}, an
#'   estimate of `alpha` at this bound indicates that the data cannot be
#'   distinguished from a Weibull sample. `lambda` is unbounded by default so
#'   that data measured in any unit can be fitted.
#' @param hessian Logical; if `TRUE`, the Hessian of the negative
#'   log-likelihood with respect to
#'   \eqn{(\ln\alpha, \ln\beta, \ln\lambda)} is computed at the estimate.
#' @param control A list of control parameters passed to [stats::optim()].
#'   The defaults are `maxit = 5000` and `reltol = 1e-10`.
#'
#' @return An object of class `"fitnmeweibull"`, a list with components
#'   \describe{
#'     \item{estimate}{Named vector of estimates of `alpha`, `beta`, and
#'       `lambda`.}
#'     \item{logLik}{Maximised log-likelihood.}
#'     \item{AIC, BIC}{Akaike and Bayesian information criteria with
#'       \eqn{k = 3} parameters.}
#'     \item{convergence}{Convergence code of the best run reported by
#'       [stats::optim()]; `0` means the algorithm stopped normally.}
#'     \item{message}{Message returned by [stats::optim()], or `NULL`.}
#'     \item{n}{Number of observations used.}
#'     \item{boundary}{Logical; `TRUE` if any estimate lies within 1\% of a
#'       finite bound in `lower` or `upper`.}
#'     \item{n_starts}{Number of starting points that produced a finite
#'       log-likelihood.}
#'     \item{hessian}{Hessian matrix if `hessian = TRUE`, otherwise
#'       `NULL`.}
#'     \item{optim}{The full [stats::optim()] result of the best run, on the
#'       internal scale described in Details.}
#'   }
#'
#' @details The log-likelihood has no closed-form maximiser, so it is
#'   maximised numerically. The parameters are estimated on the log scale,
#'   \eqn{\phi = (\ln\alpha, \ln\beta, \ln\lambda)}, which keeps them
#'   positive; by the invariance property of maximum likelihood,
#'   \eqn{\exp(\hat{\phi})} is the estimate of \eqn{\theta}. From each of 75
#'   starting points (see Equation (41) of the thesis) the Nelder-Mead
#'   algorithm is run first and its result is refined with BFGS, both via
#'   [stats::optim()]. The run with the smallest negative log-likelihood is
#'   kept.
#'
#'   Internally the data are divided by their geometric mean \eqn{m} and the
#'   rate is written as \eqn{\lambda = \lambda^* m^{-\beta}}. This changes
#'   neither the likelihood nor the estimates, but it keeps the optimisation
#'   well conditioned for data in large units such as days or currency.
#'   All reported quantities are on the original scale of `x`.
#'
#'   A convergence code of `0` only means that the algorithm stopped
#'   normally. Check `boundary` as well: an estimate of `alpha` at its upper
#'   bound means the fitted model is practically a Weibull distribution.
#'
#' @seealso [dnmeweibull()] for the density and [varnmeweibull()],
#'   [tvarnmeweibull()], and [rvarnmeweibull()] for risk measures based on
#'   the fitted parameters.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnmeweibull(200, alpha = 1, beta = 2, lambda = 0.5)
#' fit <- fitnmeweibull(x)
#' fit
#' fit$estimate
#'
#' # Fitted parameters can be passed to the other functions
#' est <- fit$estimate
#' varnmeweibull(0.95, est["alpha"], est["beta"], est["lambda"])
fitnmeweibull <- function(x,
                          start = NULL,
                          lower = c(alpha = 1e-8, beta = 1e-8, lambda = 0),
                          upper = c(alpha = 1e4, beta = 1e4, lambda = Inf),
                          hessian = FALSE,
                          control = list()) {
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }

  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    stop("x must contain at least one finite observation.", call. = FALSE)
  }

  if (any(x <= 0)) {
    stop("All observations must be positive.", call. = FALSE)
  }

  if (length(unique(x)) < 2L) {
    stop("x must contain at least two distinct values.", call. = FALSE)
  }

  required <- c("alpha", "beta", "lambda")
  lower <- .check_fit_bounds(lower, "lower")
  upper <- .check_fit_bounds(upper, "upper")

  if (any(upper <= lower)) {
    stop("Each upper bound must be greater than the lower bound.",
         call. = FALSE)
  }

  starts <- .make_nmeweibull_starts(x)

  if (!is.null(start)) {
    if (!is.numeric(start) || !all(required %in% names(start))) {
      stop("start must be a named numeric vector with names alpha, beta, ",
           "and lambda.", call. = FALSE)
    }

    start <- start[required]
    .check_nmeweibull_params(start[["alpha"]], start[["beta"]],
                             start[["lambda"]])
    starts <- rbind(start, starts)
  }

  inside <- apply(starts, 1L, function(s) all(s > lower & s < upper))
  starts <- starts[inside, , drop = FALSE]

  if (nrow(starts) == 0L) {
    stop("No starting value lies inside lower and upper.", call. = FALSE)
  }

  n <- length(x)
  log_m <- mean(log(x))
  y <- x / exp(log_m)
  log_lower <- log(lower)
  log_upper <- log(upper)

  # Internal parameters: (log alpha, log beta, log lambda*), where
  # lambda = lambda* * m^(-beta) and y = x / m. Then
  # loglik(x; theta) = loglik(y; alpha, beta, lambda*) - n * log(m).
  to_log_theta <- function(par) {
    c(par[1L], par[2L], par[3L] - exp(par[2L]) * log_m)
  }

  neg_loglik <- function(par) {
    log_theta <- to_log_theta(par)

    if (any(!is.finite(log_theta)) || any(log_theta <= log_lower) ||
        any(log_theta >= log_upper)) {
      return(1e100)
    }

    value <- -sum(.log_dnmeweibull(y, exp(par[1L]), exp(par[2L]),
                                   exp(par[3L])))

    if (!is.finite(value)) 1e100 else value
  }

  control <- utils::modifyList(list(maxit = 5000, reltol = 1e-10), control)

  best_fit <- NULL
  n_ok <- 0L

  for (i in seq_len(nrow(starts))) {
    par0 <- log(starts[i, ])
    par0[3L] <- par0[3L] + starts[i, "beta"] * log_m

    fit_nm <- tryCatch(
      stats::optim(par0, neg_loglik, method = "Nelder-Mead",
                   control = control),
      error = function(e) NULL
    )

    if (is.null(fit_nm) || fit_nm$value >= 1e100) {
      next
    }

    fit_bfgs <- tryCatch(
      stats::optim(fit_nm$par, neg_loglik, method = "BFGS",
                   control = control),
      error = function(e) NULL
    )

    candidate <- if (!is.null(fit_bfgs) && fit_bfgs$value <= fit_nm$value) {
      fit_bfgs
    } else {
      fit_nm
    }

    n_ok <- n_ok + 1L

    if (is.null(best_fit) || candidate$value < best_fit$value) {
      best_fit <- candidate
    }
  }

  if (is.null(best_fit)) {
    stop("Optimization failed for all starting values.", call. = FALSE)
  }

  log_theta <- to_log_theta(best_fit$par)
  estimate <- stats::setNames(exp(log_theta), required)

  loglik <- -best_fit$value - n * log_m
  k <- 3L

  near_upper <- is.finite(upper) & estimate >= 0.99 * upper
  near_lower <- lower > 0 & estimate <= 1.01 * lower

  hess <- NULL
  if (hessian) {
    neg_loglik_original <- function(log_par) {
      -sum(.log_dnmeweibull(x, exp(log_par[1L]), exp(log_par[2L]),
                            exp(log_par[3L])))
    }
    hess <- stats::optimHess(log_theta, neg_loglik_original)
    dimnames(hess) <- list(paste0("log_", required), paste0("log_", required))
  }

  structure(
    list(
      estimate = estimate,
      logLik = loglik,
      AIC = 2 * k - 2 * loglik,
      BIC = log(n) * k - 2 * loglik,
      convergence = best_fit$convergence,
      message = best_fit$message,
      n = n,
      boundary = any(near_upper | near_lower),
      n_starts = n_ok,
      hessian = hess,
      optim = best_fit
    ),
    class = "fitnmeweibull"
  )
}


#' Print a Fitted NME-Weibull Model
#'
#' @param x An object of class `"fitnmeweibull"` returned by
#'   [fitnmeweibull()].
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#' @export
#'
#' @examples
#' set.seed(1)
#' print(fitnmeweibull(rnmeweibull(100, 1, 2, 0.5)))
print.fitnmeweibull <- function(x, digits = 4L, ...) {
  cat("NME-Weibull maximum likelihood fit (n = ", x$n, ")\n\n", sep = "")
  print(signif(x$estimate, digits))
  cat("\nlogLik: ", format(x$logLik, digits = digits + 3L),
      "   AIC: ", format(x$AIC, digits = digits + 3L),
      "   BIC: ", format(x$BIC, digits = digits + 3L), "\n", sep = "")
  cat("Convergence code: ", x$convergence, "\n", sep = "")

  if (isTRUE(x$boundary)) {
    cat("Note: an estimate lies at a bound of the parameter space; a",
        "large alpha means the fit is practically Weibull.\n")
  }

  invisible(x)
}
