#' Log Density Function of the NME-Weibull Distribution
#'
#' Internal function to compute the log-density of the NME-Weibull distribution.
#'
#' @param x Numeric vector of positive observations.
#' @param alpha Additional parameter of the NME-Weibull distribution.
#' @param beta Shape parameter of the baseline Weibull distribution.
#' @param lambda Scale/rate parameter of the baseline Weibull distribution.
#'
#' @return A numeric vector of log-density values.
#' @noRd
.log_dnmeweibull <- function(x, alpha, beta, lambda) {
  .check_nmeweibull_params(alpha, beta, lambda)

  log_density <- rep(NA_real_, length(x))

  valid <- x > 0 & is.finite(x) & !is.na(x)

  log_x <- log(x[valid])
  x_beta <- exp(beta * log_x)
  eta <- lambda * x_beta

  exp_term <- exp(-eta)
  weibull_cdf <- 1 - exp_term

  log_density[valid] <-
    2 * log(alpha) +
    log(beta) +
    log(lambda) +
    (beta - 1) * log_x -
    eta +
    log(alpha + 2 - weibull_cdf) -
    3 * log(alpha + weibull_cdf)

  log_density[x < 0 & !is.na(x)] <- -Inf

  log_density
}


#' Generate Starting Values for NME-Weibull MLE
#'
#' Internal function to generate multiple starting values for optimization.
#'
#' @param x Numeric vector of positive observations.
#'
#' @return A matrix of starting values.
#' @noRd
.make_nmeweibull_starts <- function(x) {
  alpha_grid <- c(0.5, 1, 2)
  beta_grid <- c(0.5, 1, 1.5, 2, 3)
  lambda_multiplier <- c(0.25, 0.5, 1, 2, 4)

  starts <- list()
  index <- 1

  for (alpha in alpha_grid) {
    for (beta in beta_grid) {
      lambda_base <- 1 / (mean(x)^beta)

      for (mult in lambda_multiplier) {
        starts[[index]] <- c(
          alpha = alpha,
          beta = beta,
          lambda = lambda_base * mult
        )
        index <- index + 1
      }
    }
  }

  do.call(rbind, starts)
}


#' Fit the NME-Weibull Distribution by Maximum Likelihood Estimation
#'
#' Estimates the parameters of the New Modified Exponential-Weibull
#' distribution using maximum likelihood estimation.
#'
#' @param x Numeric vector of positive observations.
#' @param start Optional named numeric vector of starting values with names
#'   alpha, beta, and lambda.
#' @param lower Named numeric vector of lower bounds for alpha, beta, and lambda.
#' @param upper Named numeric vector of upper bounds for alpha, beta, and lambda.
#' @param hessian Logical; if TRUE, the Hessian matrix is computed at the final
#'   estimate.
#' @param control A list of control parameters passed to [stats::optim()].
#'
#' @return A list containing parameter estimates, log-likelihood, AIC, BIC,
#'   convergence code, optimization message, and the full optim object.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnmeweibull(n = 100, alpha = 1, beta = 2, lambda = 0.5)
#' fitnmeweibull(x)
fitnmeweibull <- function(x,
                          start = NULL,
                          lower = c(alpha = 1e-8, beta = 1e-8, lambda = 1e-8),
                          upper = c(alpha = 1e4, beta = 1e4, lambda = 1e4),
                          hessian = FALSE,
                          control = list()) {
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }

  x <- x[is.finite(x)]

  if (length(x) == 0) {
    stop("x must contain at least one finite observation.", call. = FALSE)
  }

  if (any(x <= 0)) {
    stop("All observations must be positive.", call. = FALSE)
  }

  required_names <- c("alpha", "beta", "lambda")

  if (!all(required_names %in% names(lower)) ||
      !all(required_names %in% names(upper))) {
    stop("lower and upper must be named vectors with names alpha, beta, and lambda.",
         call. = FALSE)
  }

  lower <- lower[required_names]
  upper <- upper[required_names]

  if (any(lower <= 0) || any(upper <= lower)) {
    stop("Bounds must be positive and upper must be greater than lower.",
         call. = FALSE)
  }

  auto_starts <- .make_nmeweibull_starts(x)

  if (!is.null(start)) {
    if (!all(required_names %in% names(start))) {
      stop("start must be a named vector with names alpha, beta, and lambda.",
           call. = FALSE)
    }

    start <- matrix(start[required_names], nrow = 1)
    colnames(start) <- required_names
    starts <- rbind(start, auto_starts)
  } else {
    starts <- auto_starts
  }

  starts <- starts[
    starts[, "alpha"] > lower["alpha"] &
      starts[, "beta"] > lower["beta"] &
      starts[, "lambda"] > lower["lambda"] &
      starts[, "alpha"] < upper["alpha"] &
      starts[, "beta"] < upper["beta"] &
      starts[, "lambda"] < upper["lambda"],
    ,
    drop = FALSE
  ]

  neg_loglik <- function(log_par) {
    if (any(!is.finite(log_par))) {
      return(1e100)
    }

    if (any(log_par < log(lower)) || any(log_par > log(upper))) {
      return(1e100)
    }

    par <- exp(log_par)

    alpha <- par[1]
    beta <- par[2]
    lambda <- par[3]

    log_density <- .log_dnmeweibull(
      x = x,
      alpha = alpha,
      beta = beta,
      lambda = lambda
    )

    if (any(!is.finite(log_density))) {
      return(1e100)
    }

    -sum(log_density)
  }

  default_control <- list(maxit = 5000, reltol = 1e-10)
  control <- utils::modifyList(default_control, control)

  best_fit <- NULL

  for (i in seq_len(nrow(starts))) {
    initial <- starts[i, required_names]

    fit_nm <- try(
      stats::optim(
        par = log(initial),
        fn = neg_loglik,
        method = "Nelder-Mead",
        control = control
      ),
      silent = TRUE
    )

    if (inherits(fit_nm, "try-error") || !is.finite(fit_nm$value)) {
      next
    }

    fit_bfgs <- try(
      stats::optim(
        par = fit_nm$par,
        fn = neg_loglik,
        method = "BFGS",
        control = control
      ),
      silent = TRUE
    )

    if (!inherits(fit_bfgs, "try-error") && is.finite(fit_bfgs$value)) {
      candidate <- fit_bfgs
    } else {
      candidate <- fit_nm
    }

    if (is.null(best_fit) || candidate$value < best_fit$value) {
      best_fit <- candidate
    }
  }

  if (is.null(best_fit)) {
    stop("Optimization failed for all starting values.", call. = FALSE)
  }

  estimate <- exp(best_fit$par)
  names(estimate) <- required_names

  loglik <- -best_fit$value
  k <- length(estimate)
  n <- length(x)

  aic <- 2 * k - 2 * loglik
  bic <- log(n) * k - 2 * loglik

  hess <- NULL

  if (hessian) {
    hess <- stats::optimHess(
      par = best_fit$par,
      fn = neg_loglik
    )
  }

  result <- list(
    estimate = estimate,
    logLik = loglik,
    AIC = aic,
    BIC = bic,
    convergence = best_fit$convergence,
    message = best_fit$message,
    n = n,
    hessian = hess,
    optim = best_fit
  )

  class(result) <- "fitnmeweibull"

  result
}
