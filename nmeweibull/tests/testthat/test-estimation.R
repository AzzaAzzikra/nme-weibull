test_that("fitnmeweibull returns the documented components", {
  set.seed(123)
  x <- rnmeweibull(300, alpha = 1, beta = 2, lambda = 0.5)
  fit <- fitnmeweibull(x, hessian = TRUE)

  expect_s3_class(fit, "fitnmeweibull")
  expect_named(fit$estimate, c("alpha", "beta", "lambda"))
  expect_true(all(fit$estimate > 0))
  expect_equal(fit$n, 300)
  expect_equal(fit$convergence, 0)
  expect_equal(fit$AIC, 6 - 2 * fit$logLik)
  expect_equal(fit$BIC, 3 * log(300) - 2 * fit$logLik)
  expect_equal(dim(fit$hessian), c(3L, 3L))
  expect_output(print(fit), "NME-Weibull maximum likelihood fit")
})

test_that("logLik equals the sum of log-densities at the estimate", {
  set.seed(42)
  x <- rnmeweibull(200, alpha = 0.5, beta = 1, lambda = 0.8)
  fit <- fitnmeweibull(x)
  est <- fit$estimate

  expect_equal(
    fit$logLik,
    sum(dnmeweibull(x, est["alpha"], est["beta"], est["lambda"], log = TRUE))
  )
})

test_that("estimates are close to the true values for a large sample", {
  set.seed(2026)
  x <- rnmeweibull(3000, alpha = 1, beta = 2, lambda = 0.5)
  est <- fitnmeweibull(x)$estimate

  expect_equal(unname(est["beta"]), 2, tolerance = 0.1)
  expect_equal(unname(est["lambda"]), 0.5, tolerance = 0.25)
})

test_that("rescaling the data only rescales lambda", {
  set.seed(7)
  x <- rnmeweibull(300, alpha = 0.8, beta = 3, lambda = 0.4)
  fit_1 <- fitnmeweibull(x)
  fit_100 <- fitnmeweibull(100 * x)

  expect_equal(fit_100$estimate[c("alpha", "beta")],
               fit_1$estimate[c("alpha", "beta")], tolerance = 1e-4)
  expect_equal(unname(fit_100$estimate["lambda"] * 100^fit_100$estimate["beta"]),
               unname(fit_1$estimate["lambda"]), tolerance = 1e-4)
  expect_equal(fit_100$logLik, fit_1$logLik - 300 * log(100),
               tolerance = 1e-6)
})

test_that("a user-supplied starting value is accepted", {
  set.seed(5)
  x <- rnmeweibull(100, alpha = 1, beta = 2, lambda = 0.5)
  fit <- fitnmeweibull(x, start = c(alpha = 1, beta = 2, lambda = 0.5))
  expect_equal(fit$n_starts, 76)
})

test_that("invalid data are rejected with clear errors", {
  expect_error(fitnmeweibull("a"), "x must be numeric")
  expect_error(fitnmeweibull(c(NA, Inf)), "at least one finite")
  expect_error(fitnmeweibull(c(1, 2, -3)), "must be positive")
  expect_error(fitnmeweibull(c(2, 2, 2)), "two distinct values")
  expect_error(fitnmeweibull(c(1, 2, 3), lower = c(alpha = 1)),
               "named numeric vector")
  expect_error(fitnmeweibull(c(1, 2, 3), start = c(1, 2, 3)),
               "named numeric vector")
})
