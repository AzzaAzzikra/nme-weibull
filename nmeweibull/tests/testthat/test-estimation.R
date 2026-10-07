test_that("MLE returns estimates close to generating parameters for simulated data", {
  set.seed(123)
  x <- rnmeweibull(n = 1000, alpha = 1, beta = 2, lambda = 0.5)

  fit <- fitnmeweibull(x)

  expect_equal(fit$convergence, 0)
  expect_true(is.finite(fit$logLik))
  expect_true(is.finite(fit$AIC))
  expect_true(is.finite(fit$BIC))

  expect_equal(unname(fit$estimate["alpha"]), 1, tolerance = 0.15)
  expect_equal(unname(fit$estimate["beta"]), 2, tolerance = 0.15)
  expect_equal(unname(fit$estimate["lambda"]), 0.5, tolerance = 0.15)
})
