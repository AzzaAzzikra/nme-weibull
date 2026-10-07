alpha <- 1
beta <- 2
lambda <- 0.5
p <- c(0.90, 0.95, 0.99)

test_that("VaR equals the quantile function and increases with p", {
  var_p <- varnmeweibull(p, alpha, beta, lambda)
  expect_equal(var_p, qnmeweibull(p, alpha, beta, lambda))
  expect_true(all(diff(var_p) > 0))
})

test_that("TVaR agrees with the survival-function representation", {
  tvar_p <- tvarnmeweibull(p, alpha, beta, lambda)
  var_p <- varnmeweibull(p, alpha, beta, lambda)
  tvar_survival <- var_p + vapply(seq_along(p), function(i) {
    integrate(snmeweibull, var_p[i], Inf, alpha = alpha, beta = beta,
              lambda = lambda)$value / (1 - p[i])
  }, numeric(1))

  expect_equal(tvar_p, tvar_survival, tolerance = 1e-8)
  expect_true(all(tvar_p >= var_p))
})

test_that("TVaR at a small level is close to the mean", {
  mean_x <- integrate(snmeweibull, 0, Inf, alpha = alpha, beta = beta,
                      lambda = lambda)$value
  expect_equal(tvarnmeweibull(1e-10, alpha, beta, lambda), mean_x,
               tolerance = 1e-6)
})

test_that("RVaR lies between the VaR at its two levels", {
  rvar <- rvarnmeweibull(c(0.90, 0.95), c(0.95, 0.99), alpha, beta, lambda)
  expect_true(all(rvar >= varnmeweibull(c(0.90, 0.95), alpha, beta, lambda)))
  expect_true(all(rvar <= varnmeweibull(c(0.95, 0.99), alpha, beta, lambda)))
  expect_length(rvarnmeweibull(0.9, c(0.95, 0.99), alpha, beta, lambda), 2)
})

test_that("RVaR approaches TVaR as the upper level approaches one", {
  expect_equal(rvarnmeweibull(0.9, 1 - 1e-12, alpha, beta, lambda),
               tvarnmeweibull(0.9, alpha, beta, lambda), tolerance = 1e-4)
})

test_that("invalid probability levels are rejected", {
  expect_error(varnmeweibull(1, alpha, beta, lambda), "strictly between")
  expect_error(tvarnmeweibull(0, alpha, beta, lambda), "strictly between")
  expect_error(rvarnmeweibull(0.95, 0.90, alpha, beta, lambda),
               "smaller than p_upper")
  expect_error(rvarnmeweibull(c(0.1, 0.2), c(0.3, 0.4, 0.5), alpha, beta,
                              lambda), "same length")
  expect_error(varnmeweibull("a", alpha, beta, lambda), "p must be numeric")
})
