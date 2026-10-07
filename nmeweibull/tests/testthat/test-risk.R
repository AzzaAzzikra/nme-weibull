test_that("VaR is consistent with the quantile function", {
  var95 <- varnmeweibull(0.95, alpha = 1, beta = 2, lambda = 0.5)

  expect_equal(
    pnmeweibull(var95, alpha = 1, beta = 2, lambda = 0.5),
    0.95,
    tolerance = 1e-8
  )
})

test_that("TVaR is greater than VaR", {
  p <- c(0.90, 0.95, 0.99)

  var <- varnmeweibull(p, alpha = 1, beta = 2, lambda = 0.5)
  tvar <- tvarnmeweibull(p, alpha = 1, beta = 2, lambda = 0.5)

  expect_true(all(tvar > var))
})

test_that("RVaR lies between lower and upper VaR", {
  p_lower <- 0.90
  p_upper <- 0.95

  var_lower <- varnmeweibull(p_lower, alpha = 1, beta = 2, lambda = 0.5)
  var_upper <- varnmeweibull(p_upper, alpha = 1, beta = 2, lambda = 0.5)
  rvar <- rvarnmeweibull(p_lower, p_upper, alpha = 1, beta = 2, lambda = 0.5)

  expect_true(var_lower <= rvar)
  expect_true(rvar <= var_upper)
})
