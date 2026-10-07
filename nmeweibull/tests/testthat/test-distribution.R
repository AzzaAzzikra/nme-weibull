test_that("CDF and survival functions are consistent", {
  x <- c(0.5, 1, 2)
  p <- pnmeweibull(x, alpha = 1, beta = 2, lambda = 0.5)
  s <- snmeweibull(x, alpha = 1, beta = 2, lambda = 0.5)

  expect_equal(p + s, rep(1, length(x)), tolerance = 1e-8)
})

test_that("hazard equals density divided by survival", {
  x <- c(0.5, 1, 2)

  f <- dnmeweibull(x, alpha = 1, beta = 2, lambda = 0.5)
  s <- snmeweibull(x, alpha = 1, beta = 2, lambda = 0.5)
  h <- hnmeweibull(x, alpha = 1, beta = 2, lambda = 0.5)

  expect_equal(h, f / s, tolerance = 1e-8)
})

test_that("quantile function is inverse of CDF", {
  p <- c(0.25, 0.5, 0.75)
  q <- qnmeweibull(p, alpha = 1, beta = 2, lambda = 0.5)

  expect_equal(
    pnmeweibull(q, alpha = 1, beta = 2, lambda = 0.5),
    p,
    tolerance = 1e-8
  )
})

test_that("random generation produces non-negative values", {
  set.seed(123)
  x <- rnmeweibull(100, alpha = 1, beta = 2, lambda = 0.5)

  expect_length(x, 100)
  expect_true(all(x >= 0))
})
