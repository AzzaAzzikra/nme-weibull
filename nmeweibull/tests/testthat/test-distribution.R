x <- c(0.1, 0.5, 1, 2, 3)

test_that("density integrates to one", {
  for (par in list(c(0.5, 1, 0.8), c(1, 2, 0.5), c(1.5, 3, 0.3),
                   c(0.1, 0.5, 2))) {
    total <- integrate(dnmeweibull, 0, Inf, alpha = par[1], beta = par[2],
                       lambda = par[3])$value
    expect_equal(total, 1, tolerance = 1e-4)
  }
})

test_that("density matches the closed form of the thesis", {
  alpha <- 1.3
  beta <- 2
  lambda <- 0.5
  z <- exp(-lambda * x^beta)
  expected <- alpha^2 * beta * lambda * x^(beta - 1) * z * (alpha + 1 + z) /
    (alpha + 1 - z)^3
  expect_equal(dnmeweibull(x, alpha, beta, lambda), expected)
  expect_equal(dnmeweibull(x, alpha, beta, lambda, log = TRUE), log(expected))
})

test_that("CDF, survival, and hazard are consistent", {
  alpha <- 0.7
  beta <- 1.5
  lambda <- 0.4
  cdf <- pnmeweibull(x, alpha, beta, lambda)
  surv <- snmeweibull(x, alpha, beta, lambda)

  expect_equal(surv, 1 - cdf)
  expect_equal(pnmeweibull(x, alpha, beta, lambda, lower.tail = FALSE), surv)
  expect_equal(pnmeweibull(x, alpha, beta, lambda, log.p = TRUE), log(cdf))
  expect_equal(snmeweibull(x, alpha, beta, lambda, log = TRUE), log(surv))
  expect_equal(hnmeweibull(x, alpha, beta, lambda),
               dnmeweibull(x, alpha, beta, lambda) / surv)
  expect_equal(hnmeweibull(x, alpha, beta, lambda, log = TRUE),
               log(dnmeweibull(x, alpha, beta, lambda) / surv))
})

test_that("functions handle the boundaries of the support", {
  expect_equal(pnmeweibull(c(-1, 0, Inf), 1, 2, 0.5), c(0, 0, 1))
  expect_equal(snmeweibull(c(-1, 0, Inf), 1, 2, 0.5), c(1, 1, 0))
  expect_equal(dnmeweibull(c(-1, Inf), 1, 2, 0.5), c(0, 0))
  expect_true(is.na(dnmeweibull(NA_real_, 1, 2, 0.5)))
})

test_that("hazard limits at zero follow the shape parameter", {
  expect_equal(hnmeweibull(0, 1, 0.5, 2), Inf)
  expect_equal(hnmeweibull(0, 1, 1, 2), 2 * (1 + 2) / 1)
  expect_equal(hnmeweibull(0, 1, 2, 2), 0)
  expect_equal(dnmeweibull(0, 1, 1, 2), 2 * (1 + 2) / 1)
})

test_that("quantile function inverts the CDF", {
  p <- c(0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999)
  for (par in list(c(0.1, 0.5, 2), c(1, 2, 0.5), c(5, 3, 1))) {
    q <- qnmeweibull(p, par[1], par[2], par[3])
    expect_equal(pnmeweibull(q, par[1], par[2], par[3]), p)
  }

  expect_equal(qnmeweibull(c(0, 1), 1, 2, 0.5), c(0, Inf))
  expect_equal(qnmeweibull(0.2, 1, 2, 0.5, lower.tail = FALSE),
               qnmeweibull(0.8, 1, 2, 0.5))
  expect_equal(qnmeweibull(log(0.3), 1, 2, 0.5, log.p = TRUE),
               qnmeweibull(0.3, 1, 2, 0.5))
})

test_that("the Weibull distribution is the limit for large alpha", {
  beta <- 2
  lambda <- 0.5
  expect_equal(pnmeweibull(x, 1e8, beta, lambda),
               pweibull(x, shape = beta, scale = lambda^(-1 / beta)),
               tolerance = 1e-6)
})

test_that("random generation is reproducible and has the right length", {
  set.seed(1)
  a <- rnmeweibull(50, 1, 2, 0.5)
  set.seed(1)
  b <- rnmeweibull(50, 1, 2, 0.5)

  expect_identical(a, b)
  expect_length(a, 50)
  expect_true(all(a > 0))
  expect_length(rnmeweibull(1:7, 1, 2, 0.5), 7)
  expect_length(rnmeweibull(0, 1, 2, 0.5), 0)
})

test_that("random generation follows the distribution", {
  set.seed(2026)
  sample <- rnmeweibull(5000, 0.5, 1, 0.8)
  test <- ks.test(sample, pnmeweibull, alpha = 0.5, beta = 1, lambda = 0.8)
  expect_gt(test$p.value, 0.001)
})

test_that("invalid input is rejected with clear errors", {
  expect_error(dnmeweibull(1, -1, 2, 0.5), "alpha must be positive")
  expect_error(dnmeweibull(1, 1, 0, 0.5), "beta must be positive")
  expect_error(dnmeweibull(1, 1, 2, Inf), "lambda must be finite")
  expect_error(dnmeweibull(1, 1, 2, NA_real_), "lambda must be finite")
  expect_error(dnmeweibull(1, c(1, 2), 2, 0.5), "single numeric value")
  expect_error(dnmeweibull("a", 1, 2, 0.5), "x must be numeric")
  expect_error(pnmeweibull("a", 1, 2, 0.5), "q must be numeric")
  expect_error(qnmeweibull(1.5, 1, 2, 0.5), "between 0 and 1")
  expect_error(qnmeweibull(-0.1, 1, 2, 0.5), "between 0 and 1")
  expect_error(rnmeweibull(-1, 1, 2, 0.5), "non-negative")
})
