test_that("is_unimodal recognises simple unimodal patterns", {
  expect_true (is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1)))
  expect_true (is_unimodal(c(0.5, 0.3, 0.1)))           # strictly decreasing
  expect_true (is_unimodal(c(0.1, 0.3, 0.5)))           # strictly increasing
  expect_true (is_unimodal(c(0.1, 0.2, 0.2, 0.2, 0.1))) # plateau allowed
  expect_false(is_unimodal(c(0.1, 0.4, 0.2, 0.4, 0.1))) # bimodal
  expect_false(is_unimodal(c(0.5, 0.1, 0.4, 0.1)))      # multiple modes
})

test_that("is_logconcave recognises log-concave densities", {
  # Discrete normal-like: log-concave
  pi <- dnorm(seq(-3, 3, length.out = 31))
  pi <- pi / sum(pi)
  expect_true(is_logconcave(pi))
  # Bimodal mixture: not log-concave
  x <- seq(-5, 5, length.out = 50)
  pi <- 0.5 * dnorm(x, mean = -2) + 0.5 * dnorm(x, mean = 2)
  pi <- pi / sum(pi)
  expect_false(is_logconcave(pi))
})

test_that("is_monotonic respects direction", {
  expect_true (is_monotonic(c(0.5, 0.3, 0.1)))
  expect_true (is_monotonic(c(0.5, 0.3, 0.1), direction = "decreasing"))
  expect_false(is_monotonic(c(0.5, 0.3, 0.1), direction = "increasing"))
  expect_true (is_monotonic(c(0.1, 0.3, 0.5), direction = "increasing"))
  expect_false(is_monotonic(c(0.1, 0.4, 0.2)))  # neither
})
