# is_unimodal recognises simple unimodal patterns
expect_true (is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1)))
expect_true (is_unimodal(c(0.5, 0.3, 0.1)))
expect_true (is_unimodal(c(0.1, 0.3, 0.5)))
expect_true (is_unimodal(c(0.1, 0.2, 0.2, 0.2, 0.1)))
expect_false(is_unimodal(c(0.1, 0.4, 0.2, 0.4, 0.1)))
expect_false(is_unimodal(c(0.5, 0.1, 0.4, 0.1)))

# is_logconcave recognises log-concave densities
pi <- dnorm(seq(-3, 3, length.out = 31))
pi <- pi / sum(pi)
expect_true(is_logconcave(pi))
x <- seq(-5, 5, length.out = 50)
pi <- 0.5 * dnorm(x, mean = -2) + 0.5 * dnorm(x, mean = 2)
pi <- pi / sum(pi)
expect_false(is_logconcave(pi))

# is_monotonic respects direction
expect_true (is_monotonic(c(0.5, 0.3, 0.1)))
expect_true (is_monotonic(c(0.5, 0.3, 0.1), direction = "decreasing"))
expect_false(is_monotonic(c(0.5, 0.3, 0.1), direction = "increasing"))
expect_true (is_monotonic(c(0.1, 0.3, 0.5), direction = "increasing"))
expect_false(is_monotonic(c(0.1, 0.4, 0.2)))
