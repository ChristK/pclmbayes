# =============================================================================
# Tests for the Bayesian sampler bpclm()
# Integration tests with small chains for CRAN-friendly runtime.
# =============================================================================

# bpclm runs and returns valid posterior draws
set.seed(2025)
n   <- 800
y   <- rgamma(n, shape = 5, rate = 1)
brk_w <- c(0, 2, 4, 6, 8, 10, 12, 16)
m   <- as.numeric(table(cut(y, brk_w, include.lowest = TRUE)))

fit <- bpclm(m = m, wide_breaks = brk_w,
             a = 0, b = 16, ngrid = 60,
             ndx = 13L, degree = 3L, penalty_order = 3L,
             niter = 1000L, burnin = 200L, adapt = 200L,
             seed = 1)

expect_inherits(fit, "bpclm")
expect_equal(nrow(fit$phi_chain), 800L)
expect_equal(ncol(fit$phi_chain), fit$basis$K)
expect_equal(length(fit$tau_chain), 800L)

expect_equal(sum(fit$pi), 1, tolerance = 1e-8)
expect_true(all(abs(rowSums(fit$pi_chain) - 1) < 1e-8))

expect_true(fit$accept > 0.20)
expect_true(fit$accept < 0.95)

mom <- pclmbayes:::.moments_from_pi(fit$pi, fit$grid_mid)
expect_true(abs(mom["mean"] - 5) < 0.4)

# bpclm with shape = 'unimodal' produces only unimodal draws
data(bloodlead, envir = environment(), package = "pclmbayes")
fit <- bpclm(m = bloodlead$count,
             wide_breaks = with(bloodlead, cbind(lower, upper)),
             a = 0, b = 80, ngrid = 60,
             ndx = 13L, degree = 3L, penalty_order = 3L,
             niter = 600L, burnin = 100L, adapt = 200L,
             shape = "unimodal", seed = 7)
uni_flags <- apply(fit$pi_chain, 1L, is_unimodal)
expect_true(all(uni_flags))

# bpclm with shape = 'logconcave' produces only log-concave draws
fit_lc <- bpclm(m = bloodlead$count,
                wide_breaks = with(bloodlead, cbind(lower, upper)),
                a = 0, b = 80, ngrid = 60,
                ndx = 13L, degree = 3L, penalty_order = 3L,
                niter = 400L, burnin = 100L, adapt = 150L,
                shape = "logconcave", seed = 8)
lc_flags <- apply(fit_lc$pi_chain, 1L, is_logconcave)
expect_true(all(lc_flags))

# bpclm with shape = 'monotonic' (decreasing) on monotonic data
set.seed(99)
y_dec <- rexp(800, rate = 0.5)
brk_dec <- c(0, 1, 2, 3, 4, 6, 8, 12)
m_dec   <- as.numeric(table(cut(y_dec, brk_dec, include.lowest = TRUE)))
fit_mn  <- bpclm(m = m_dec, wide_breaks = brk_dec,
                 a = 0, b = 12, ngrid = 50,
                 ndx = 11L, degree = 3L, penalty_order = 3L,
                 niter = 400L, burnin = 100L, adapt = 150L,
                 shape = "monotonic",
                 shape_args = list(direction = "decreasing"),
                 seed = 21)
mn_flags <- apply(fit_mn$pi_chain, 1L,
                  function(p) is_monotonic(p, direction = "decreasing"))
expect_true(all(mn_flags))
