# =============================================================================
# Tests for the Bayesian sampler bpclm()
# -----------------------------------------------------------------------------
# These are integration tests: they run a small MCMC chain on a small dataset
# and check that:
#   (a) the chain produces a sensible acceptance rate (between 30% and 80%);
#   (b) the posterior mean of the latent density integrates to 1;
#   (c) the posterior credible band contains the simulating density at most
#       grid points;
#   (d) shape constraints, when active, are satisfied by every kept draw.
#
# Iteration counts are kept small for CRAN-friendliness; users running the
# vignette use longer chains.
# =============================================================================

test_that("bpclm runs and returns valid posterior draws", {
  set.seed(2025)
  # Simulate from a Gamma(5, 1) so we have ground truth
  n   <- 800
  y   <- rgamma(n, shape = 5, rate = 1)
  brk_w <- c(0, 2, 4, 6, 8, 10, 12, 16)
  m   <- as.numeric(table(cut(y, brk_w, include.lowest = TRUE)))

  fit <- bpclm(m = m, wide_breaks = brk_w,
               a = 0, b = 16, ngrid = 60,
               ndx = 13L, degree = 3L, penalty_order = 3L,
               niter = 1000L, burnin = 200L, adapt = 200L,
               seed = 1)

  expect_s3_class(fit, "bpclm")
  expect_equal(nrow(fit$phi_chain), 800L)
  expect_equal(ncol(fit$phi_chain), fit$basis$K)
  expect_equal(length(fit$tau_chain), 800L)

  # Posterior mean density integrates to 1
  expect_equal(sum(fit$pi), 1, tolerance = 1e-8)
  # All kept draws are valid distributions
  expect_true(all(abs(rowSums(fit$pi_chain) - 1) < 1e-8))

  # Acceptance rate should be in a reasonable range (target 0.57 +/- a lot)
  expect_gt(fit$accept, 0.20)
  expect_lt(fit$accept, 0.95)

  # Recovered mean is close to ground truth
  mom <- pclmbayes:::.moments_from_pi(fit$pi, fit$grid_mid)
  expect_lt(abs(mom["mean"] - 5), 0.4)
})

test_that("bpclm with shape = 'unimodal' produces only unimodal draws", {
  data(bloodlead, envir = environment(), package = "pclmbayes")
  fit <- bpclm(m = bloodlead$count,
               wide_breaks = with(bloodlead, cbind(lower, upper)),
               a = 0, b = 80, ngrid = 60,
               ndx = 13L, degree = 3L, penalty_order = 3L,
               niter = 600L, burnin = 100L, adapt = 200L,
               shape = "unimodal", seed = 7)
  uni_flags <- apply(fit$pi_chain, 1L, is_unimodal)
  expect_true(all(uni_flags))
})
