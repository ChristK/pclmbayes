# =============================================================================
# Tests for the frequentist fit pclm()
# -----------------------------------------------------------------------------
# These tests check that:
#   (a) the analytic gradient produced by .pclm_eval() agrees with finite
#       differences (this is essential because the same gradient is used by
#       the Bayesian sampler);
#   (b) on a simulated grouped sample drawn from a known density, the
#       frequentist fit recovers the underlying density well; and
#   (c) summary statistics (mean, sd, quantiles) match the simulating
#       distribution within sampling tolerance.
# =============================================================================

test_that("the gradient of the log-likelihood matches finite differences", {
  set.seed(2009)
  # Tiny problem so that finite differences are quick
  fine_breaks <- seq(0, 10, length.out = 21)
  mids        <- (head(fine_breaks, -1) + tail(fine_breaks, -1)) / 2
  bs          <- bspline_basis(mids, a = 0, b = 10, ndx = 5L, degree = 3L)
  K           <- bs$K
  wide_breaks <- seq(0, 10, by = 2)
  C           <- bin_matrix(wide_breaks, fine_breaks)
  # Generate plausible counts from a target density to avoid degenerate tests
  true_phi <- (1:K) - mean(1:K)        # increasing then re-centred
  pi_true  <- pclmbayes:::.softmax(bs$B %*% true_phi)
  ga_true  <- as.numeric(C %*% pi_true)
  m        <- as.integer(round(500 * ga_true))

  phi <- rnorm(K, sd = 0.3)
  phi <- phi - mean(phi)
  ev  <- pclmbayes:::.pclm_eval(phi, B = bs$B, C = C, m = m, compute_FI = FALSE)

  # Finite difference gradient. Use central differences but project the
  # perturbation onto the Σphi = 0 hyperplane so we compare like-for-like.
  eps <- 1e-5
  grad_fd <- numeric(K)
  for (k in seq_len(K)) {
    # Direction in the constrained space: e_k - mean(e_k)
    dvec <- numeric(K); dvec[k] <- 1; dvec <- dvec - mean(dvec)
    f_plus  <- pclmbayes:::.pclm_eval(phi + eps * dvec, B = bs$B, C = C,
                                      m = m, compute_FI = FALSE)$logL
    f_minus <- pclmbayes:::.pclm_eval(phi - eps * dvec, B = bs$B, C = C,
                                      m = m, compute_FI = FALSE)$logL
    grad_fd[k] <- (f_plus - f_minus) / (2 * eps)
  }
  # Compare against the analytic gradient projected onto the same hyperplane
  grad_proj <- ev$grad - mean(ev$grad)
  expect_lt(max(abs(grad_proj - grad_fd)), 1e-4)
})

test_that("pclm() recovers a gamma density from grouped counts", {
  set.seed(1)
  shape <- 5; rate <- 1
  n <- 5000
  y <- rgamma(n, shape = shape, rate = rate)
  brk <- seq(0, 20, by = 1)              # 1-unit-wide bins => mild grouping
  m   <- as.numeric(table(cut(y, brk, include.lowest = TRUE)))
  # Coarser bins for a tougher test
  brk_w <- c(0, 2, 4, 6, 8, 10, 12, 16, 20)
  m_w   <- as.numeric(table(cut(y, brk_w, include.lowest = TRUE)))

  fit <- pclm(m = m_w, wide_breaks = brk_w,
              a = 0, b = 20, ngrid = 80,
              ndx = 17L, degree = 3L, penalty_order = 3L)
  expect_s3_class(fit, "pclm")
  expect_true(fit$converged)

  # Recovered moments should be close to those of Gamma(5, 1): mean = 5, sd = sqrt(5)
  mom_hat <- pclmbayes:::.moments_from_pi(fit$pi, fit$grid_mid)
  expect_lt(abs(mom_hat["mean"] - shape / rate), 0.2)
  expect_lt(abs(mom_hat["sd"]   - sqrt(shape / rate^2)), 0.3)

  # Median should be close to qgamma(0.5, shape, rate)
  q50 <- quantile(fit, probs = 0.5)
  expect_lt(abs(q50 - qgamma(0.5, shape, rate)), 0.3)
})

test_that("pclm fit on the bloodlead data gives sensible values", {
  data(bloodlead, envir = environment(), package = "pclmbayes")
  fit <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 80,
              ndx = 17L, degree = 3L, penalty_order = 3L)
  expect_s3_class(fit, "pclm")
  expect_equal(sum(fit$pi), 1, tolerance = 1e-8)
  # Expected counts should sum to total observed counts
  expect_equal(sum(fit$fitted_counts), sum(bloodlead$count),
               tolerance = 1e-6)
  # Fitted mean ≈ paper's reported posterior mean of 21.8 (frequentist may
  # differ slightly but should be in the right ballpark).
  mom <- pclmbayes:::.moments_from_pi(fit$pi, fit$grid_mid)
  expect_gt(mom["mean"], 18)
  expect_lt(mom["mean"], 26)
})
