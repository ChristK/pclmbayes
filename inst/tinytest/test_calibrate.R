# =============================================================================
# Tests for calibrate.pclm and calibrate.bpclm
# =============================================================================

set.seed(2026)
data(bloodlead, envir = environment(), package = "pclmbayes")

# ---- calibrate.pclm ----------------------------------------------------
fit_p <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 80,
              ndx = 17L, degree = 3L, penalty_order = 3L)

# Pre-calibration: bin totals not exact
expect_true(max(abs(fit_p$fitted_counts - fit_p$m)) > 1e-3)

fit_pc <- calibrate(fit_p)
expect_inherits(fit_pc, "pclm")
expect_true(isTRUE(fit_pc$calibrated))
# Bin totals reproduced exactly to numerical precision
expect_true(max(abs(fit_pc$fitted_counts - fit_pc$m)) < 1e-6)
expect_equal(sum(fit_pc$pi), 1, tolerance = 1e-8)

# ---- calibrate.bpclm ---------------------------------------------------
fit_b <- bpclm(m = bloodlead$count,
               wide_breaks = with(bloodlead, cbind(lower, upper)),
               a = 0, b = 80, ngrid = 80,
               ndx = 13L, degree = 3L, penalty_order = 3L,
               niter = 500L, burnin = 100L, adapt = 150L,
               seed = 13)
fit_bc <- calibrate(fit_b)
expect_inherits(fit_bc, "bpclm")
expect_true(isTRUE(fit_bc$calibrated))
# Each posterior draw exactly preserves bin totals
m_plus <- sum(fit_bc$m)
target <- fit_bc$m / m_plus
g_chain <- fit_bc$pi_chain %*% t(fit_bc$C)
expect_true(max(abs(sweep(g_chain, 2L, target))) < 1e-6)
# pi_lower/upper still aligned in length
expect_equal(length(fit_bc$pi_lower), length(fit_bc$pi))
expect_equal(length(fit_bc$pi_upper), length(fit_bc$pi))

# ---- pclm_exact: exercises SQP iteration -------------------------------
set.seed(42)
fit_e <- pclm_exact(m = bloodlead$count,
                    wide_breaks = with(bloodlead, cbind(lower, upper)),
                    a = 0, b = 80, ngrid = 80,
                    ndx = 13L, degree = 3L, penalty_order = 3L,
                    max_iter = 200L, tol = 1e-9)
expect_inherits(fit_e, "pclm_exact")
expect_inherits(fit_e, "pclm")
# Polish step guarantees exact preservation when bins partition
# the support; allow a generous tolerance to absorb numerical noise.
expect_true(max(abs(fit_e$fitted_counts - fit_e$m)) < 1e-3)
