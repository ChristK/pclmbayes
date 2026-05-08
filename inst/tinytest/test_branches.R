# =============================================================================
# Targeted tests for branches that the main test files don't exercise.
# These are not feature tests; they exist to keep coverage honest by hitting
# optional-arg paths and one graceful-degradation recovery.
# =============================================================================

set.seed(2027)
data(bloodlead, envir = environment(), package = "pclmbayes")
wb <- with(bloodlead, cbind(lower, upper))

# ---- pclm: auto a/b derivation when both NULL ---------------------------
fit_auto <- pclm(m = bloodlead$count, wide_breaks = wb,
                 ngrid = 50L, ndx = 11L, degree = 3L, penalty_order = 3L)
expect_inherits(fit_auto, "pclm")
expect_equal(fit_auto$grid[1L], min(wb[, 1L]))
expect_equal(fit_auto$grid[length(fit_auto$grid)], max(wb[, 2L]))

# ---- pclm: verbose path emits messages ----------------------------------
v_msg <- capture.output(
  pclm(m = bloodlead$count, wide_breaks = wb,
       a = 0, b = 80, ngrid = 50L, ndx = 11L, degree = 3L,
       penalty_order = 3L, tau = 1, verbose = TRUE),
  type = "message"
)
expect_true(any(grepl("\\[pclm\\] iter", v_msg)))

# ---- bpclm: phi_init, Sigma, delta_init non-NULL paths ------------------
warm <- pclm(m = bloodlead$count, wide_breaks = wb,
             a = 0, b = 80, ngrid = 50L,
             ndx = 11L, degree = 3L, penalty_order = 3L, tau = 1)
fit_init <- bpclm(m = bloodlead$count, wide_breaks = wb,
                  a = 0, b = 80, ngrid = 50L,
                  ndx = 11L, degree = 3L, penalty_order = 3L,
                  niter = 200L, burnin = 50L, adapt = 100L,
                  phi_init = warm$phi,
                  Sigma = warm$vcov,
                  delta_init = 0.5,
                  tau_init = 1,
                  seed = 1)
expect_inherits(fit_init, "bpclm")

# ---- bpclm: verbose path emits progress messages ------------------------
v_msg_b <- capture.output(
  bpclm(m = bloodlead$count, wide_breaks = wb,
        a = 0, b = 80, ngrid = 50L,
        ndx = 11L, degree = 3L, penalty_order = 3L,
        niter = 200L, burnin = 50L, adapt = 100L,
        verbose = TRUE, seed = 2),
  type = "message"
)
expect_true(any(grepl("\\[bpclm\\] iter", v_msg_b)))

# ---- bpclm: feasibility recovery Strategy 1 (oversmoothed warm-start) ---
# Bimodal mixture so the unconstrained warm-start fit is bimodal; requesting
# shape = "unimodal" then triggers the oversmoothing-recovery branch.
set.seed(2028)
y_bi <- c(rnorm(400, mean = 2, sd = 0.3),
          rnorm(400, mean = 7, sd = 0.4))
brk_bi <- seq(0, 10, by = 1)
m_bi   <- as.numeric(table(cut(y_bi, brk_bi, include.lowest = TRUE)))
fit_rec <- suppressWarnings(
  bpclm(m = m_bi, wide_breaks = brk_bi,
        a = 0, b = 10, ngrid = 50L,
        ndx = 11L, degree = 3L, penalty_order = 3L,
        niter = 300L, burnin = 100L, adapt = 100L,
        shape = "unimodal", seed = 3)
)
expect_inherits(fit_rec, "bpclm")
# Every retained draw must satisfy the constraint
expect_true(all(apply(fit_rec$pi_chain, 1L, is_unimodal)))

# ---- plot.bpclm with cred != x$cred_level (recomputes band) -------------
fit_b_small <- bpclm(m = bloodlead$count, wide_breaks = wb,
                     a = 0, b = 80, ngrid = 50L,
                     ndx = 11L, degree = 3L, penalty_order = 3L,
                     niter = 200L, burnin = 50L, adapt = 100L, seed = 4)
pdf(NULL); on.exit(dev.off(), add = TRUE)
out <- plot(fit_b_small, cred = 0.5)   # default is 0.9 -> recomputation path
expect_true(is.list(out))
expect_true(all(c("x", "y", "lower", "upper") %in% names(out)))
