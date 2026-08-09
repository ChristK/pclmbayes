# =============================================================================
# Basis-dimension default (.default_ndx) and weak-identification guard
# (.check_basis_dim).
#
# Background: the composite-link matrix C is J x I, so the wide-bin counts
# identify at most J - 1 directions in phi. The remaining K - (J - 1) are set
# by tau * P alone. When K is close to J those directions are smooth, the
# difference penalty barely resists them, and the fit develops large spurious
# excursions at the edges of the support. See ?pclm, section "Choosing the
# basis dimension".
# =============================================================================

# -----------------------------------------------------------------------------
# .default_ndx: the rule K = min(max(J + 7, 20), ngrid, 200)
# -----------------------------------------------------------------------------
K_of <- function(ngrid, degree, J) {
  pclmbayes:::.default_ndx(ngrid, degree, J) + degree
}

# Small J keeps K = 20, the Lambert & Eilers (2009) basis, so results for
# problems that already worked are unchanged.
expect_equal(K_of(100L, 3L,  7L), 20L)
expect_equal(K_of(100L, 3L, 12L), 20L)
expect_equal(K_of(100L, 3L, 13L), 20L)

# Beyond J = 13 the basis grows as J + 7, staying clear of the failure band.
expect_equal(K_of(100L, 3L, 14L), 21L)
expect_equal(K_of(100L, 3L, 19L), 26L)
expect_equal(K_of(100L, 3L, 30L), 37L)

# ngrid caps K: more B-splines than fine-grid intervals leaves the basis
# rank-deficient in its row space.
expect_equal(K_of(20L, 3L, 19L), 20L)
expect_true(K_of(15L, 3L, 40L) <= 15L)

# The 200 cap bounds the O(K^3) scoring cost.
expect_equal(K_of(1000L, 3L, 400L), 200L)

# ndx is always at least 1, and honours a non-default degree.
expect_true(pclmbayes:::.default_ndx(5L, 3L, 2L) >= 1L)
expect_equal(K_of(100L, 2L, 19L), 26L)

# The default always clears the guard whenever the grid allows it.
for (J in c(5L, 13L, 19L, 40L, 80L)) {
  expect_true(K_of(200L, 3L, J) >= J + 4L)
}

# -----------------------------------------------------------------------------
# .check_basis_dim: warns inside the band, silent outside
# -----------------------------------------------------------------------------
chk <- function(K, J, ndx_supplied = TRUE, ngrid = 100L) {
  pclmbayes:::.check_basis_dim(K = K, J = J, degree = 3L,
                               ndx_supplied = ndx_supplied, ngrid = ngrid)
}

# K = J + 4 is the documented safe threshold: silent at and above it.
expect_silent(chk(K = 23L, J = 19L))
expect_silent(chk(K = 40L, J = 19L))
# ... and warns below it, on both sides of K = J.
expect_warning(chk(K = 22L, J = 19L))
expect_warning(chk(K = 20L, J = 19L))   # the old default on 19 wide bins
expect_warning(chk(K = 19L, J = 19L))
expect_warning(chk(K = 15L, J = 19L))

# The message names the binding constraint.
w_ndx <- tryCatch(chk(K = 20L, J = 19L, ndx_supplied = TRUE),
                  warning = conditionMessage)
expect_true(grepl("Set `ndx` to at least 20", w_ndx, fixed = TRUE))

w_grid <- tryCatch(chk(K = 20L, J = 19L, ndx_supplied = FALSE, ngrid = 20L),
                   warning = conditionMessage)
expect_true(grepl("Increase `ngrid` to at least 23", w_grid, fixed = TRUE))

# When ngrid is ample, a defaulted ndx that still falls short (the 200 cap)
# points at ndx, not ngrid.
w_cap <- tryCatch(chk(K = 200L, J = 300L, ndx_supplied = FALSE, ngrid = 1000L),
                  warning = conditionMessage)
expect_true(grepl("Set `ndx`", w_cap, fixed = TRUE))

# The warning reports K, J and the count of identified directions.
expect_true(grepl("K = 20", w_ndx, fixed = TRUE))
expect_true(grepl("J = 19", w_ndx, fixed = TRUE))
expect_true(grepl("J - 1 = 18", w_ndx, fixed = TRUE))

# -----------------------------------------------------------------------------
# Wiring: pclm(), bpclm() and pclm_exact() default and guard correctly
# -----------------------------------------------------------------------------
data(bloodlead, envir = environment(), package = "pclmbayes")
bl_wb <- with(bloodlead, cbind(lower, upper))   # J = 7

# Defaulted ndx reproduces the historical K = 20 for these small-J examples,
# so previously published fits are unchanged.
fit_def <- pclm(m = bloodlead$count, wide_breaks = bl_wb,
                a = 0, b = 80, ngrid = 80)
fit_17  <- pclm(m = bloodlead$count, wide_breaks = bl_wb,
                a = 0, b = 80, ngrid = 80, ndx = 17L)
expect_equal(fit_def$basis$ndx, 17L)
expect_equal(fit_def$basis$K,   20L)
expect_equal(fit_def$pi, fit_17$pi, tolerance = 1e-12)

# A defaulted fit never trips its own guard when the grid is ample.
expect_silent(pclm(m = bloodlead$count, wide_breaks = bl_wb,
                   a = 0, b = 80, ngrid = 80))

# Many wide bins: the default grows the basis past the failure band.
set.seed(11)
brk_many <- seq(0, 80, by = 4)                    # J = 20 wide bins
m_many   <- as.numeric(table(cut(rgamma(20000, 6, 0.25), brk_many,
                                 include.lowest = TRUE)))
fit_many <- pclm(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                 ngrid = 100)
expect_equal(fit_many$basis$K, 27L)               # max(J + 7, 20) = 27
expect_true(fit_many$basis$K >= length(m_many) + 4L)

# An explicit ndx inside the band warns, and check_basis = FALSE silences it.
expect_warning(pclm(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                    ngrid = 100, ndx = 17L))
expect_silent(pclm(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                   ngrid = 100, ndx = 17L, check_basis = FALSE))

# pclm_exact() shares the rule and the guard.
fit_ex <- pclm_exact(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                     ngrid = 100)
expect_equal(fit_ex$basis$K, 27L)
expect_warning(pclm_exact(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                          ngrid = 100, ndx = 17L))
expect_silent(pclm_exact(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                         ngrid = 100, ndx = 17L, check_basis = FALSE))

# bpclm() inherits both through its warm start, and warns exactly once.
fit_b <- bpclm(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
               ngrid = 100, niter = 200L, burnin = 50L, adapt = 20L,
               seed = 4)
expect_equal(fit_b$basis$K, 27L)
expect_warning(bpclm(m = m_many, wide_breaks = brk_many, a = 0, b = 80,
                     ngrid = 100, ndx = 17L, niter = 200L, burnin = 50L,
                     adapt = 20L, seed = 4))

# -----------------------------------------------------------------------------
# Regression: the failure this guard exists to prevent
# -----------------------------------------------------------------------------
# Ungrouping the ONS mid-2024 UK population from 5-year bands (plus the open
# 90+ bin, J = 19) back to single years of age. With K = 20 the fit put
# ~2.37 million people at age 0 against a truth of 668 thousand; the default
# basis must keep the error small at the boundary.
data(ons_pop_uk_2024, envir = environment(), package = "pclmbayes")
sy_truth <- as.numeric(
  ons_pop_uk_2024$count[ons_pop_uk_2024$age_label != "90+"])
m_open <- c(as.numeric(tapply(sy_truth, rep(seq(0, 85, 5), each = 5), sum)),
            as.numeric(ons_pop_uk_2024$count[
              ons_pop_uk_2024$age_label == "90+"]))
wb_open <- rbind(cbind(seq(0, 85, 5), seq(5, 90, 5)), c(90, 100))

fit_ons <- pclm(m = m_open, wide_breaks = wb_open, a = 0, b = 100,
                ngrid = 100)
sy_fit <- as.numeric(bin_matrix(seq(0, 90, 1), fit_ons$grid) %*%
                       fit_ons$pi)[1:90] * sum(fit_ons$m)
err <- abs(sy_fit - sy_truth) / sy_truth

expect_true(fit_ons$basis$K >= length(m_open) + 4L)
expect_true(max(err[1:21]) < 0.05)      # was 2.54 with K = 20
expect_true(mean(err[1:21]) < 0.02)     # was 0.30 with K = 20
expect_true(max(abs(fit_ons$phi)) < 25) # was 157 with K = 20

# The unconstrained coefficients stay bounded because the unidentified
# directions are now high-frequency, where the difference penalty bites.
fit_bad <- suppressWarnings(
  pclm(m = m_open, wide_breaks = wb_open, a = 0, b = 100, ngrid = 100,
       ndx = 17L))
expect_true(max(abs(fit_bad$phi)) > 10 * max(abs(fit_ons$phi)))
