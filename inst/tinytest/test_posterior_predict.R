# =============================================================================
# Tests for posterior_predict() and plot.pclm_posterior_predict()
# =============================================================================

set.seed(2026)
data(bloodlead, envir = environment(), package = "pclmbayes")

# Use a small bpclm fit (also serves as the parent of pclm via warmstart)
fit_b <- bpclm(m = bloodlead$count,
               wide_breaks = with(bloodlead, cbind(lower, upper)),
               a = 0, b = 80, ngrid = 80,
               ndx = 13L, degree = 3L, penalty_order = 3L,
               niter = 500L, burnin = 100L, adapt = 150L,
               seed = 5)

fit_p <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 80,
              ndx = 13L, degree = 3L, penalty_order = 3L)

# ---- posterior_predict.bpclm: type = "predictive" ----------------------
pp_b_pred <- posterior_predict(fit_b, type = "predictive", seed = 1)
expect_inherits(pp_b_pred, "pclm_posterior_predict")
expect_true(is.matrix(pp_b_pred$draws))
expect_equal(nrow(pp_b_pred$draws), nrow(fit_b$pi_chain))
expect_equal(ncol(pp_b_pred$draws), length(fit_b$grid_mid))
# Per-draw bin totals are preserved exactly (multinomial)
collapse <- function(x) pclmbayes:::.collapse_to_bins(x,
                          fit_b$wide_breaks, fit_b$grid)
for (s in c(1L, nrow(pp_b_pred$draws) %/% 2L, nrow(pp_b_pred$draws))) {
  expect_equal(collapse(pp_b_pred$draws[s, ]), fit_b$m,
               tolerance = 1e-8)
}
expect_equal(length(pp_b_pred$mean), length(fit_b$grid_mid))
expect_equal(length(pp_b_pred$lower), length(fit_b$grid_mid))
expect_equal(length(pp_b_pred$upper), length(fit_b$grid_mid))
expect_true(all(pp_b_pred$lower <= pp_b_pred$upper))

# ---- posterior_predict.bpclm: type = "rate" -----------------------------
pp_b_rate <- posterior_predict(fit_b, type = "rate")
expect_inherits(pp_b_rate, "pclm_posterior_predict")
expect_equal(dim(pp_b_rate$draws),
             c(nrow(fit_b$pi_chain), length(fit_b$grid_mid)))
# Rate draws should be non-negative
expect_true(all(pp_b_rate$draws >= 0))

# ---- posterior_predict.pclm: type = "rate" ------------------------------
pp_p_rate <- posterior_predict(fit_p, type = "rate")
expect_inherits(pp_p_rate, "pclm_posterior_predict")
expect_equal(nrow(pp_p_rate$draws), 1L)
expect_equal(ncol(pp_p_rate$draws), length(fit_p$grid_mid))

# ---- posterior_predict.pclm: type = "predictive" ------------------------
pp_p_pred <- posterior_predict(fit_p, type = "predictive",
                                n_draws = 50L, seed = 2)
expect_inherits(pp_p_pred, "pclm_posterior_predict")
expect_equal(dim(pp_p_pred$draws), c(50L, length(fit_p$grid_mid)))
collapse_p <- function(x) pclmbayes:::.collapse_to_bins(x,
                            fit_p$wide_breaks, fit_p$grid)
# Bin totals preserved on every draw
for (s in c(1L, 25L, 50L)) {
  expect_equal(collapse_p(pp_p_pred$draws[s, ]), fit_p$m,
               tolerance = 1e-8)
}

# ---- print and summary methods -----------------------------------------
out_print <- capture.output(print(pp_b_pred))
expect_true(length(out_print) > 0L)
out_print_rate <- capture.output(print(pp_b_rate))
expect_true(length(out_print_rate) > 0L)
out_summ <- capture.output(summary(pp_b_pred))
expect_true(length(out_summ) > 0L)
out_summ_rate <- capture.output(summary(pp_b_rate))
expect_true(length(out_summ_rate) > 0L)

# ---- plot.pclm_posterior_predict ---------------------------------------
# Use a null PDF device so nothing actually renders.
pdf(file = NULL)
res <- plot(pp_b_pred)
expect_true(is.list(res))
expect_equal(length(res$mean), length(fit_b$grid_mid))
# Also exercise show_bins = FALSE
res2 <- plot(pp_b_pred, show_bins = FALSE)
expect_true(is.list(res2))
# Rate flavour (different default title)
res3 <- plot(pp_b_rate)
expect_true(is.list(res3))
dev.off()
