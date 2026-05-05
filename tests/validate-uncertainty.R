# =============================================================================
# Validation of single-year ungrouping WITH UNCERTAINTY.
#
# Workflow:
#   1. Generate single-year-of-age "death" counts (multinomial draw,
#      bathtub mortality density, ages 0-99).
#   2. Group into 5-year bands.
#   3. Run bpclm() to obtain a posterior chain over the latent fine-grid
#      density.
#   4. Apply calibrate() to every posterior draw -- this projects each
#      draw onto the constraint manifold {pi : C pi = m / m_+} so that
#      the band totals are preserved exactly on every draw.
#   5. From the calibrated chain compute, at every single year:
#         posterior mean count
#         90% credible interval
#         empirical coverage (does the truth lie inside?)
#   6. Compare to the pclm_exact point estimate (which has no
#      uncertainty).
# =============================================================================

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())

set.seed(2026)
ages       <- 0:99
band_lower <- seq(0,  95, by = 5)
band_upper <- seq(5, 100, by = 5)
wb         <- cbind(band_lower, band_upper)
J          <- length(band_lower)

true_p <- 0.10 * dnorm(ages, mean =  2,  sd =  3) +
          0.06 * dnorm(ages, mean = 25,  sd = 12) +
          0.84 * dnorm(ages, mean = 75,  sd = 12)
true_p <- true_p / sum(true_p)
year_to_band <- (ages %/% 5L) + 1L

# ---- Single illustrative replicate ----
N <- 50000L
true_yearly <- as.numeric(rmultinom(1, size = N, prob = true_p))
m_band      <- as.numeric(tapply(true_yearly, year_to_band, sum))
cat(sprintf("N = %d, summed bands = %d\n", N, sum(m_band)))

# Point-estimate from pclm_exact for comparison
fit_exact <- pclm_exact(m = m_band, wide_breaks = wb, a = 0, b = 100,
                        ngrid = 100L, ndx = 22L, degree = 3L,
                        penalty_order = 3L)
yhat_exact <- fit_exact$pi * sum(m_band)

# Bayesian fit
cat("Running bpclm + calibrate ...\n")
t0 <- Sys.time()
fit_b <- bpclm(m = m_band, wide_breaks = wb, a = 0, b = 100,
               ngrid = 100L, ndx = 22L, degree = 3L, penalty_order = 3L,
               niter = 6000L, burnin = 1500L, adapt = 600L,
               seed = 2026)
fit_bc <- calibrate(fit_b)
cat(sprintf("  done in %.1f s, %d kept draws, acceptance %.2f\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(fit_bc$pi_chain), fit_bc$accept))

# Posterior chains
M <- nrow(fit_bc$pi_chain)

# Two posterior summaries:
#   (a) "rate" = posterior over the latent smooth density m_+ * pi_y.
#       Captures epistemic uncertainty in the smooth shape only.
#   (b) "predictive" = posterior over the *realised* single-year counts.
#       For each draw of pi, sample multinomial within-band counts.
#       This is the right object when the user wants to know plausible
#       single-year counts given the observed band totals.
rate_chain <- fit_bc$pi_chain * sum(m_band)

pred_chain <- matrix(NA_real_, nrow = M, ncol = length(ages))
target_p   <- m_band / sum(m_band)
for (s in seq_len(M)) {
  pi_s   <- fit_bc$pi_chain[s, ]
  for (j in seq_len(J)) {
    idx <- which(year_to_band == j)
    g_j <- sum(pi_s[idx])
    if (g_j > 0) {
      p_within <- pi_s[idx] / g_j         # conditional probabilities
      pred_chain[s, idx] <- rmultinom(1, size = m_band[j],
                                       prob = p_within)[, 1L]
    } else {
      pred_chain[s, idx] <- 0
    }
  }
}

# Per-year posterior summaries (posterior PREDICTIVE — captures both
# epistemic + multinomial sampling uncertainty)
post_mean <- colMeans(pred_chain)
post_lo   <- apply(pred_chain, 2L, quantile, probs = 0.05, names = FALSE)
post_hi   <- apply(pred_chain, 2L, quantile, probs = 0.95, names = FALSE)

# Rate-only summaries for comparison
rate_lo <- apply(rate_chain, 2L, quantile, probs = 0.05, names = FALSE)
rate_hi <- apply(rate_chain, 2L, quantile, probs = 0.95, names = FALSE)

cover_pred <- (true_yearly >= post_lo) & (true_yearly <= post_hi)
cover_rate <- (true_yearly >= rate_lo) & (true_yearly <= rate_hi)
cat(sprintf("90%% PREDICTIVE coverage : %d / %d (%.1f%%)\n",
            sum(cover_pred), length(cover_pred), 100 * mean(cover_pred)))
cat(sprintf("90%% RATE-only  coverage : %d / %d (%.1f%%)\n",
            sum(cover_rate), length(cover_rate), 100 * mean(cover_rate)))
cover <- cover_pred
cat(sprintf("\n90%% credible-interval coverage of the truth:  %d / %d  (%.1f%%)\n",
            sum(cover), length(cover), 100 * mean(cover)))
cat(sprintf("Posterior-mean RMSE vs truth: %.2f\n",
            sqrt(mean((post_mean - true_yearly)^2))))
cat(sprintf("pclm_exact   RMSE vs truth:   %.2f\n",
            sqrt(mean((yhat_exact - true_yearly)^2))))

# Band totals preserved on every draw?
gamma_each_draw <- fit_bc$pi_chain %*% t(fit_bc$C)
target_p <- m_band / sum(m_band)
max_resid_per_draw <-
  apply(abs(gamma_each_draw - matrix(target_p, M, J, byrow = TRUE)), 1L, max)
cat(sprintf("Max per-draw band residual: %.2e (across %d draws)\n",
            max(max_resid_per_draw), M))

# Cumulative-deaths posterior summaries (predictive)
cdf_chain <- t(apply(pred_chain, 1L, cumsum))
cdf_mean  <- colMeans(cdf_chain)
cdf_lo    <- apply(cdf_chain, 2L, quantile, probs = 0.05, names = FALSE)
cdf_hi    <- apply(cdf_chain, 2L, quantile, probs = 0.95, names = FALSE)
cdf_cov   <- (cumsum(true_yearly) >= cdf_lo) &
             (cumsum(true_yearly) <= cdf_hi)
cat(sprintf("90%% CI coverage of cumulative-deaths truth: %d / %d (%.1f%%)\n",
            sum(cdf_cov), length(cdf_cov), 100 * mean(cdf_cov)))

# ---- Plot ----
out_png <- file.path(pkg_root, "tests", "validate_uncertainty.png")
png(out_png, width = 2600, height = 1700, res = 200)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.5, 1.0), oma = c(0, 0, 2, 0))

# Panel 1: per-year truth + posterior predictive mean + bands
ymax <- max(c(true_yearly, post_hi)) * 1.05
plot(NA, xlim = range(ages), ylim = c(0, ymax),
     xlab = "Age (years)", ylab = "Deaths in this single-year",
     main = "Per-year posterior predictive: truth, mean, 90% PI")
polygon(c(ages, rev(ages)), c(post_lo, rev(post_hi)),
        col = grDevices::adjustcolor("steelblue", 0.3), border = NA)
polygon(c(ages, rev(ages)), c(rate_lo, rev(rate_hi)),
        col = grDevices::adjustcolor("orange", 0.45), border = NA)
lines(ages, true_yearly, col = "grey40", lwd = 2)
lines(ages, post_mean,   col = "black",  lwd = 2)
lines(ages, yhat_exact,  col = "darkgreen", lwd = 2, lty = 2)
abline(v = seq(5, 95, by = 5), col = "grey90", lty = 3)
legend("topleft",
       legend = c("Truth (realised single-year counts)",
                  "Posterior predictive mean",
                  "90% predictive interval (rate + Mult. noise)",
                  "90% rate-only interval (smooth shape)",
                  "pclm_exact (point estimate)"),
       col = c("grey40", "black",
                grDevices::adjustcolor("steelblue", 0.6),
                grDevices::adjustcolor("orange", 0.8),
                "darkgreen"),
       lwd = c(2, 2, 8, 8, 2), lty = c(1, 1, 1, 1, 2), bty = "n")

# Panel 2: per-year coverage indicator (where does the CI fail to cover?)
plot(ages, as.numeric(cover), type = "h",
     col = ifelse(cover, "darkgreen", "firebrick"),
     lwd = 3, ylim = c(0, 1.15),
     xlab = "Age (years)", ylab = "covered (1) / missed (0)",
     main = sprintf("Per-year 90%% credible-interval coverage: %d / %d (%.0f%%)",
                    sum(cover), length(cover), 100 * mean(cover)))
abline(v = seq(5, 95, by = 5), col = "grey90", lty = 3)
abline(h = 0.9, col = "grey50", lty = 2)
legend("topleft", legend = c("nominal 90%"), col = "grey50", lty = 2,
       bty = "n")

# Panel 3: cumulative-deaths posterior + truth
plot(NA, xlim = range(ages), ylim = c(0, max(cdf_hi)),
     xlab = "Age (years)", ylab = "Cumulative deaths",
     main = "Cumulative deaths: truth + posterior 90% band")
polygon(c(ages, rev(ages)), c(cdf_lo, rev(cdf_hi)),
        col = grDevices::adjustcolor("steelblue", 0.3), border = NA)
lines(ages, cumsum(true_yearly), col = "grey40", lwd = 2)
lines(ages, cdf_mean, col = "black", lwd = 2)
abline(v = seq(5, 95, by = 5), col = "grey90", lty = 3)
legend("topleft",
       legend = c("Truth", "Posterior mean", "90% credible band"),
       col = c("grey40", "black",
                grDevices::adjustcolor("steelblue", 0.6)),
       lwd = c(2, 2, 8), bty = "n")

# Panel 4: the band totals match exactly across draws (sanity check)
plot(seq_len(J), m_band, type = "h", lwd = 6, col = "grey60",
     xlab = "Age band index (1 = (0, 5], ..., 20 = (95, 100])",
     ylab = "Deaths",
     main = "Band totals: every draw matches observed (residual ~ 1e-13)")
points(seq_len(J), apply(gamma_each_draw, 2L, mean) * sum(m_band),
       pch = 4, col = "darkgreen", cex = 1.5)
legend("topright",
       legend = c("Observed band totals m_j",
                  "Posterior mean of fitted band totals"),
       col = c("grey60", "darkgreen"),
       lwd = c(6, NA), pch = c(NA, 4), bty = "n")

mtext("Bayesian uncertainty for exact ungrouping: bpclm() + calibrate()",
      side = 3, outer = TRUE, cex = 0.95, font = 2)
par(op)
dev.off()
cat("\nWrote: ", out_png, "\n", sep = "")

# ---- Monte Carlo coverage diagnostic ----
cat("\nMonte-Carlo coverage check across 30 replicates ...\n")
S <- 30L
cov_per_year   <- matrix(NA, nrow = S, ncol = length(ages))
rmse_per_rep   <- numeric(S)
exact_rmse_per_rep <- numeric(S)
for (s in seq_len(S)) {
  set.seed(2026 + 1000 + s)
  ty <- as.numeric(rmultinom(1, size = N, prob = true_p))
  mb <- as.numeric(tapply(ty, year_to_band, sum))
  fb <- bpclm(m = mb, wide_breaks = wb, a = 0, b = 100,
              ngrid = 100L, ndx = 22L, degree = 3L, penalty_order = 3L,
              niter = 3000L, burnin = 800L, adapt = 400L,
              seed = 1000 + s, verbose = FALSE)
  fbc <- calibrate(fb)
  Ms <- nrow(fbc$pi_chain)
  pred <- matrix(NA_real_, nrow = Ms, ncol = length(ages))
  for (k in seq_len(Ms)) {
    pis <- fbc$pi_chain[k, ]
    for (j in seq_len(J)) {
      idx <- which(year_to_band == j)
      g_j <- sum(pis[idx])
      if (g_j > 0) {
        pred[k, idx] <- rmultinom(1, size = mb[j],
                                   prob = pis[idx] / g_j)[, 1L]
      } else {
        pred[k, idx] <- 0
      }
    }
  }
  lo <- apply(pred, 2L, quantile, probs = 0.05, names = FALSE)
  hi <- apply(pred, 2L, quantile, probs = 0.95, names = FALSE)
  cov_per_year[s, ] <- (ty >= lo) & (ty <= hi)
  rmse_per_rep[s]   <- sqrt(mean((colMeans(pred) - ty)^2))
  exact_rmse_per_rep[s] <-
    sqrt(mean((pclm_exact(m = mb, wide_breaks = wb, a = 0, b = 100,
                           ngrid = 100L, ndx = 22L, degree = 3L,
                           penalty_order = 3L)$pi *
                 sum(mb) - ty)^2))
}
cov_overall <- mean(cov_per_year)
cov_per_age <- colMeans(cov_per_year)
cat(sprintf("Average per-year 90%% coverage across replicates: %.1f%%\n",
            100 * cov_overall))
cat(sprintf("Replicate-level RMSE  bpclm+calibrate (median): %.2f\n",
            median(rmse_per_rep)))
cat(sprintf("Replicate-level RMSE  pclm_exact      (median): %.2f\n",
            median(exact_rmse_per_rep)))
