# Exact preservation of wide-bin totals: compare three approaches.
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())
load(file.path(pkg_root, "data", "bloodlead.rda"))

ok <- function(label, cond, info = NULL) {
  cat(if (isTRUE(cond)) "OK    " else "FAIL  ", label, "\n", sep = "")
  if (!isTRUE(cond) && !is.null(info)) cat("        ", info, "\n", sep = "")
  invisible(isTRUE(cond))
}

cat("============================================================\n")
cat(" Bloodlead -- exact preservation of wide-bin totals\n")
cat("============================================================\n")

fit_norm <- pclm(m = bloodlead$count,
                 wide_breaks = with(bloodlead, cbind(lower, upper)),
                 a = 0, b = 80, ngrid = 80L, ndx = 17L, degree = 3L,
                 penalty_order = 3L)

fit_cal <- calibrate(fit_norm)

fit_exact <- pclm_exact(m = bloodlead$count,
                        wide_breaks = with(bloodlead, cbind(lower, upper)),
                        a = 0, b = 80, ngrid = 80L, ndx = 17L, degree = 3L,
                        penalty_order = 3L)

cmp <- data.frame(
  lower = bloodlead$lower, upper = bloodlead$upper,
  m_obs   = bloodlead$count,
  pclm    = round(fit_norm$fitted_counts,  3L),
  pclm_calibrated = round(fit_cal$fitted_counts,  3L),
  pclm_exact      = round(fit_exact$fitted_counts, 3L)
)
print(cmp, row.names = FALSE)

cat(sprintf("\nMax |fitted - m|:  pclm  = %.3f\n",
            max(abs(fit_norm$fitted_counts  - fit_norm$m))))
cat(sprintf("                   calib = %.3e\n",
            max(abs(fit_cal$fitted_counts   - fit_cal$m))))
cat(sprintf("                   exact = %.3e   (residual = %.2e)\n",
            max(abs(fit_exact$fitted_counts - fit_exact$m)),
            fit_exact$constraint_residual))

ok("calibrate() preserves totals to 1e-6",
   max(abs(fit_cal$fitted_counts - fit_cal$m)) < 1e-6)
ok("pclm_exact() preserves totals to 1e-6",
   max(abs(fit_exact$fitted_counts - fit_exact$m)) < 1e-6)
ok("pclm_exact() converged",
   isTRUE(fit_exact$converged))

# ---- Bayesian calibration ----
fit_b   <- bpclm(m = bloodlead$count,
                 wide_breaks = with(bloodlead, cbind(lower, upper)),
                 a = 0, b = 80, ngrid = 80L, ndx = 17L, degree = 3L,
                 penalty_order = 3L,
                 niter = 2000L, burnin = 500L, adapt = 300L,
                 shape = "unimodal", seed = 99)
fit_b_c <- calibrate(fit_b)
gamma_each_draw <- fit_b_c$pi_chain %*% t(fit_b_c$C)   # nsim x J
max_per_draw    <- apply(abs(gamma_each_draw -
                               matrix(bloodlead$count / sum(bloodlead$count),
                                      nrow = nrow(gamma_each_draw),
                                      ncol = length(bloodlead$count),
                                      byrow = TRUE)), 1L, max)
ok("bpclm calibrate(): every draw preserves bin totals",
   max(max_per_draw) < 1e-6,
   sprintf("max residual across draws = %.2e", max(max_per_draw)))

# ============================================================================
# Realistic ungrouping scenario: deaths in 5-year age bands -> single years
# ============================================================================
cat("\n============================================================\n")
cat(" Synthetic ungrouping: 5-year bands -> single years of age\n")
cat("============================================================\n")
set.seed(2026)
# Simulate 10000 single-year-of-age "death" counts with a bathtub shape
ages <- 0:99
true_density <- 0.10 * dnorm(ages, mean =  2, sd =  3) +     # infant peak
                0.10 * dnorm(ages, mean = 25, sd = 12) +     # young adults
                0.80 * dnorm(ages, mean = 75, sd = 12)       # main mortality
true_density <- true_density / sum(true_density)
N_total      <- 10000L
single_year_counts <- as.numeric(rmultinom(1, size = N_total,
                                            prob = true_density))
cat("simulated total deaths  =", N_total,
    " (sum of single-year counts)\n")

# Group into 5-year age bands
band_lower <- seq(0,  95, by = 5)
band_upper <- seq(5, 100, by = 5)
J <- length(band_lower)
m <- integer(J)
for (j in seq_len(J))
  m[j] <- sum(single_year_counts[band_lower[j]:(band_upper[j] - 1L) + 1L])
cat("total deaths in 5-year bands =", sum(m), "\n")
ok("5-year band totals sum to N_total", sum(m) == N_total)

# Fit each method
fit_g     <- pclm(m = m, wide_breaks = cbind(band_lower, band_upper),
                  a = 0, b = 100, ngrid = 100L,
                  ndx = 22L, degree = 3L, penalty_order = 3L)
fit_g_cal <- calibrate(fit_g)
fit_g_x   <- pclm_exact(m = m, wide_breaks = cbind(band_lower, band_upper),
                        a = 0, b = 100, ngrid = 100L,
                        ndx = 22L, degree = 3L, penalty_order = 3L)

# Imputed single-year counts under each method
imputed <- function(fit) round(fit$pi * sum(fit$m), 3L)
band_check <- function(yearly_counts) {
  out <- integer(J)
  for (j in seq_len(J))
    out[j] <- sum(yearly_counts[band_lower[j]:(band_upper[j] - 1L) + 1L])
  out
}

cat("\nMax |fitted_band - m|:\n")
cat(sprintf("  pclm                                  = %.3f\n",
            max(abs(fit_g$fitted_counts - m))))
cat(sprintf("  calibrate(pclm)                       = %.3e\n",
            max(abs(fit_g_cal$fitted_counts - m))))
cat(sprintf("  pclm_exact                            = %.3e\n",
            max(abs(fit_g_x$fitted_counts - m))))

# After rounding the imputed single-year counts, do they still sum to m?
ok("pclm_exact: rounded single-year counts sum to bands (within +-1)",
   max(abs(band_check(round(fit_g_x$pi * sum(m))) - m)) <= 1)
ok("calibrate(pclm): rounded single-year counts sum to bands (within +-1)",
   max(abs(band_check(round(fit_g_cal$pi * sum(m))) - m)) <= 1)

# Visualise
out_png <- file.path(pkg_root, "tests", "exact_totals.png")
png(out_png, width = 2400, height = 1500, res = 200)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.5, 1.0), oma = c(0, 0, 2, 0))

plot_fit <- function(fit, title, true_d = NULL) {
  delta  <- diff(fit$grid)
  fitted <- fit$pi / delta
  bin_d  <- (fit$m / sum(fit$m)) / (fit$wide_breaks[, 2] - fit$wide_breaks[, 1])
  ymax   <- max(c(fitted, bin_d, true_d), na.rm = TRUE) * 1.05
  plot(NA, xlim = range(fit$grid), ylim = c(0, ymax),
       xlab = "Age (years)", ylab = "Density", main = title)
  for (j in seq_len(nrow(fit$wide_breaks)))
    rect(fit$wide_breaks[j, 1], 0, fit$wide_breaks[j, 2], bin_d[j],
         col = grDevices::adjustcolor("grey80", 0.5), border = "grey50")
  if (!is.null(true_d))
    lines(ages + 0.5, true_d, col = "firebrick", lwd = 2, lty = 2)
  lines(fit$grid_mid, fitted, col = "black", lwd = 2)
}
plot_fit(fit_g,     "pclm (BIC) -- bands NOT preserved",
         true_d = true_density)
plot_fit(fit_g_cal, "calibrate(pclm) -- bands preserved",
         true_d = true_density)
plot_fit(fit_g_x,   "pclm_exact -- bands preserved (no kinks)",
         true_d = true_density)

# Per-bin discrepancy
plot(seq_len(J), fit_g$fitted_counts - m, type = "h", lwd = 6,
     col = "steelblue", xlab = "Age band",
     ylab = "fitted_count - observed",
     main = "Smoothing-induced shrinkage of band totals (pclm only)")
abline(h = 0, lty = 2, col = "grey50")
points(seq_len(J), fit_g_cal$fitted_counts - m, pch = 19, col = "firebrick")
points(seq_len(J), fit_g_x$fitted_counts   - m, pch = 4,  col = "darkgreen")
legend("topright",
       legend = c("pclm", "calibrate", "pclm_exact"),
       col = c("steelblue", "firebrick", "darkgreen"),
       lwd = c(6, NA, NA), pch = c(NA, 19, 4), bty = "n")
mtext("Exact preservation of wide-bin totals: three approaches",
      side = 3, outer = TRUE, cex = 0.95, font = 2)
par(op)
dev.off()
cat("\nWrote: ", out_png, "\n", sep = "")
