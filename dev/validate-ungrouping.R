# =============================================================================
# Validation of single-year ungrouping from 5-year bands.
#
# Procedure:
#   1. Generate single-year-of-age "death" counts on ages 0..99 from a
#      realistic bathtub mortality density (multinomial draw, total N).
#   2. Group these counts into 5-year age bands (0,5], (5,10], ...
#   3. Re-ungroup with three methods:
#         pclm       (default BIC fit, not exact)
#         calibrate  (post-hoc, exact bin totals)
#         pclm_exact (constrained MAP, exact bin totals)
#   4. For each method, compare to the TRUE single-year counts:
#         per-year residuals
#         RMS error per year
#         maximum absolute error per year
#         cumulative-sum CDF error
#         L1, L2, L_inf norms of (estimate - truth)
#   5. Repeat the entire procedure S times to characterise sampling
#      variability.  Report distribution of errors across replicates.
#
# A second pass uses a SMOOTH true density that matches the algorithm's
# B-spline class so we can see the irreducible recovery error from
# grouping alone, isolated from the model-mismatch error.
# =============================================================================

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())

# ---- Setup ---------------------------------------------------------------
ages         <- 0:99
band_lower   <- seq(0,  95, by = 5)
band_upper   <- seq(5, 100, by = 5)
J            <- length(band_lower)            # 20 bands
wb           <- cbind(band_lower, band_upper)

# Realistic bathtub mortality probability per single year of age.
true_p <- 0.10 * dnorm(ages, mean =  2,  sd =  3) +    # infant peak
          0.06 * dnorm(ages, mean = 25,  sd = 12) +    # young-adult peak
          0.84 * dnorm(ages, mean = 75,  sd = 12)      # main mortality peak
true_p <- true_p / sum(true_p)

# Bin assignment: which 5-year band does each single year of age belong to?
year_to_band <- (ages %/% 5L) + 1L              # 1..J
band_truth_p <- as.numeric(tapply(true_p, year_to_band, sum))

# Helper: ungroup with each method and return single-year counts
ungroup_one <- function(m, ngrid = 100L, ndx = 22L) {
  fit_default <- pclm(m = m, wide_breaks = wb, a = 0, b = 100,
                      ngrid = ngrid, ndx = ndx, degree = 3L,
                      penalty_order = 3L)
  fit_cal     <- calibrate(fit_default)
  fit_exact   <- pclm_exact(m = m, wide_breaks = wb, a = 0, b = 100,
                            ngrid = ngrid, ndx = ndx, degree = 3L,
                            penalty_order = 3L)

  # Convert pi (length ngrid) to single-year counts (length 100).
  # When ngrid = 100 and a=0, b=100, the fine grid IS single years.
  total <- sum(m)
  list(
    pclm        = fit_default$pi * total,
    calibrate   = fit_cal$pi * total,
    pclm_exact  = fit_exact$pi * total,
    fits        = list(default = fit_default, cal = fit_cal, exact = fit_exact)
  )
}

# Helper: compute summary error metrics
err_metrics <- function(est, truth) {
  list(
    L1     = sum(abs(est - truth)),
    L2     = sqrt(sum((est - truth) ^ 2)),
    Linf   = max(abs(est - truth)),
    rmse   = sqrt(mean((est - truth) ^ 2)),
    bias   = mean(est - truth),
    cdf_Linf = max(abs(cumsum(est) - cumsum(truth))),
    cdf_L2   = sqrt(mean((cumsum(est) - cumsum(truth)) ^ 2))
  )
}

# =============================================================================
# Single illustrative replicate
# =============================================================================
cat("============================================================\n")
cat(" Single illustrative replicate (N = 50,000)\n")
cat("============================================================\n")
set.seed(2026)
N <- 50000L
true_yearly <- as.numeric(rmultinom(1, size = N, prob = true_p))
m_band      <- as.numeric(tapply(true_yearly, year_to_band, sum))
cat("True N =", N, ",  band totals sum =", sum(m_band), "\n")

ug <- ungroup_one(m_band)

# Per-method error summary
methods <- c("pclm", "calibrate", "pclm_exact")
err_tbl <- do.call(rbind, lapply(methods, function(mn) {
  m_err <- err_metrics(ug[[mn]], true_yearly)
  data.frame(method = mn,
             L1   = round(m_err$L1, 1),
             L2   = round(m_err$L2, 2),
             Linf = round(m_err$Linf, 2),
             rmse = round(m_err$rmse, 3),
             bias = signif(m_err$bias, 3),
             cdf_Linf = round(m_err$cdf_Linf, 2),
             cdf_L2   = round(m_err$cdf_L2,   3),
             stringsAsFactors = FALSE)
}))
cat("\nError metrics for the single replicate:\n")
print(err_tbl, row.names = FALSE)

# Were the band totals preserved?
band_check <- function(yearly) as.numeric(tapply(yearly, year_to_band, sum))
cat("\nBand-total preservation (max |fitted_band - observed_band|):\n")
for (mn in methods)
  cat(sprintf("  %-12s = %.3e\n", mn,
              max(abs(band_check(ug[[mn]]) - m_band))))

# =============================================================================
# Replicates: characterise sampling variability of the recovery error
# =============================================================================
cat("\n============================================================\n")
cat(" Monte Carlo over 100 replicates (N = 50,000 each)\n")
cat("============================================================\n")
S <- 100L
collect <- list(pclm = list(), calibrate = list(), pclm_exact = list())
for (s in seq_len(S)) {
  set.seed(2026 + s)
  ty <- as.numeric(rmultinom(1, size = N, prob = true_p))
  mb <- as.numeric(tapply(ty, year_to_band, sum))
  ug_s <- ungroup_one(mb)
  for (mn in methods) collect[[mn]][[s]] <- err_metrics(ug_s[[mn]], ty)
}

agg <- function(field, mn) {
  vals <- vapply(collect[[mn]], `[[`, FUN.VALUE = numeric(1L), field)
  vals
}
mc_summary <- do.call(rbind, lapply(methods, function(mn) {
  data.frame(
    method   = mn,
    rmse_med = round(median(agg("rmse",     mn)),  3),
    Linf_med = round(median(agg("Linf",     mn)),  2),
    cdfLinf_med = round(median(agg("cdf_Linf", mn)), 2),
    cdfL2_med   = round(median(agg("cdf_L2",   mn)), 3),
    stringsAsFactors = FALSE
  )
}))
cat("\nMonte Carlo medians over", S, "replicates:\n")
print(mc_summary, row.names = FALSE)

# The "oracle" reference: the irreducible recovery error from grouping
# alone equals the multinomial sampling noise *within* each band.  For a
# 5-year band of expected count n_j and uniform within-band p, the per-
# year sampling noise is roughly sqrt(n_j * (1/5) * (4/5)) ~ sqrt(n_j / 6).
# That's the floor any ungrouping method has to live with.
cat(sprintf("\nReference: typical band has E[n_j] ~ %.0f deaths, ",
            mean(N * band_truth_p)))
cat(sprintf("so within-band single-year sampling SD ~ %.1f deaths/year.\n",
            sqrt(mean(N * band_truth_p) / 6)))

# =============================================================================
# Plot: single illustrative replicate
# =============================================================================
out_png <- file.path(pkg_root, "tests", "validate_ungrouping.png")
png(out_png, width = 2600, height = 1700, res = 200)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.5, 1.0), oma = c(0, 0, 2, 0))

# Panel 1: per-year truth vs recovered counts
ymax <- max(c(true_yearly, ug$pclm, ug$calibrate, ug$pclm_exact)) * 1.05
plot(ages, true_yearly, type = "h", col = "grey60", lwd = 5,
     xlab = "Age (years)", ylab = "Deaths in this single-year",
     main = "Single-year counts: truth (grey) vs ungrouped",
     ylim = c(0, ymax))
lines(ages, ug$pclm,       col = "steelblue", lwd = 2)
lines(ages, ug$calibrate,  col = "firebrick", lwd = 2, lty = 1)
lines(ages, ug$pclm_exact, col = "darkgreen", lwd = 2, lty = 2)
legend("topleft",
       legend = c("Truth (single-year counts)",
                  "pclm (BIC, no constraint)",
                  "calibrate(pclm)",
                  "pclm_exact"),
       col = c("grey60", "steelblue", "firebrick", "darkgreen"),
       lwd = c(5, 2, 2, 2), lty = c(1, 1, 1, 2), bty = "n")

# Panel 2: per-year residuals
res_pclm  <- ug$pclm       - true_yearly
res_cal   <- ug$calibrate  - true_yearly
res_exact <- ug$pclm_exact - true_yearly
yrange <- range(c(res_pclm, res_cal, res_exact))
plot(ages, res_pclm, type = "l", col = "steelblue", lwd = 2,
     xlab = "Age (years)", ylab = "estimated - true",
     main = "Per-year residuals", ylim = yrange)
lines(ages, res_cal,   col = "firebrick", lwd = 2)
lines(ages, res_exact, col = "darkgreen", lwd = 2, lty = 2)
abline(h = 0, lty = 3, col = "grey50")
abline(v = seq(5, 95, by = 5), lty = 3, col = "grey85")
legend("topright",
       legend = c("pclm", "calibrate", "pclm_exact"),
       col = c("steelblue", "firebrick", "darkgreen"),
       lwd = 2, lty = c(1, 1, 2), bty = "n")

# Panel 3: cumulative-sum residual at every age boundary
cdf_t  <- cumsum(true_yearly)
cdf_p  <- cumsum(ug$pclm)
cdf_c  <- cumsum(ug$calibrate)
cdf_e  <- cumsum(ug$pclm_exact)
yrange <- range(c(cdf_p - cdf_t, cdf_c - cdf_t, cdf_e - cdf_t))
plot(ages, cdf_p - cdf_t, type = "l", col = "steelblue", lwd = 2,
     xlab = "Age (years)", ylab = "Cumulative deaths: estimated - true",
     main = "CDF residual",
     ylim = yrange)
lines(ages, cdf_c - cdf_t, col = "firebrick", lwd = 2)
lines(ages, cdf_e - cdf_t, col = "darkgreen", lwd = 2, lty = 2)
abline(h = 0, lty = 3, col = "grey50")
abline(v = seq(5, 95, by = 5), lty = 3, col = "grey85")
legend("topleft",
       legend = c("pclm", "calibrate", "pclm_exact",
                  "5-year band boundaries"),
       col = c("steelblue", "firebrick", "darkgreen", "grey85"),
       lwd = c(2, 2, 2, 1), lty = c(1, 1, 2, 3), bty = "n")

# Panel 4: distribution of RMSE across replicates
rmse_pclm  <- agg("rmse", "pclm")
rmse_cal   <- agg("rmse", "calibrate")
rmse_exact <- agg("rmse", "pclm_exact")
boxplot(list(pclm = rmse_pclm, calibrate = rmse_cal, pclm_exact = rmse_exact),
        ylab = "RMSE per single year (deaths)",
        main = sprintf("Distribution of RMSE across %d Monte Carlo replicates",
                       S),
        col = grDevices::adjustcolor(c("steelblue", "firebrick", "darkgreen"),
                                      alpha.f = 0.3),
        border = c("steelblue", "firebrick", "darkgreen"))
abline(h = sqrt(mean(N * band_truth_p) / 6), lty = 2, col = "grey40")
mtext(sprintf("Within-band sampling SD ~ %.1f deaths/year",
              sqrt(mean(N * band_truth_p) / 6)),
      side = 3, line = -1.4, cex = 0.7, col = "grey40", adj = 0.98)

mtext("Validation of single-year ungrouping from 5-year mortality bands",
      side = 3, outer = TRUE, cex = 0.95, font = 2)
par(op)
dev.off()

cat("\nWrote: ", out_png, "\n", sep = "")
