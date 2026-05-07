# Test posterior_predict() on bpclm and pclm_exact / pclm.
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())

ok <- function(label, cond) {
  cat(if (isTRUE(cond)) "OK    " else "FAIL  ", label, "\n", sep = "")
  invisible(isTRUE(cond))
}

set.seed(1)
ages <- 0:99
true_p <- 0.10 * dnorm(ages, 2, 3) +
          0.06 * dnorm(ages, 25, 12) +
          0.84 * dnorm(ages, 75, 12)
true_p <- true_p / sum(true_p)
N <- 50000L
true_yearly <- as.numeric(rmultinom(1, size = N, prob = true_p))
year_to_band <- (ages %/% 5L) + 1L
m_band <- as.numeric(tapply(true_yearly, year_to_band, sum))
wb <- cbind(seq(0, 95, 5), seq(5, 100, 5))

# === bpclm + calibrate + posterior_predict ===
fit_b <- bpclm(m = m_band, wide_breaks = wb, a = 0, b = 100,
               ngrid = 100L, ndx = 22L, degree = 3L, penalty_order = 3L,
               niter = 2000L, burnin = 500L, adapt = 300L, seed = 7)
fit_b <- calibrate(fit_b)

pp_pred <- posterior_predict(fit_b, type = "predictive", level = 0.9, seed = 7)
pp_rate <- posterior_predict(fit_b, type = "rate",       level = 0.9)

cat("\nposterior_predict(bpclm, predictive):\n")
print(pp_pred, n = 4)

cat("\nposterior_predict(bpclm, rate):\n")
print(pp_rate, n = 4)

# Check exact band-total preservation under predictive
band_check <- function(yearly) as.numeric(tapply(yearly, year_to_band, sum))
band_per_draw <- apply(pp_pred$draws, 1L, band_check)   # J x nsim
ok("predictive draws preserve band totals exactly",
   max(abs(band_per_draw - matrix(m_band, nrow = length(m_band),
                                  ncol = ncol(band_per_draw)))) == 0)

# Per-year coverage
cov <- (true_yearly >= pp_pred$lower) & (true_yearly <= pp_pred$upper)
cat(sprintf("Per-year predictive 90%% coverage: %d / %d (%.1f%%)\n",
            sum(cov), length(cov), 100 * mean(cov)))
ok("predictive coverage is reasonable",
   mean(cov) >= 0.75 && mean(cov) <= 0.95)

# === pclm_exact + posterior_predict ===
fit_e <- pclm_exact(m = m_band, wide_breaks = wb, a = 0, b = 100,
                    ngrid = 100L, ndx = 22L, degree = 3L, penalty_order = 3L)
pp_e  <- posterior_predict(fit_e, type = "predictive",
                            n_draws = 1000L, level = 0.9, seed = 7)
cat("\nposterior_predict(pclm_exact, predictive):\n")
print(pp_e, n = 4)
band_per_draw <- apply(pp_e$draws, 1L, band_check)
ok("pclm_exact predictive draws preserve band totals exactly",
   max(abs(band_per_draw - matrix(m_band, nrow = length(m_band),
                                  ncol = ncol(band_per_draw)))) == 0)

# === Plot ===
out_png <- file.path(pkg_root, "tests", "posterior_predict_demo.png")
png(out_png, width = 2400, height = 1200, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
plot(pp_pred, show_bins = FALSE,
     main = "bpclm + calibrate, type = 'predictive' (90% PI)",
     xlab = "Age (years)", ylab = "Deaths in this single-year")
lines(ages, true_yearly, col = "firebrick", lwd = 2, lty = 2)
legend("topleft",
       legend = c("Posterior predictive mean", "90% PI", "Truth"),
       col = c("black", grDevices::adjustcolor("steelblue", 0.6), "firebrick"),
       lwd = c(2, 8, 2), lty = c(1, 1, 2), bty = "n")

plot(pp_rate, show_bins = FALSE,
     main = "bpclm + calibrate, type = 'rate' (90% CI)",
     xlab = "Age (years)", ylab = "Expected deaths in this single-year")
lines(ages, true_yearly, col = "firebrick", lwd = 2, lty = 2)
legend("topleft",
       legend = c("Posterior mean rate", "90% CI", "Realised truth"),
       col = c("black", grDevices::adjustcolor("steelblue", 0.6), "firebrick"),
       lwd = c(2, 8, 2), lty = c(1, 1, 2), bty = "n")
dev.off()
cat("\nWrote: ", out_png, "\n", sep = "")
