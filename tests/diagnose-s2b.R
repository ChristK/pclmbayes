# Re-run S2 using the SAME RNG state as the synthetic-eval script
# (i.e. consume the state advanced by S1 before drawing y2), and
# compare three chain lengths.

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())

# === Replicate synthetic-eval RNG state ============================
set.seed(42)
n1 <- 2000
y1 <- rgamma(n1, shape = 5, rate = 1)   # consumes RNG state
ml <- 2; sl <- 0.5
n2 <- 2000
y2 <- rlnorm(n2, meanlog = ml, sdlog = sl)
cat("Sample summary (synth-eval RNG state):\n")
cat("  empirical mean =", round(mean(y2), 3), "  true =",
    round(exp(ml + sl^2/2), 3), "\n")
cat("  empirical sd   =", round(sd(y2), 3), "  true =",
    round(sqrt((exp(sl^2)-1)*exp(2*ml+sl^2)), 3), "\n\n")

bf <- seq(0, 30, by = 2)
wb <- cbind(c(head(bf, -1), 30), c(tail(bf, -1), 80))
m  <- integer(nrow(wb))
for (j in seq_len(nrow(wb))) {
  if (j == nrow(wb)) m[j] <- sum(y2 >= wb[j, 1])
  else m[j] <- sum(y2 >= wb[j, 1] & y2 < wb[j, 2])
}
cat("bin counts:\n")
print(data.frame(lower = wb[, 1], upper = wb[, 2], count = m))

# Frequentist fit and unimodality check
fit_f <- pclm(m = m, wide_breaks = wb, a = 0, b = 80,
              ngrid = 100L, ndx = 19L, degree = 3L, penalty_order = 3L)
cat("\nFreq fit unimodal?", is_unimodal(fit_f$pi),
    " | tau =", fit_f$tau, " | edf =", round(fit_f$edf, 2),
    " | logL =", round(fit_f$logL, 2), "\n")

# Three chain lengths, all using the synth-eval seed = 2
fits <- list()
for (cfg in list(list(niter = 4000L,  burnin = 1000L, label = "short  4k"),
                 list(niter = 8000L,  burnin = 2000L, label = "med    8k"),
                 list(niter = 16000L, burnin = 4000L, label = "long  16k"))) {
  cat("\nFitting:", cfg$label, "...\n")
  t0 <- Sys.time()
  fit <- bpclm(m = m, wide_breaks = wb, a = 0, b = 80,
               ngrid = 100L, ndx = 19L, degree = 3L, penalty_order = 3L,
               niter = cfg$niter, burnin = cfg$burnin, adapt = 500L,
               shape = "unimodal", seed = 2)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fits[[cfg$label]] <- fit
  mu_chain <- apply(fit$pi_chain, 1L, function(p) sum(fit$grid_mid * p))
  pk <- max(fit$pi / diff(fit$grid))   # posterior-mean peak height
  cat(sprintf("  acc = %.2f, t = %.1fs, mu = %.3f (CI %.3f, %.3f), peak = %.4f\n",
              fit$accept, dt,
              mean(mu_chain),
              quantile(mu_chain, 0.05, names = FALSE),
              quantile(mu_chain, 0.95, names = FALSE),
              pk))
}

out_png <- file.path(pkg_root, "tests", "diagnose_s2b.png")
png(out_png, width = 2400, height = 800, res = 200)
op <- par(mfrow = c(1, 3), mar = c(4.0, 4.0, 2.5, 1.0), oma = c(0, 0, 2, 0))
for (label in names(fits)) {
  fit <- fits[[label]]
  delta <- diff(fit$grid)
  fitted <- fit$pi / delta
  lo <- fit$pi_lower / delta
  hi <- fit$pi_upper / delta
  bin_w  <- fit$wide_breaks[, 2] - fit$wide_breaks[, 1]
  bin_d  <- (fit$m / sum(fit$m)) / bin_w
  ymax   <- 0.16
  plot(NA, xlim = c(0, 35), ylim = c(0, ymax),
       xlab = "y", ylab = "Density",
       main = sprintf("S2 (synth-eval sample) | %s | acc=%.2f",
                      label, fit$accept))
  for (j in seq_len(nrow(fit$wide_breaks))) {
    rect(fit$wide_breaks[j, 1], 0, fit$wide_breaks[j, 2], bin_d[j],
         col = grDevices::adjustcolor("grey80", 0.5), border = "grey50")
  }
  polygon(c(fit$grid_mid, rev(fit$grid_mid)), c(lo, rev(hi)),
          col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  lines(fit$grid_mid, fitted, col = "black", lwd = 2)
  curve(dlnorm(x, ml, sl), from = 0.01, to = 35,
        add = TRUE, col = "firebrick", lwd = 2, lty = 2, n = 401)
}
mtext("Effect of chain length on S2 with the synth-eval RNG sample",
      side = 3, outer = TRUE, cex = 0.9, font = 2)
par(op)
dev.off()
cat("\nWrote:", out_png, "\n")
