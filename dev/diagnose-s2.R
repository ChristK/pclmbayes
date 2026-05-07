# Diagnostic for the S2 (lognormal with open last bin) fit.
# Variants tried:
#   A. Original setup: open bin (30, 80], shape = "unimodal"
#   B. Shape = NULL (no constraint), open bin (30, 80]
#   C. Tighter upper bound: open bin (30, 35], shape = "unimodal"
#   D. As B but with much smaller b (= 35)
#   E. As A but a longer chain
# We compare the freq + Bayesian fits and report observed vs expected
# bin counts, the warmstart freq fit, and the posterior-mean density.

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())

set.seed(42)

ml <- 2; sl <- 0.5
n2 <- 2000
y2 <- rlnorm(n2, meanlog = ml, sdlog = sl)

cat("Sample summary:\n")
cat("  n =", n2, "\n")
cat("  empirical mean =", round(mean(y2), 3),
    " | true mean =", round(exp(ml + sl^2/2), 3), "\n")
cat("  max(y) =", round(max(y2), 3), "\n")
cat("  q99 of data =", round(quantile(y2, 0.99), 3),
    " | true q99 =", round(qlnorm(0.99, ml, sl), 3), "\n")

# Build wide bins for variant A: (0,2], ..., (28,30], (30, 80]
bf      <- seq(0, 30, by = 2)
wb_A    <- cbind(c(head(bf, -1), 30), c(tail(bf, -1), 80))
m_A     <- integer(nrow(wb_A))
for (j in seq_len(nrow(wb_A))) {
  if (j == nrow(wb_A))
    m_A[j] <- sum(y2 >= wb_A[j, 1])
  else
    m_A[j] <- sum(y2 >= wb_A[j, 1] & y2 < wb_A[j, 2])
}
cat("\nA: bin counts (open bin extended to 80):\n")
print(data.frame(lower = wb_A[, 1], upper = wb_A[, 2], count = m_A))
cat("  total counts in last (open) bin =", tail(m_A, 1L), "\n")

# Variant C: tighter upper bound b = 35
wb_C    <- cbind(c(head(bf, -1), 30), c(tail(bf, -1), 35))
m_C     <- m_A   # identical counts; only the bin upper bound changes

# ---- Frequentist fits at each setup ----
fitA_f <- pclm(m = m_A, wide_breaks = wb_A, a = 0, b = 80,
               ngrid = 100L, ndx = 19L, degree = 3L, penalty_order = 3L)
fitC_f <- pclm(m = m_C, wide_breaks = wb_C, a = 0, b = 35,
               ngrid = 70L,  ndx = 17L, degree = 3L, penalty_order = 3L)
fitC2_f <- pclm(m = m_C, wide_breaks = wb_C, a = 0, b = 35,
               ngrid = 140L,  ndx = 27L, degree = 3L, penalty_order = 3L)

# Whether each freq fit's pi is unimodal (= would the unimodality
# fallback be triggered?)
cat("\nFreq fit unimodal?\n")
cat("  A (b=80) -> ", is_unimodal(fitA_f$pi), "  tau =", fitA_f$tau, "\n")
cat("  C (b=35) -> ", is_unimodal(fitC_f$pi), "  tau =", fitC_f$tau, "\n")
cat("  C* finer grid (b=35, ngrid=140) -> ",
    is_unimodal(fitC2_f$pi), "  tau =", fitC2_f$tau, "\n")

# ---- Bayesian fits ----
fitA <- bpclm(m = m_A, wide_breaks = wb_A, a = 0, b = 80,
              ngrid = 100L, ndx = 19L, degree = 3L, penalty_order = 3L,
              niter = 6000L, burnin = 1500L, adapt = 500L,
              shape = "unimodal", seed = 21)
fitB <- bpclm(m = m_A, wide_breaks = wb_A, a = 0, b = 80,
              ngrid = 100L, ndx = 19L, degree = 3L, penalty_order = 3L,
              niter = 6000L, burnin = 1500L, adapt = 500L,
              shape = NULL, seed = 22)
fitC <- bpclm(m = m_C, wide_breaks = wb_C, a = 0, b = 35,
              ngrid = 70L,  ndx = 17L, degree = 3L, penalty_order = 3L,
              niter = 6000L, burnin = 1500L, adapt = 500L,
              shape = "unimodal", seed = 23)
fitD <- bpclm(m = m_C, wide_breaks = wb_C, a = 0, b = 35,
              ngrid = 70L,  ndx = 17L, degree = 3L, penalty_order = 3L,
              niter = 6000L, burnin = 1500L, adapt = 500L,
              shape = NULL, seed = 24)

cat("\nAcceptance rates:\n")
cat("  A (b=80, unimodal)      acc = ", round(fitA$accept, 3), "\n")
cat("  B (b=80, no constraint) acc = ", round(fitB$accept, 3), "\n")
cat("  C (b=35, unimodal)      acc = ", round(fitC$accept, 3), "\n")
cat("  D (b=35, no constraint) acc = ", round(fitD$accept, 3), "\n")

# Posterior mean of mu_Y, sd_Y for each
report <- function(label, fit) {
  mu_chain <- apply(fit$pi_chain, 1L,
                    function(p) sum(fit$grid_mid * p))
  sd_chain <- apply(fit$pi_chain, 1L,
                    function(p) {
                      mu <- sum(fit$grid_mid * p)
                      sqrt(sum((fit$grid_mid - mu)^2 * p))
                    })
  cat(sprintf("  %s : mu = %.3f (CI %.3f, %.3f), sd = %.3f (CI %.3f, %.3f)\n",
              label,
              mean(mu_chain),
              quantile(mu_chain, 0.05, names = FALSE),
              quantile(mu_chain, 0.95, names = FALSE),
              mean(sd_chain),
              quantile(sd_chain, 0.05, names = FALSE),
              quantile(sd_chain, 0.95, names = FALSE)))
}
cat("\nPosterior summaries (true mean = ", round(exp(ml + sl^2/2), 3),
    ", true sd = ", round(sqrt((exp(sl^2)-1)*exp(2*ml+sl^2)), 3), ")\n", sep="")
report("A: b=80, unimodal     ", fitA)
report("B: b=80, no constraint", fitB)
report("C: b=35, unimodal     ", fitC)
report("D: b=35, no constraint", fitD)

# ---- Plot a comparison ----
out_png <- file.path(pkg_root, "tests", "diagnose_s2.png")
png(out_png, width = 2400, height = 1600, res = 200)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.5, 1.0), oma = c(0, 0, 2, 0))
plot_one <- function(fit, title, xlim = NULL) {
  delta <- diff(fit$grid)
  fitted <- fit$pi / delta
  lo <- fit$pi_lower / delta
  hi <- fit$pi_upper / delta
  if (is.null(xlim)) xlim <- range(fit$grid)
  bin_w  <- fit$wide_breaks[, 2] - fit$wide_breaks[, 1]
  bin_d  <- (fit$m / sum(fit$m)) / bin_w
  ymax   <- max(c(fitted, hi, dlnorm(fit$grid_mid, ml, sl), bin_d), na.rm = TRUE) * 1.05
  plot(NA, xlim = xlim, ylim = c(0, ymax),
       xlab = "y", ylab = "Density", main = title)
  for (j in seq_len(nrow(fit$wide_breaks))) {
    rect(fit$wide_breaks[j, 1], 0, fit$wide_breaks[j, 2], bin_d[j],
         col = grDevices::adjustcolor("grey80", 0.5), border = "grey50")
  }
  polygon(c(fit$grid_mid, rev(fit$grid_mid)), c(lo, rev(hi)),
          col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  lines(fit$grid_mid, fitted, col = "black", lwd = 2)
  curve(dlnorm(x, ml, sl),
        from = max(0.01, fit$grid[1L]), to = fit$grid[length(fit$grid)],
        add = TRUE, col = "firebrick", lwd = 2, lty = 2, n = 401)
}
plot_one(fitA, "A: open bin -> (30, 80], unimodal prior")
plot_one(fitB, "B: open bin -> (30, 80], NO constraint")
plot_one(fitC, "C: open bin -> (30, 35], unimodal prior", xlim = c(0, 35))
plot_one(fitD, "D: open bin -> (30, 35], NO constraint", xlim = c(0, 35))
mtext("Diagnostic for S2: effect of open-bin upper bound and unimodality constraint",
      side = 3, outer = TRUE, cex = 0.9, font = 2)
par(op)
dev.off()
cat("\nWrote: ", out_png, "\n", sep = "")
