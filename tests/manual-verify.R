# =============================================================================
# Manual verification script.
# Run from the package root with:
#   Rscript tests/manual-verify.R
#
# Verifies that:
#   1. The B-spline basis is a partition of unity.
#   2. The difference penalty has the expected rank and null space.
#   3. The bin matrix correctly handles aligned and partially overlapping bins.
#   4. The analytic gradient agrees with central finite differences.
#   5. The frequentist fit recovers a Gamma(5, 1) density from grouped counts.
#   6. The Bayesian fit on `bloodlead` reproduces the posterior summaries
#      published in Lambert and Eilers (2009), Section 6.1.
#   7. With shape = "unimodal", every kept draw of the latent density is
#      unimodal.
# =============================================================================

# --- locate package root (works from any cwd) ---------------------------------
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L])
  else if (sys.nframe() >= 1L && !is.null(sys.frame(1)$ofile))
    sys.frame(1)$ofile
  else NA_character_
})()
if (!is.na(script_path) && file.exists(script_path)) {
  pkg_root <- normalizePath(file.path(dirname(script_path), ".."),
                            mustWork = FALSE)
} else {
  pkg_root <- normalizePath(".", mustWork = FALSE)
}
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  stop("Cannot find package root from cwd '", getwd(),
       "' (script_path = '", script_path, "').")
}
cat("Package root: ", pkg_root, "\n", sep = "")

# Source every R file directly so we don't need to install the package.
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) {
  source(f, local = globalenv())
}

# Load datasets
load(file.path(pkg_root, "data", "bloodlead.rda"))
load(file.path(pkg_root, "data", "tbdeaths1907.rda"))

ok <- function(label, cond, info = NULL) {
  status <- if (isTRUE(cond)) "OK    " else "FAIL  "
  cat(status, label, "\n", sep = "")
  if (!isTRUE(cond) && !is.null(info)) cat("        ", info, "\n", sep = "")
  invisible(isTRUE(cond))
}
hr <- function(s) cat("\n---- ", s, " ", strrep("-", max(0, 60 - nchar(s))),
                      "\n", sep = "")

# 1 -------------------------------------------------------- B-spline basis ---
hr("(1) B-spline partition of unity")
mids <- seq(0.05, 9.95, by = 0.1)
bs   <- bspline_basis(mids, a = 0, b = 10, ndx = 17L, degree = 3L)
ok("dim(B) == (100, 20)", all(dim(bs$B) == c(100L, 20L)))
ok("rowSums(B) == 1 to 1e-9",
   max(abs(rowSums(bs$B) - 1)) < 1e-9,
   sprintf("max deviation = %.2e", max(abs(rowSums(bs$B) - 1))))

# 2 ----------------------------------------------- difference penalty ---
hr("(2) Difference penalty matrix")
P3 <- diff_penalty(K = 8L, r = 3L)
ok("rank(P3) == K - r = 5", qr(P3)$rank == 5L)
k <- 1:8
ok("constants in null space",     max(abs(P3 %*% rep(1, 8))) < 1e-10)
ok("linear seq in null space",    max(abs(P3 %*% k))         < 1e-10)
ok("quadratic seq in null space", max(abs(P3 %*% k^2))       < 1e-10)

# 3 ---------------------------------------------------------- bin matrix ---
hr("(3) Bin matrix")
fb <- seq(0, 10, by = 0.5)
wb <- seq(0, 10, by = 2)
C  <- bin_matrix(wb, fb)
ok("dim(C) == (5, 20)",         all(dim(C) == c(5L, 20L)))
# When wide bins are aligned partitions of the support, each fine interval
# is fully inside exactly one wide bin, so colSums(C) == 1.  rowSums(C)
# equals (wide bin width) / (fine interval width), here 2/0.5 = 4, which
# is correct because gamma_j = sum_i C[j,i] * pi_i and pi_i is a probability
# mass, not a density.
ok("rowSums(C) == 4 (aligned: 2/0.5)", max(abs(rowSums(C) - 4)) < 1e-12)
ok("colSums(C) == 1 (aligned)",        max(abs(colSums(C) - 1)) < 1e-12)
# Cross-check: gamma_j from a uniform pi should equal bin-width / total-width
pi_uni <- rep(1/20, 20)
ga_uni <- as.numeric(C %*% pi_uni)
ok("uniform pi -> gamma_j = bin_width / total_width",
   max(abs(ga_uni - 0.2)) < 1e-12)
fb2 <- seq(0, 1, length.out = 11)
C2  <- bin_matrix(matrix(c(0.05, 0.55), nrow = 1), fb2)
ok("partial overlap C2[1, 1]  == 0.5",  abs(C2[1, 1] - 0.5)  < 1e-12)
ok("partial overlap C2[1, 6]  == 0.5",  abs(C2[1, 6] - 0.5)  < 1e-12)
ok("interior cells == 1",
   max(abs(C2[1, 2:5] - 1)) < 1e-12)

# 4 ----------------------------------- analytic gradient vs finite diffs ---
hr("(4) Analytic gradient vs central finite differences")
set.seed(2009)
fb   <- seq(0, 10, length.out = 21)
mids <- (head(fb, -1) + tail(fb, -1)) / 2
bs   <- bspline_basis(mids, a = 0, b = 10, ndx = 5L, degree = 3L)
K    <- bs$K
wb   <- seq(0, 10, by = 2)
Cmat <- bin_matrix(wb, fb)
true_phi <- (1:K) - mean(1:K)
pi_true  <- .softmax(bs$B %*% true_phi)
ga_true  <- as.numeric(Cmat %*% pi_true)
m_obs    <- as.integer(round(500 * ga_true))

phi <- rnorm(K, sd = 0.3); phi <- phi - mean(phi)
ev  <- .pclm_eval(phi, B = bs$B, C = Cmat, m = m_obs, compute_FI = FALSE)

eps     <- 1e-5
grad_fd <- numeric(K)
for (kk in seq_len(K)) {
  dvec <- numeric(K); dvec[kk] <- 1; dvec <- dvec - mean(dvec)
  fp <- .pclm_eval(phi + eps * dvec, B = bs$B, C = Cmat, m = m_obs,
                   compute_FI = FALSE)$logL
  fm <- .pclm_eval(phi - eps * dvec, B = bs$B, C = Cmat, m = m_obs,
                   compute_FI = FALSE)$logL
  grad_fd[kk] <- (fp - fm) / (2 * eps)
}
grad_proj <- ev$grad - mean(ev$grad)
err <- max(abs(grad_proj - grad_fd))
cat(sprintf("    max |analytic - FD| = %.3e\n", err))
ok("analytic gradient ~= finite difference (< 1e-4)", err < 1e-4)

# 5 -------------------------- pclm() recovers Gamma(5, 1) density ---
hr("(5) Frequentist fit on simulated grouped Gamma(5, 1) data")
set.seed(1)
y    <- rgamma(5000, shape = 5, rate = 1)
brk  <- c(0, 2, 4, 6, 8, 10, 12, 16, 20)
m_w  <- as.numeric(table(cut(y, brk, include.lowest = TRUE)))
fit  <- pclm(m = m_w, wide_breaks = brk,
             a = 0, b = 20, ngrid = 80L,
             ndx = 17L, degree = 3L, penalty_order = 3L)
mu_hat <- sum(fit$grid_mid * fit$pi)
sd_hat <- sqrt(sum((fit$grid_mid - mu_hat)^2 * fit$pi))
cat(sprintf("    fitted mean = %.3f (target 5.000), sd = %.3f (target %.3f)\n",
            mu_hat, sd_hat, sqrt(5)))
ok("|fitted mean - 5| < 0.2",        abs(mu_hat - 5)        < 0.2)
ok("|fitted sd - sqrt(5)| < 0.3",    abs(sd_hat - sqrt(5))  < 0.3)
ok("fit converged",                  isTRUE(fit$converged))

# 6 --- bpclm() on bloodlead vs paper's published posterior summaries ---
hr("(6) Bayesian fit on bloodlead vs paper's published summaries")
set.seed(2009)
fitB <- bpclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80,
              ngrid = 80L, ndx = 17L, degree = 3L, penalty_order = 3L,
              niter = 4000L, burnin = 1000L, adapt = 500L,
              shape = "unimodal", seed = 2009)
acc <- fitB$accept
cat(sprintf("    acceptance rate = %.2f (target 0.57)\n", acc))
ok("accept in [0.30, 0.85]", acc >= 0.30 && acc <= 0.85,
   sprintf("acceptance = %.2f", acc))

mu_chain <- apply(fitB$pi_chain, 1L,
                  function(p) sum(fitB$grid_mid * p))
sd_chain <- apply(fitB$pi_chain, 1L,
                  function(p) {
                    mu <- sum(fitB$grid_mid * p)
                    sqrt(sum((fitB$grid_mid - mu)^2 * p))
                  })
prob_above_30_chain <- apply(fitB$pi_chain, 1L, function(p) {
  1 - approx(x = fitB$grid, y = c(0, cumsum(p)),
             xout = 30, rule = 2)$y
})
q20_chain <- apply(fitB$pi_chain, 1L,
                   function(p) approx(x = c(0, cumsum(p)),
                                      y = fitB$grid,
                                      xout = 0.20, ties = mean, rule = 2)$y)
q80_chain <- apply(fitB$pi_chain, 1L,
                   function(p) approx(x = c(0, cumsum(p)),
                                      y = fitB$grid,
                                      xout = 0.80, ties = mean, rule = 2)$y)

paper <- list(
  mu            = list(mean = 21.8, lo = 20.6, hi = 23.0),
  sd            = list(mean =  8.3, lo =  7.3, hi =  9.6),
  prob_above_30 = list(mean = 0.14, lo = 0.10, hi = 0.19),
  q20           = list(mean = 14.6, lo = 13.1, hi = 15.9),
  q80           = list(mean = 27.8, lo = 26.1, hi = 29.5)
)
report <- function(label, chain, targ) {
  m  <- mean(chain)
  lo <- quantile(chain, 0.05, names = FALSE)
  hi <- quantile(chain, 0.95, names = FALSE)
  cat(sprintf("    %-15s  ours: %6.3f  (%.3f, %.3f)   |   paper: %6.3f  (%.3f, %.3f)\n",
              label, m, lo, hi, targ$mean, targ$lo, targ$hi))
  half <- (targ$hi - targ$lo) / 2
  ok(sprintf("posterior mean of %s within 1.5 x paper CI half-width", label),
     abs(m - targ$mean) < 1.5 * half,
     sprintf("ours=%.3f, paper=%.3f", m, targ$mean))
}
report("mean(Y)",       mu_chain,            paper$mu)
report("sd(Y)",         sd_chain,            paper$sd)
report("Pr(Y > 30)",    prob_above_30_chain, paper$prob_above_30)
report("Q_0.20",        q20_chain,           paper$q20)
report("Q_0.80",        q80_chain,           paper$q80)

# 7 ---- shape constraint truly enforces unimodality on every draw ----
hr("(7) Unimodality of every kept draw")
uni <- apply(fitB$pi_chain, 1L, is_unimodal)
ok(sprintf("100%% of %d kept draws are unimodal",
           nrow(fitB$pi_chain)), all(uni))

cat("\n", strrep("=", 65), "\n",
    "Verification complete.\n", sep = "")
