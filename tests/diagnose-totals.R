# How exactly does the fit preserve wide-bin totals?
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())
load(file.path(pkg_root, "data", "bloodlead.rda"))

fit <- pclm(m = bloodlead$count,
            wide_breaks = with(bloodlead, cbind(lower, upper)),
            a = 0, b = 80, ngrid = 80L,
            ndx = 17L, degree = 3L, penalty_order = 3L)
cat("\nLeading-tau case (BIC selected): tau =", fit$tau, "\n")
cmp <- data.frame(lower = fit$wide_breaks[, 1],
                  upper = fit$wide_breaks[, 2],
                  m_obs = fit$m,
                  m_fit = round(fit$fitted_counts, 3),
                  diff  = round(fit$fitted_counts - fit$m, 3))
print(cmp, row.names = FALSE)
cat(sprintf("Total observed = %d, total fitted = %.3f, |diff|_max = %.3f\n",
            sum(fit$m), sum(fit$fitted_counts),
            max(abs(fit$fitted_counts - fit$m))))

# Try with vanishingly small tau (no smoothing): does it preserve totals?
fit0 <- pclm(m = bloodlead$count,
             wide_breaks = with(bloodlead, cbind(lower, upper)),
             a = 0, b = 80, ngrid = 80L,
             ndx = 17L, degree = 3L, penalty_order = 3L,
             tau = 1e-8)
cat("\ntau = 1e-8 (essentially unpenalised):\n")
cmp0 <- data.frame(lower = fit0$wide_breaks[, 1],
                   upper = fit0$wide_breaks[, 2],
                   m_obs = fit0$m,
                   m_fit = round(fit0$fitted_counts, 3),
                   diff  = round(fit0$fitted_counts - fit0$m, 3))
print(cmp0, row.names = FALSE)
cat(sprintf("|diff|_max = %.3f\n",
            max(abs(fit0$fitted_counts - fit0$m))))
