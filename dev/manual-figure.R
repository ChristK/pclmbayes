# Reproduces Fig. 4 of Lambert & Eilers (2009): bloodlead density with
# 90% credible band. Saves to outputs/bloodlead_fit.png.

pkg_root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)
                                              [grep("^--file=",
                                                    commandArgs(trailingOnly = FALSE))[1]]
                                            |> sub(pattern = "^--file=", replacement = "")),
                                    ".."))
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())
load(file.path(pkg_root, "data", "bloodlead.rda"))

set.seed(2009)
fit_f <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 80L, ndx = 17L, degree = 3L,
              penalty_order = 3L)
fit_b <- bpclm(m = bloodlead$count,
               wide_breaks = with(bloodlead, cbind(lower, upper)),
               a = 0, b = 80, ngrid = 80L, ndx = 17L, degree = 3L,
               penalty_order = 3L,
               niter = 5000L, burnin = 1000L, adapt = 500L,
               shape = "unimodal", seed = 2009)

out <- file.path(pkg_root, "tests", "bloodlead_fit.png")
png(out, width = 1100, height = 700, res = 150)
par(mar = c(4.2, 4.2, 2.0, 1.0))
plot(fit_b, xlab = expression(paste("Blood lead concentration  ",
                                      mu, "g/dl")),
     ylab = "Density",
     main = "Bayesian PCLM fit to the bloodlead data\n(reproducing Lambert & Eilers 2009, Fig. 4)")
lines(fit_f$grid_mid, fit_f$pi / diff(fit_f$grid),
      col = "firebrick", lwd = 2, lty = 2)
abline(v = 30, col = "darkgreen", lty = 3)
legend("topright",
       legend = c("Bayesian posterior mean", "Bayesian 90% credible band",
                  "Frequentist (BIC)", "30 μg/dl risk threshold"),
       col = c("black", "steelblue", "firebrick", "darkgreen"),
       lty = c(1, NA, 2, 3),
       pch = c(NA, 15, NA, NA),
       lwd = c(2, NA, 2, 1),
       bty = "n")
dev.off()
cat("Wrote: ", out, "\n", sep = "")
