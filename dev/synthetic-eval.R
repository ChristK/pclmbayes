# =============================================================================
# Performance evaluation on synthetic data.
# Run from the package root with:
#   Rscript dev/synthetic-eval.R
#
# Four scenarios:
#   (S1) Unimodal Gamma(5, 1), equal-width bins (width = 0.5*sigma).
#   (S2) Unimodal heavy-tailed Lognormal(meanlog = 2, sdlog = 0.5), equal-
#        width bins with the last bin "open" (e.g. "30+").  We replace
#        the open bin by a finite upper bound (here b = 80) following
#        Lambert and Eilers (2009, Section 4).
#   (S3) Bimodal mixture 0.5 N(2, 0.5^2) + 0.5 N(7, 0.7^2), UNEQUAL-width
#        bins.  No shape constraint -- bimodality should be recovered.
#   (S4) Trimodal mixture 0.4 N(1, 0.5^2) + 0.3 N(5, 0.6^2) +
#        0.3 N(9, 0.8^2), equal-width bins.
#
# Metrics for each scenario:
#   * acceptance rate of the MCMC step
#   * Integrated squared error (ISE) between true and posterior-mean density
#   * |bias| in mean(Y), sd(Y), median(Y) and 0.05 / 0.95 quantiles
#   * coverage indicator for the 90% credible interval of mean(Y) and sd(Y)
#
# Output:
#   dev/synthetic_eval.png      -- 2x2 panel of true vs fitted densities
#   dev/synthetic_eval.csv      -- numerical summary table
# =============================================================================

# ---- locate package root ----
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
if (!is.na(script_path) && file.exists(script_path)) {
  pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  pkg_root <- normalizePath(".", mustWork = FALSE)
}
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) {
  source(f, local = globalenv())
}

# ---- helpers ----
ise <- function(true_dens, fitted_dens, grid_mids, delta) {
  # Riemann approximation of the integrated squared error
  sum((true_dens - fitted_dens) ^ 2 * delta)
}

# Truncate-and-renormalise a true density on the support [a, b], returning
# its values at the fine-grid midpoints and the area under it (which should
# be ~ 1 if the support contains essentially all of the mass).
true_density_on_grid <- function(dfun, a, b, grid_mid, delta, ...) {
  d <- dfun(grid_mid, ...)
  d[grid_mid < a | grid_mid > b] <- 0
  d
}

# Quantile function of a discrete grid density via inverse CDF
qgrid <- function(p, dens_at_mid, breaks) {
  delta <- diff(breaks)
  pi_grid <- dens_at_mid * delta
  pi_grid <- pi_grid / sum(pi_grid)        # normalise (in case of mass loss)
  Fg      <- c(0, cumsum(pi_grid))
  approx(x = Fg, y = breaks, xout = p, ties = mean, rule = 2)$y
}

# Group continuous observations into wide bins; if the right-tail bin is
# "open", everything above the largest finite upper bound is added to it.
# Returns counts m and a J x 2 wide-bin limit matrix.
group_obs <- function(y, breaks_lower, breaks_upper, last_open = FALSE) {
  J <- length(breaks_lower)
  m <- integer(J)
  for (j in seq_len(J)) {
    if (j == J && last_open) {
      m[j] <- sum(y >= breaks_lower[j])
    } else {
      m[j] <- sum(y >= breaks_lower[j] & y < breaks_upper[j])
    }
  }
  m
}

# Fit summary that compares the posterior to a known truth
summarise_fit <- function(fit, true_mean, true_sd, true_q05, true_q50,
                          true_q95, true_dens_at_grid_mid, label) {
  delta <- diff(fit$grid)
  pi_mean <- fit$pi
  fitted_dens <- pi_mean / delta
  ise_val <- ise(true_dens_at_grid_mid, fitted_dens, fit$grid_mid, delta)

  # Posterior of mean and sd
  mu_chain <- apply(fit$pi_chain, 1L,
                    function(p) sum(fit$grid_mid * p))
  sd_chain <- apply(fit$pi_chain, 1L,
                    function(p) {
                      mu <- sum(fit$grid_mid * p)
                      sqrt(sum((fit$grid_mid - mu) ^ 2 * p))
                    })
  q05_chain <- apply(fit$pi_chain, 1L,
                     function(p) approx(x = c(0, cumsum(p)),
                                        y = fit$grid, xout = 0.05,
                                        ties = mean, rule = 2)$y)
  q50_chain <- apply(fit$pi_chain, 1L,
                     function(p) approx(x = c(0, cumsum(p)),
                                        y = fit$grid, xout = 0.50,
                                        ties = mean, rule = 2)$y)
  q95_chain <- apply(fit$pi_chain, 1L,
                     function(p) approx(x = c(0, cumsum(p)),
                                        y = fit$grid, xout = 0.95,
                                        ties = mean, rule = 2)$y)
  cover_mu <- (true_mean >= quantile(mu_chain, 0.05) &
               true_mean <= quantile(mu_chain, 0.95))
  cover_sd <- (true_sd   >= quantile(sd_chain, 0.05) &
               true_sd   <= quantile(sd_chain, 0.95))

  data.frame(
    scenario      = label,
    accept        = round(fit$accept,                 3L),
    ISE           = signif(ise_val,                   3L),
    bias_mean     = signif(mean(mu_chain) - true_mean, 3L),
    bias_sd       = signif(mean(sd_chain) - true_sd,   3L),
    bias_median   = signif(mean(q50_chain) - true_q50, 3L),
    bias_q05      = signif(mean(q05_chain) - true_q05, 3L),
    bias_q95      = signif(mean(q95_chain) - true_q95, 3L),
    cover_mean_90 = cover_mu,
    cover_sd_90   = cover_sd,
    stringsAsFactors = FALSE
  )
}

set.seed(42)

# =============================================================================
# Scenario 1: Unimodal Gamma(5, 1), equal-width bins
# =============================================================================
n1   <- 2000
y1   <- rgamma(n1, shape = 5, rate = 1)
sig1 <- sqrt(5)
brk1 <- seq(0, 16, by = sig1 * 0.5)              # equal-width 0.5 sigma bins
lo1  <- head(brk1, -1); up1 <- tail(brk1, -1)
m1   <- group_obs(y1, lo1, up1, last_open = TRUE)
m1[length(m1)] <- m1[length(m1)] + sum(y1 >= max(brk1))
wb1  <- cbind(lo1, up1)
fit1 <- bpclm(m = m1, wide_breaks = wb1,
              a = 0, b = max(brk1), ngrid = 80L,
              ndx = 17L, degree = 3L, penalty_order = 3L,
              niter = 4000L, burnin = 1000L, adapt = 500L,
              shape = "unimodal", seed = 1)
true_d1 <- function(x) dgamma(x, shape = 5, rate = 1)
truth1 <- list(
  mean = 5,                  sd  = sqrt(5),
  q05  = qgamma(0.05, 5, 1),
  q50  = qgamma(0.50, 5, 1),
  q95  = qgamma(0.95, 5, 1)
)
td1 <- true_density_on_grid(true_d1, 0, max(brk1), fit1$grid_mid,
                            diff(fit1$grid))
sum1 <- summarise_fit(fit1, truth1$mean, truth1$sd,
                      truth1$q05, truth1$q50, truth1$q95,
                      td1, "S1: Gamma(5,1) | equal | unimodal prior")

# =============================================================================
# Scenario 2: Lognormal (heavy right tail), open last bin
# =============================================================================
ml <- 2; sl <- 0.5
n2 <- 2000
y2 <- rlnorm(n2, meanlog = ml, sdlog = sl)
# Wide bins (0,2], (2,4], ..., (28,30], 30+
breaks_finite <- seq(0, 30, by = 2)
lo2 <- breaks_finite                              # 0, 2, ..., 30
# treat the last bin as "30+" and replace by (30, 80]
wb2 <- cbind(c(head(lo2, -1), 30), c(tail(lo2, -1), 80))
m2  <- group_obs(y2, wb2[, 1], wb2[, 2], last_open = TRUE)
fit2 <- bpclm(m = m2, wide_breaks = wb2,
              a = 0, b = 80, ngrid = 100L,
              ndx = 19L, degree = 3L, penalty_order = 3L,
              niter = 8000L, burnin = 2000L, adapt = 800L,
              shape = "unimodal", seed = 2)
true_d2 <- function(x) dlnorm(x, meanlog = ml, sdlog = sl)
truth2 <- list(
  mean = exp(ml + sl^2 / 2),
  sd   = sqrt((exp(sl^2) - 1) * exp(2 * ml + sl^2)),
  q05  = qlnorm(0.05, ml, sl),
  q50  = qlnorm(0.50, ml, sl),
  q95  = qlnorm(0.95, ml, sl)
)
td2 <- true_density_on_grid(true_d2, 0, 80, fit2$grid_mid, diff(fit2$grid))
sum2 <- summarise_fit(fit2, truth2$mean, truth2$sd,
                      truth2$q05, truth2$q50, truth2$q95,
                      td2, "S2: Lognormal | equal | OPEN last bin (30+)")

# =============================================================================
# Scenario 3: Bimodal mixture, UNEQUAL-width bins, no shape constraint
# =============================================================================
mu_m <- c(2, 7); sd_m <- c(0.5, 0.7); w_m <- c(0.5, 0.5)
n3 <- 4000
component <- sample(seq_along(w_m), n3, replace = TRUE, prob = w_m)
y3 <- rnorm(n3, mean = mu_m[component], sd = sd_m[component])
y3 <- y3[y3 >= 0 & y3 <= 12]                     # truncate to support
# Unequal-width bins: narrow around the modes, wide where density is low
wb3 <- rbind(c(0,    1.0),  c(1.0,  1.5),  c(1.5, 2.0),
             c(2.0,  2.5),  c(2.5,  3.5),  c(3.5, 5.0),
             c(5.0,  6.0),  c(6.0,  6.5),  c(6.5, 7.0),
             c(7.0,  7.5),  c(7.5,  8.5),  c(8.5, 12.0))
m3 <- group_obs(y3, wb3[, 1], wb3[, 2], last_open = FALSE)
fit3 <- bpclm(m = m3, wide_breaks = wb3,
              a = 0, b = 12, ngrid = 120L,
              ndx = 23L, degree = 3L, penalty_order = 3L,
              niter = 5000L, burnin = 1500L, adapt = 600L,
              shape = NULL, seed = 3)         # NO unimodality constraint!
true_d3 <- function(x) {
  w_m[1] * dnorm(x, mu_m[1], sd_m[1]) + w_m[2] * dnorm(x, mu_m[2], sd_m[2])
}
truth3 <- list(
  mean = sum(w_m * mu_m),
  sd   = sqrt(sum(w_m * (sd_m^2 + mu_m^2)) - (sum(w_m * mu_m))^2)
)
xx <- seq(0, 12, length.out = 5001); dx <- xx[2] - xx[1]
F_xx <- cumsum(true_d3(xx)) * dx
truth3$q05 <- approx(F_xx / max(F_xx), xx, 0.05)$y
truth3$q50 <- approx(F_xx / max(F_xx), xx, 0.50)$y
truth3$q95 <- approx(F_xx / max(F_xx), xx, 0.95)$y
td3 <- true_density_on_grid(true_d3, 0, 12, fit3$grid_mid, diff(fit3$grid))
sum3 <- summarise_fit(fit3, truth3$mean, truth3$sd,
                      truth3$q05, truth3$q50, truth3$q95,
                      td3, "S3: Bimodal | UNequal | unconstrained")

# =============================================================================
# Scenario 4: Trimodal mixture, equal-width bins, unconstrained
# =============================================================================
mu_t <- c(1, 5, 9); sd_t <- c(0.5, 0.6, 0.8); w_t <- c(0.4, 0.3, 0.3)
n4 <- 4000
comp4 <- sample(seq_along(w_t), n4, replace = TRUE, prob = w_t)
y4 <- rnorm(n4, mean = mu_t[comp4], sd = sd_t[comp4])
y4 <- y4[y4 >= 0 & y4 <= 14]
brk4 <- seq(0, 14, by = 0.5)
wb4 <- cbind(head(brk4, -1), tail(brk4, -1))
m4  <- group_obs(y4, wb4[, 1], wb4[, 2], last_open = FALSE)
fit4 <- bpclm(m = m4, wide_breaks = wb4,
              a = 0, b = 14, ngrid = 140L,
              ndx = 27L, degree = 3L, penalty_order = 3L,
              niter = 5000L, burnin = 1500L, adapt = 600L,
              shape = NULL, seed = 4)
true_d4 <- function(x) {
  w_t[1] * dnorm(x, mu_t[1], sd_t[1]) +
  w_t[2] * dnorm(x, mu_t[2], sd_t[2]) +
  w_t[3] * dnorm(x, mu_t[3], sd_t[3])
}
truth4 <- list(
  mean = sum(w_t * mu_t),
  sd   = sqrt(sum(w_t * (sd_t^2 + mu_t^2)) - (sum(w_t * mu_t))^2)
)
xx <- seq(0, 14, length.out = 5001); dx <- xx[2] - xx[1]
F_xx <- cumsum(true_d4(xx)) * dx
truth4$q05 <- approx(F_xx / max(F_xx), xx, 0.05)$y
truth4$q50 <- approx(F_xx / max(F_xx), xx, 0.50)$y
truth4$q95 <- approx(F_xx / max(F_xx), xx, 0.95)$y
td4 <- true_density_on_grid(true_d4, 0, 14, fit4$grid_mid, diff(fit4$grid))
sum4 <- summarise_fit(fit4, truth4$mean, truth4$sd,
                      truth4$q05, truth4$q50, truth4$q95,
                      td4, "S4: Trimodal | equal | unconstrained")

# ---- Plot ----
out_png <- file.path(pkg_root, "tests", "synthetic_eval.png")
png(out_png, width = 2400, height = 1600, res = 200)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.5, 1.0), oma = c(0, 0, 2, 0))
plot_one <- function(fit, true_dfun, title, xlab, xlim = NULL) {
  delta <- diff(fit$grid)
  fitted <- fit$pi / delta
  lo <- fit$pi_lower / delta
  hi <- fit$pi_upper / delta
  bin_w  <- fit$wide_breaks[, 2] - fit$wide_breaks[, 1]
  bin_d  <- (fit$m / sum(fit$m)) / bin_w
  if (is.null(xlim)) xlim <- range(fit$grid)
  ymax   <- max(c(fitted, hi, true_dfun(fit$grid_mid),
                  bin_d), na.rm = TRUE) * 1.05
  plot(NA, xlim = xlim, ylim = c(0, ymax),
       xlab = xlab, ylab = "Density", main = title)
  for (j in seq_len(nrow(fit$wide_breaks))) {
    rect(fit$wide_breaks[j, 1], 0,
         fit$wide_breaks[j, 2], bin_d[j],
         col = grDevices::adjustcolor("grey80", 0.5), border = "grey50")
  }
  polygon(c(fit$grid_mid, rev(fit$grid_mid)), c(lo, rev(hi)),
          col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  lines(fit$grid_mid, fitted, col = "black", lwd = 2)
  curve(true_dfun, from = fit$grid[1L],
        to = fit$grid[length(fit$grid)],
        add = TRUE, col = "firebrick", lwd = 2, lty = 2, n = 401)
}
plot_one(fit1, true_d1,
         "S1: Gamma(5, 1)  |  equal bins  |  unimodal prior",  "y")
plot_one(fit2, true_d2,
         "S2: Lognormal    |  equal bins  |  OPEN last bin (30+, fit on [0, 80])",
         "y", xlim = c(0, 35))   # zoom past the empty open-bin tail
plot_one(fit3, true_d3,
         "S3: Bimodal       |  UNequal bins |  unconstrained",  "y")
plot_one(fit4, true_d4,
         "S4: Trimodal      |  equal bins  |  unconstrained",  "y")
mtext("Performance check on synthetic data: histogram, posterior mean (black), 90% credible band (blue), truth (red dashed)",
      side = 3, outer = TRUE, cex = 0.9, font = 2)
par(op)
dev.off()

# ---- Combined summary ----
all_summary <- rbind(sum1, sum2, sum3, sum4)
out_csv <- file.path(pkg_root, "tests", "synthetic_eval.csv")
write.csv(all_summary, out_csv, row.names = FALSE)

cat("\n=== Synthetic-data performance summary ===\n")
print(all_summary, row.names = FALSE)
cat("\nWrote:\n  ", out_png, "\n  ", out_csv, "\n", sep = "")
