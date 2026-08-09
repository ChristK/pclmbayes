# =============================================================================
# 01_applied_and_simulation.R
# -----------------------------------------------------------------------------
# Reproduces Figures 1-2 and the supporting CSVs for the paper:
#
#   * tbdeaths_table.csv  -- observed vs fitted counts for the Netherlands 1907
#                            tuberculosis-deaths applied example, under three
#                            methods (pclm, calibrated pclm, pclm_exact).
#   * figure2_tbdeaths.png -- corresponding density plot.
#   * sim_results.csv      -- long-format per-replicate accuracy/discrepancy
#                            metrics across {targets} x {grouping schemes}.
#   * figure1_simulation.png -- four-panel grid summarising the simulation.
#   * results.json         -- numeric summaries used in the manuscript text.
#
# Method labels used throughout:
#   * pclm  -- frequentist penalised composite link model (Lambert & Eilers
#              2009; pclmbayes::pclm), BIC-selected tau.
#   * calib -- post-hoc bin-conditional renormalisation of the pclm fit
#              (pclmbayes::calibrate) so that the implied wide-bin counts
#              equal the observed counts exactly.
#   * exact -- the constrained MAP fit (pclmbayes::pclm_exact) that
#              minimises (1/2) phi'P phi subject to gamma(phi) = m / m_+;
#              also exact, but with no discontinuities at bin boundaries.
#
# Running this script
# -------------------
# From the package root:
#   Rscript inst/paper/01_applied_and_simulation.R [out_dir]
# If `out_dir` is omitted it defaults to "paper_results/" in the current
# working directory. The script depends only on base R, the splines
# package (Recommended), and the installed pclmbayes package.
#
# Author: pclmbayes paper code, reconstructed 2026-05-21.
# =============================================================================

suppressPackageStartupMessages({
  library(pclmbayes)
  library(splines)
})


# ----- 0. Output directory ---------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
out_dir   <- if (length(args) >= 1L) args[[1L]] else "paper_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(out_dir, "run.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
on.exit({ sink(type = "output"); close(log_con) }, add = TRUE)

set.seed(20250507L)


# =============================================================================
# 1. Applied example: tuberculosis deaths, The Netherlands 1907
# =============================================================================
# Fits the three methods to the wide-bin counts and reports observed vs
# fitted counts. The pclm fit is the BIC-selected smoothing parameter on the
# default tau grid; calibrate() rescales the fit so that wide-bin totals
# match exactly; pclm_exact() solves the equality-constrained problem
# directly.

cat("== Applied example: TB Netherlands 1907 ==\n")
data(tbdeaths1907, package = "pclmbayes")
m_tb  <- tbdeaths1907$count
wb_tb <- cbind(tbdeaths1907$lower, tbdeaths1907$upper)

# Single year of age on (0, 120) -> 120 fine intervals of width 1.
a_tb       <- 0
b_tb       <- 120
ngrid_tb   <- 120
# `ndx` is left at its default, so the basis dimension is derived from the
# number of wide bins as K = min(max(J + 7, 20), ngrid, 200). With J = 12
# bands here that gives K = 20 cubic B-splines, the value used by Lambert
# and Eilers (2009). See ?pclm, "Choosing the basis dimension".
degree_tb  <- 3L
penord_tb  <- 3L       # 3rd-order penalty (paper recommendation)

# --- (a) frequentist PCLM, BIC-selected tau --------------------------------
t_pclm <- system.time(
  fit_pclm <- pclm(m = m_tb, wide_breaks = wb_tb,
                   a = a_tb, b = b_tb,
                   ngrid = ngrid_tb,
                   degree = degree_tb, penalty_order = penord_tb,
                   select = "BIC")
)["elapsed"]

# --- (b) calibrated PCLM ----------------------------------------------------
t_calib <- system.time(fit_calib <- calibrate(fit_pclm))["elapsed"]

# --- (c) pclm_exact ---------------------------------------------------------
t_exact <- system.time(
  fit_exact <- pclm_exact(m = m_tb, wide_breaks = wb_tb,
                          a = a_tb, b = b_tb,
                          ngrid = ngrid_tb,
                          degree = degree_tb, penalty_order = penord_tb)
)["elapsed"]

# Build the per-band table (rounded to 2 dp for the pclm fit; calib & exact
# are constructed to match m exactly so they are reported as integers).
tb_table <- data.frame(
  lower             = tbdeaths1907$lower,
  upper             = tbdeaths1907$upper,
  observed          = m_tb,
  fitted_pclm       = round(fit_pclm$fitted_counts,  2L),
  fitted_calibrated = round(fit_calib$fitted_counts, 2L),
  fitted_exact      = round(fit_exact$fitted_counts, 2L)
)
print(tb_table)
write.csv(tb_table,
          file = file.path(out_dir, "tbdeaths_table.csv"),
          row.names = FALSE)

# Per-method maximum absolute discrepancy between fitted and observed counts.
disc_pclm  <- max(abs(fit_pclm$fitted_counts  - m_tb))
disc_calib <- max(abs(fit_calib$fitted_counts - m_tb))
disc_exact <- max(abs(fit_exact$fitted_counts - m_tb))
cat(sprintf("max |fitted - observed|: pclm=%.4g  calib=%.3e  exact=%.3e\n",
            disc_pclm, disc_calib, disc_exact))


# ----- Figure 2: the three fitted densities on the same axes ---------------
# Densities are pi / Delta with Delta = (b - a) / ngrid.
delta_tb  <- (b_tb - a_tb) / ngrid_tb
mids_tb   <- fit_pclm$grid_mid
dens_pclm <- fit_pclm$pi  / delta_tb
dens_cal  <- fit_calib$pi / delta_tb
dens_ext  <- fit_exact$pi / delta_tb

# Histogram-style densities for the wide bins, so the reader can see the
# (irregular) age bands the smoothing has to bridge.
bin_w     <- wb_tb[, 2L] - wb_tb[, 1L]
bin_dens  <- (m_tb / sum(m_tb)) / bin_w

png(file.path(out_dir, "figure2_tbdeaths.png"),
    width = 8, height = 5, units = "in", res = 200)
op <- par(mar = c(4.2, 4.2, 2.6, 1.2))
on.exit(par(op), add = TRUE)
ymax <- max(c(dens_pclm, dens_cal, dens_ext, bin_dens)) * 1.08
plot(NA, xlim = c(0, 100), ylim = c(0, ymax),
     xlab = "Age (years)", ylab = "Density",
     main = "TB deaths, The Netherlands 1907")
# Wide-bin rectangles
for (j in seq_along(bin_dens)) {
  rect(xleft = wb_tb[j, 1L], ybottom = 0,
       xright = wb_tb[j, 2L], ytop = bin_dens[j],
       col = adjustcolor("grey80", alpha.f = 0.55),
       border = "grey50")
}
lines(mids_tb, dens_pclm, col = "black",      lwd = 2)
lines(mids_tb, dens_cal,  col = "steelblue",  lwd = 2, lty = 2)
lines(mids_tb, dens_ext,  col = "firebrick",  lwd = 2, lty = 3)
legend("topright",
       legend = c("pclm (smooth, approx)",
                  "calibrated",
                  "pclm_exact"),
       col = c("black", "steelblue", "firebrick"),
       lty = c(1, 2, 3), lwd = 2, bty = "n")
dev.off()


# =============================================================================
# 2. Simulation study
# =============================================================================
# Generates `nrep` synthetic datasets per cell of {target} x {scheme}, fits
# the three methods, and records:
#   * ise_*       -- integrated squared error between fitted density and
#                    truth (Riemann sum on the fine grid).
#   * max_disc_*  -- max |fitted_counts - observed_counts|.
#   * tv_*        -- total variation distance between fitted density and
#                    truth (0.5 * sum |f_hat - f| * Delta).
#   * time_*      -- wall-clock fit time (pclm and pclm_exact; calibrate
#                    sits on top of the pclm fit and is essentially free).
#
# Targets (truth densities on the common support [0, 100]):
#   lognormal -- meanlog = log(40), sdlog = 0.45 (typical adult mortality
#                shape).
#   bimodal   -- 0.6 N(20, 5^2) + 0.4 N(70, 8^2), truncated to [0, 100].
#   bathtub   -- mixture peaking at infancy and again at the oldest ages
#                (infant + senescent mortality).
#   uniform   -- a flat reference target on [10, 90].
#
# Grouping schemes (wide bins on [0, 100]):
#   five_year  -- 0, 5, 10, ..., 100  (J = 20)
#   ten_year   -- 0, 10, 20, ..., 100 (J = 10)
#   irregular  -- 0, 1, 5, 15, 25, 45, 65, 75, 85, 100  (J = 9; mixes
#                 narrow infant bins with wide adult bins, mimicking many
#                 historical mortality tabulations).

cat("\n== Simulation study ==\n")

# ----- Target density evaluators on the common support [0, 100] ------------
a_sim    <- 0
b_sim    <- 100
n_grid   <- 200L                              # fine evaluation grid
fine_brk <- seq(a_sim, b_sim, length.out = n_grid + 1L)
mids_sim <- (head(fine_brk, -1L) + tail(fine_brk, -1L)) / 2
delta    <- (b_sim - a_sim) / n_grid

truncate_to_density <- function(d) {
  # Helper: discretise a continuous density on the fine grid and renormalise
  # so it sums to 1 (i.e. so it is a discrete probability vector pi_i).
  p <- d * delta
  p / sum(p)
}

target_lognormal <- function() {
  truncate_to_density(dlnorm(mids_sim, meanlog = log(40), sdlog = 0.45))
}
target_bimodal <- function() {
  d <- 0.6 * dnorm(mids_sim, mean = 20, sd = 5) +
       0.4 * dnorm(mids_sim, mean = 70, sd = 8)
  truncate_to_density(d)
}
target_bathtub <- function() {
  # An infant-mortality peak (sharp, narrow) plus a senescent peak.
  d <- 0.20 * dnorm(mids_sim, mean =  1, sd =  2) +
       0.05 * dunif(mids_sim, min  =  5, max = 60) +
       0.75 * dnorm(mids_sim, mean = 75, sd = 10)
  truncate_to_density(d)
}
target_uniform <- function() {
  truncate_to_density(dunif(mids_sim, min = 10, max = 90))
}

truth_makers <- list(
  lognormal = target_lognormal,
  bimodal   = target_bimodal,
  bathtub   = target_bathtub,
  uniform   = target_uniform
)

# ----- Grouping schemes -----------------------------------------------------
schemes <- list(
  five_year = matrix(c(head(seq(0, 100, by =  5), -1L),
                       tail(seq(0, 100, by =  5), -1L)),
                     ncol = 2L),
  ten_year  = matrix(c(head(seq(0, 100, by = 10), -1L),
                       tail(seq(0, 100, by = 10), -1L)),
                     ncol = 2L),
  irregular = local({
    brks <- c(0, 1, 5, 15, 25, 45, 65, 75, 85, 100)
    matrix(c(head(brks, -1L), tail(brks, -1L)), ncol = 2L)
  })
)

# Pre-compute one bin matrix per scheme on the common fine grid.
C_per_scheme <- lapply(schemes, bin_matrix, fine_breaks = fine_brk)

# ----- Per-replicate metrics -----------------------------------------------
N_per   <- 1000L     # individuals per simulated dataset
nrep    <- 200L      # replicates per cell (matches results.json)
n_cells <- length(truth_makers) * length(schemes)

# Helper: convert observed wide-bin counts to a pclm/calib/exact fit and
# return the per-method metrics relative to a given truth.
fit_methods <- function(m_obs, wide_breaks, truth_pi, Delta) {

  ngrid_sim <- length(truth_pi)
  t_pclm <- system.time(
    f_p <- pclm(m = m_obs, wide_breaks = wide_breaks,
                a = a_sim, b = b_sim,
                ngrid = ngrid_sim,
                degree = 3L, penalty_order = 3L)
  )["elapsed"]
  f_c <- calibrate(f_p)
  t_ex <- system.time(
    f_e <- pclm_exact(m = m_obs, wide_breaks = wide_breaks,
                      a = a_sim, b = b_sim,
                      ngrid = ngrid_sim,
                      degree = 3L, penalty_order = 3L)
  )["elapsed"]

  # ISE between fitted density and truth (densities are pi / Delta).
  ise <- function(p_hat) sum(((p_hat - truth_pi) / Delta) ^ 2 * Delta)
  # TV  = 0.5 * integral |f_hat - f|
  tv  <- function(p_hat) 0.5 * sum(abs(p_hat - truth_pi))
  disc <- function(f) max(abs(f$fitted_counts - m_obs))

  list(
    ise_pclm       = ise(f_p$pi),
    ise_calib      = ise(f_c$pi),
    ise_exact      = ise(f_e$pi),
    max_disc_pclm  = disc(f_p),
    max_disc_calib = disc(f_c),
    max_disc_exact = disc(f_e),
    tv_pclm        = tv(f_p$pi),
    tv_calib       = tv(f_c$pi),
    tv_exact       = tv(f_e$pi),
    time_pclm      = unname(t_pclm),
    time_exact     = unname(t_ex),
    # Returned densities for later plotting (one example per cell).
    pi_pclm  = f_p$pi,
    pi_calib = f_c$pi,
    pi_exact = f_e$pi
  )
}

# Storage for the long-format CSV.
sim_rows <- vector("list", length(truth_makers) * length(schemes) * nrep)
row_i    <- 0L

# Storage for representative densities (one per cell, first replicate).
example_fits <- vector("list", n_cells)
example_idx  <- 0L

for (tg_name in names(truth_makers)) {
  truth_pi <- truth_makers[[tg_name]]()
  for (sc_name in names(schemes)) {
    wb_sc <- schemes[[sc_name]]
    # Wide-bin probabilities under the truth.
    gamma_truth <- as.numeric(C_per_scheme[[sc_name]] %*% truth_pi)
    if (any(gamma_truth < 0) || abs(sum(gamma_truth) - 1) > 1e-8) {
      stop(sprintf("Truth probabilities for %s / %s do not normalise.",
                   tg_name, sc_name))
    }
    example_idx <- example_idx + 1L
    for (r in seq_len(nrep)) {
      # Draw N_per individuals: multinomial counts in wide bins.
      m_obs <- as.integer(rmultinom(1L, size = N_per, prob = gamma_truth))
      fit   <- fit_methods(m_obs, wb_sc, truth_pi, delta)

      row_i <- row_i + 1L
      sim_rows[[row_i]] <- data.frame(
        target          = tg_name,
        scheme          = sc_name,
        rep             = r,
        ise_pclm        = fit$ise_pclm,
        ise_calib       = fit$ise_calib,
        ise_exact       = fit$ise_exact,
        max_disc_pclm   = fit$max_disc_pclm,
        max_disc_calib  = fit$max_disc_calib,
        max_disc_exact  = fit$max_disc_exact,
        tv_pclm         = fit$tv_pclm,
        tv_calib        = fit$tv_calib,
        tv_exact        = fit$tv_exact,
        time_pclm       = fit$time_pclm,
        time_exact      = fit$time_exact,
        stringsAsFactors = FALSE
      )
      if (r == 1L) {
        example_fits[[example_idx]] <- list(
          target = tg_name, scheme = sc_name,
          truth  = truth_pi,
          pi_pclm  = fit$pi_pclm,
          pi_calib = fit$pi_calib,
          pi_exact = fit$pi_exact
        )
      }
    }
    cat(sprintf("  done: %s / %s\n", tg_name, sc_name))
  }
}

sim_df <- do.call(rbind, sim_rows)
write.csv(sim_df,
          file = file.path(out_dir, "sim_results.csv"),
          row.names = FALSE)


# ----- Figure 1: four-panel summary of the simulation ----------------------
# Panel A: one representative density per target (truth + three fits) under
#          the irregular scheme (the hardest of the three).
# Panel B: boxplots of ISE by target, faceted by method.
# Panel C: boxplots of max bin discrepancy by target (log scale; the point
#          of the figure is that pclm has visible discrepancy while calib
#          and exact are at machine precision).
# Panel D: boxplots of total-variation distance by target, faceted by method.

png(file.path(out_dir, "figure1_simulation.png"),
    width = 10, height = 8, units = "in", res = 200)
op <- par(mfrow = c(2L, 2L), mar = c(4.0, 4.2, 2.2, 1.0))
on.exit(par(op), add = TRUE)

# ---- Panel A: representative densities (irregular scheme) ----------------
panelA_examples <- Filter(function(e) e$scheme == "irregular", example_fits)
plot(NA, xlim = c(0, 100),
     ylim = c(0, max(vapply(panelA_examples,
                            function(e) max(e$truth / delta), numeric(1))) * 1.1),
     xlab = "y", ylab = "Density",
     main = "A. Example fits under the irregular scheme")
cols   <- c("black", "steelblue", "firebrick", "darkgreen")
ltys   <- c(1, 2, 3, 4)
for (k in seq_along(panelA_examples)) {
  e <- panelA_examples[[k]]
  lines(mids_sim, e$truth / delta, col = cols[k], lwd = 2, lty = ltys[k])
}
legend("topright", legend = vapply(panelA_examples, `[[`, "", "target"),
       col = cols[seq_along(panelA_examples)],
       lty = ltys[seq_along(panelA_examples)],
       lwd = 2, bty = "n", cex = 0.9)

# ---- Panel B: ISE by target across all schemes ---------------------------
ise_long <- with(sim_df, data.frame(
  target = rep(target, 3L),
  method = rep(c("pclm", "calib", "exact"), each = nrow(sim_df)),
  ise    = c(ise_pclm, ise_calib, ise_exact)
))
ise_long$cell <- with(ise_long, paste(target, method, sep = "\n"))
ise_long$cell <- factor(ise_long$cell,
                        levels = as.vector(outer(
                          c("pclm", "calib", "exact"),
                          c("lognormal", "bimodal", "bathtub", "uniform"),
                          function(m, t) paste(t, m, sep = "\n"))))
boxplot(ise ~ cell, data = ise_long, log = "y",
        las = 2, cex.axis = 0.7,
        xlab = "", ylab = "Integrated squared error (log scale)",
        main = "B. ISE by target x method")

# ---- Panel C: max bin discrepancy (log scale, +epsilon for plot stability)
eps <- 1e-14
disc_long <- with(sim_df, data.frame(
  target = rep(target, 3L),
  method = rep(c("pclm", "calib", "exact"), each = nrow(sim_df)),
  disc   = c(max_disc_pclm, max_disc_calib, max_disc_exact) + eps
))
disc_long$cell <- with(disc_long, paste(target, method, sep = "\n"))
disc_long$cell <- factor(disc_long$cell,
                         levels = levels(ise_long$cell))
boxplot(disc ~ cell, data = disc_long, log = "y",
        las = 2, cex.axis = 0.7,
        xlab = "", ylab = expression(max[j] ~ "|" ~ hat(m)[j] - m[j] ~ "|"
                                     ~ "(log scale)"),
        main = "C. Bin-total discrepancy")

# ---- Panel D: TV distance ------------------------------------------------
tv_long <- with(sim_df, data.frame(
  target = rep(target, 3L),
  method = rep(c("pclm", "calib", "exact"), each = nrow(sim_df)),
  tv     = c(tv_pclm, tv_calib, tv_exact)
))
tv_long$cell <- with(tv_long, paste(target, method, sep = "\n"))
tv_long$cell <- factor(tv_long$cell, levels = levels(ise_long$cell))
boxplot(tv ~ cell, data = tv_long,
        las = 2, cex.axis = 0.7,
        xlab = "", ylab = "Total variation distance",
        main = "D. TV(fit, truth) by target x method")

dev.off()


# =============================================================================
# 3. Numeric summaries -> results.json
# =============================================================================
# These are the numbers the manuscript text refers to. We reproduce the
# exact field names already used in paper_results/results.json so that
# downstream Word/LaTeX templates do not need to change.

# Worst pclm band on the TB example
disc_tb_byband <- abs(fit_pclm$fitted_counts - m_tb)
worst <- which.max(disc_tb_byband)

# Medians across the simulation
med <- function(x) median(x)
q25 <- function(x) unname(quantile(x, 0.25))
q75 <- function(x) unname(quantile(x, 0.75))

results <- list(
  tb = list(
    total_deaths              = sum(m_tb),
    n_bands                   = nrow(tbdeaths1907),
    max_disc_pclm             = unname(disc_pclm),
    max_disc_calib            = unname(disc_calib),
    max_disc_exact            = unname(disc_exact),
    worst_band_pclm_lower     = tbdeaths1907$lower[worst],
    worst_band_pclm_upper     = tbdeaths1907$upper[worst],
    worst_band_pclm_disc      = unname(max(disc_tb_byband)),
    selected_tau_pclm         = unname(fit_pclm$tau),
    edf_pclm                  = unname(fit_pclm$edf),
    iter_exact                = fit_exact$iter,
    constraint_residual       = fit_exact$constraint_residual,
    time_pclm_s               = unname(t_pclm),
    time_calib_s              = unname(t_calib),
    time_exact_s              = unname(t_exact)
  ),
  sim = list(
    nrep_per_cell             = nrep,
    N_per_dataset             = N_per,
    n_targets                 = length(truth_makers),
    n_schemes                 = length(schemes),
    median_ise_pclm           = med(sim_df$ise_pclm),
    median_ise_calib          = med(sim_df$ise_calib),
    median_ise_exact          = med(sim_df$ise_exact),
    median_max_disc_pclm      = med(sim_df$max_disc_pclm),
    q25_max_disc_pclm         = q25(sim_df$max_disc_pclm),
    q75_max_disc_pclm         = q75(sim_df$max_disc_pclm),
    median_max_disc_calib     = med(sim_df$max_disc_calib),
    median_max_disc_exact     = med(sim_df$max_disc_exact),
    median_tv_pclm            = med(sim_df$tv_pclm),
    median_tv_calib           = med(sim_df$tv_calib),
    median_tv_exact           = med(sim_df$tv_exact),
    tv_calib_over_pclm        = med(sim_df$tv_calib) / med(sim_df$tv_pclm),
    median_time_pclm_ms       = 1000 * med(sim_df$time_pclm),
    median_time_exact_ms      = 1000 * med(sim_df$time_exact),
    pct_ise_calib_vs_pclm     = 100 *
      (med(sim_df$ise_calib) - med(sim_df$ise_pclm)) / med(sim_df$ise_pclm),
    pct_ise_exact_vs_pclm     = 100 *
      (med(sim_df$ise_exact) - med(sim_df$ise_pclm)) / med(sim_df$ise_pclm)
  )
)

# Minimal JSON writer (avoids depending on jsonlite). Handles only nested
# named lists of scalars (which is all we need here).
to_json_scalar <- function(x) {
  if (is.character(x)) shQuote(x, type = "cmd")
  else if (is.logical(x)) tolower(x)
  else if (is.numeric(x)) formatC(x, digits = 6L, format = "g")
  else stop("Unsupported scalar type")
}
to_json <- function(x, indent = "") {
  if (is.list(x)) {
    nm <- names(x)
    body <- vapply(seq_along(x), function(i) {
      sprintf('%s  "%s": %s', indent, nm[i],
              to_json(x[[i]], indent = paste0(indent, "  ")))
    }, character(1L))
    paste0("{\n", paste(body, collapse = ",\n"), "\n", indent, "}")
  } else {
    to_json_scalar(x)
  }
}
writeLines(to_json(results), con = file.path(out_dir, "results.json"))


# ----- Provenance ----------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo.txt"))
file.create(file.path(out_dir, "run.done"))

cat("\n=== DONE ===\n")
cat("Outputs written to:", normalizePath(out_dir), "\n")
cat("Files: figure1_simulation.png, figure2_tbdeaths.png, results.json, ",
    "run.log, sessionInfo.txt, sim_results.csv, tbdeaths_table.csv\n",
    sep = "")
