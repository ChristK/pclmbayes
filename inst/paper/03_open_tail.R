# =============================================================================
# 03_open_tail.R
# -----------------------------------------------------------------------------
# Reproduces Figure 4 and open_tail_results.csv for the paper.
#
# Question: how well do the three methods (pclm, calibrated pclm,
# pclm_exact) recover the right-tail of the age distribution when the
# published data closes the top with an open terminal bin?
#
# Truth: UK mid-2024 population by single year of age (ages 0..104) plus
# a "105+" open bin. We reconstruct the single-year-of-age truth by
# splicing two ONS sources, both bundled with the package:
#   * pclmbayes::ons_pop_uk_2024          single-year ages 0..89 plus
#                                          a "90+" aggregated bin.
#   * pclmbayes::ons_centenarians_uk_2024 single-year ages 90..104 plus
#                                          a "105+" open bin.
# The 90+ MYE2 total and the 90..104 + 105+ centenarian total may differ
# by a few units because the two ONS releases are rounded independently;
# we keep the centenarian breakdown as the truth at ages >= 90.
#
# Experiment: for each candidate terminal age X in {80, 85, 90, 95}, build
# a 9-bin grouping scheme
#       0-<1, 1-4, 5-14, 15-24, 25-44, 45-64, 65-74, 75-(X-1), X+
# and aggregate the truth into wide-bin counts under that scheme. Fit
# the three methods, and measure how well they recover the single-year
# pattern across ages 80..104 (where the open bin's width changes from
# 25 down to 10 years).
#
# Per-method metrics (matching open_tail_results.csv):
#   * open_observed -- observed wide-bin count in the terminal bin.
#   * open_fitted   -- fitted wide-bin count in the terminal bin.
#   * open_disc     -- |open_fitted - open_observed|.
#   * max_abs_80_104 -- max |fitted_count_a - truth_count_a| over
#                       a = 80..104.
#   * sum_abs_80_104 -- sum |fitted_count_a - truth_count_a|.
#   * ise_tail      -- integrated squared error of the single-year
#                      proportion estimates over ages 80..104.
#
# Running this script
# -------------------
# From the package root:
#   Rscript inst/paper/03_open_tail.R [out_dir]
# Default: out_dir = "paper_results". This script has no external data
# dependencies (the ONS extracts are bundled with the package).
#
# Author: pclmbayes paper code, reconstructed 2026-05-21.
# =============================================================================

suppressPackageStartupMessages({
  library(pclmbayes)
})


# ----- 0. I/O ---------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "paper_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(out_dir, "tail.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
on.exit({ sink(type = "output"); close(log_con) }, add = TRUE)


# ----- 1. Build the single-year-of-age truth -------------------------------
# 1a. Ages 0..89 from MYE2 - Persons; we deliberately discard the MYE2
#     "90+" aggregate because the centenarian release gives finer detail.
# 1b. Ages 90..104 + "105+" open bin from the centenarian release.
data(ons_pop_uk_2024,          package = "pclmbayes")
data(ons_centenarians_uk_2024, package = "pclmbayes")

# 1a. ages 0..89 from the main mid-year estimate.
mye_singles  <- ons_pop_uk_2024[ons_pop_uk_2024$age_label != "90+", ]
stopifnot(nrow(mye_singles) == 90L, all(mye_singles$age_lower == 0:89))
pop_0_89     <- mye_singles$count

# 1b. ages 90..104 + 105+ from the centenarian release.
cen_singles  <- ons_centenarians_uk_2024[
                  ons_centenarians_uk_2024$age_label != "105+", ]
cen_open     <- ons_centenarians_uk_2024[
                  ons_centenarians_uk_2024$age_label == "105+", "count"]
stopifnot(nrow(cen_singles) == 15L,
          all(cen_singles$age_lower == 90:104))
ages_90_104  <- cen_singles$count
open_105     <- as.integer(cen_open)

# Truth: ages 0..104 + 105+ open bin (length 106).
truth_counts <- c(pop_0_89, ages_90_104, open_105)
cat("Truth: ages 0..104 plus 105+ open bucket.\n")
cat(sprintf("  Total ages 0..104: %s\n",
            format(sum(truth_counts[1L:105L]), big.mark = ",")))
cat(sprintf("  105+ open count  : %d\n", open_105))
cat(sprintf("  Grand total (ages 0..105+): %s\n",
            format(sum(truth_counts), big.mark = ",")))
cat(sprintf("  Ages 80..104 sub-total: %s\n",
            format(sum(truth_counts[81L:105L]), big.mark = ",")))


# ----- 2. Set up the fine grid and the wide-bin schemes --------------------
# Fine support [0, 105): one fine interval per integer year of age.
# `truth_fine` holds the single-year truth for ages 0..104 (length 105);
# the 105+ open count is kept separately and added to the open terminal
# wide bin in the observation step (see below), so the methods see the
# real ONS observation -- an open bin whose count includes 105+ -- while
# the per-single-year metric over ages 80..104 is computed against the
# uncontaminated single-year truth.
a <- 0
b <- 105
ngrid       <- 105L
fine_breaks <- seq(a, b, by = 1L)        # integer breaks 0, 1, ..., 105
truth_fine  <- truth_counts[1L:105L]      # ages 0..104

make_scheme_open <- function(X) {
  brks <- c(0, 1, 5, 15, 25, 45, 65, 75, X, 105)
  cbind(head(brks, -1L), tail(brks, -1L))
}

aggregate_to_wide <- function(counts_fine, fine_breaks, wb) {
  C <- bin_matrix(wb, fine_breaks)
  as.numeric(C %*% counts_fine)
}


# ----- 3. Run the experiment ------------------------------------------------
terminal_grid <- c(80L, 85L, 90L, 95L)
all_rows      <- list()
density_store <- list()

for (X in terminal_grid) {
  wb    <- make_scheme_open(X)
  # Wide-bin counts: aggregate the single-year truth, then ADD the 105+
  # open count to the last (terminal) wide bin -- exactly what an analyst
  # working from the published grouped data would observe.
  m_obs <- aggregate_to_wide(truth_fine, fine_breaks, wb)
  m_obs[nrow(wb)] <- m_obs[nrow(wb)] + open_105

  cat(sprintf("\n=== Terminal at %d+: %d bins; m =\n", X, nrow(wb)))
  print(data.frame(lower = wb[, 1L], upper = wb[, 2L], count = m_obs))
  cat("  band widths:", paste(wb[, 2L] - wb[, 1L], collapse = " "), "\n")
  cat("  open-bin width (terminal..104+open):",
      wb[nrow(wb), 2L] - wb[nrow(wb), 1L], "\n")

  # The basis dimension is derived from the number of wide bins (J = 9 in
  # every terminal scheme here). See ?pclm, "Choosing the basis dimension".
  f_p <- pclm(m = m_obs, wide_breaks = wb,
              a = a, b = b, ngrid = ngrid,
              degree = 3L, penalty_order = 3L,
              select = "BIC")
  f_c <- calibrate(f_p)
  f_e <- pclm_exact(m = m_obs, wide_breaks = wb,
                    a = a, b = b, ngrid = ngrid,
                    degree = 3L, penalty_order = 3L)

  N        <- sum(m_obs)         # = sum(truth_fine) + open_105
  fit_p_i  <- N * f_p$pi
  fit_c_i  <- N * f_c$pi
  fit_e_i  <- N * f_e$pi

  # Look at the open bin (last wide bin) and the disaggregated fit over
  # ages 80..104 (fine indices 81..105 in 1-based R indexing).
  idx_tail <- 81L:105L

  per_method <- function(fit_i) {
    open_fitted <- sum(fit_i[(X + 1L):105L])    # fitted total in last wide bin
    list(
      open_observed = m_obs[nrow(wb)],
      open_fitted   = open_fitted,
      open_disc     = abs(open_fitted - m_obs[nrow(wb)]),
      max_abs_80_104 = max(abs(fit_i[idx_tail] - truth_fine[idx_tail])),
      sum_abs_80_104 = sum(abs(fit_i[idx_tail] - truth_fine[idx_tail])),
      ise_tail = sum(((fit_i[idx_tail] - truth_fine[idx_tail]) / N) ^ 2)
    )
  }

  res_p <- per_method(fit_p_i)
  res_c <- per_method(fit_c_i)
  res_e <- per_method(fit_e_i)

  rows <- data.frame(
    terminal       = X,
    method         = c("pclm", "calib", "exact"),
    open_observed  = c(res_p$open_observed, res_c$open_observed, res_e$open_observed),
    open_fitted    = c(res_p$open_fitted,   res_c$open_fitted,   res_e$open_fitted),
    open_disc      = c(res_p$open_disc,     res_c$open_disc,     res_e$open_disc),
    max_abs_80_104 = c(res_p$max_abs_80_104, res_c$max_abs_80_104, res_e$max_abs_80_104),
    sum_abs_80_104 = c(res_p$sum_abs_80_104, res_c$sum_abs_80_104, res_e$sum_abs_80_104),
    ise_tail       = c(res_p$ise_tail,      res_c$ise_tail,      res_e$ise_tail),
    stringsAsFactors = FALSE
  )
  all_rows[[length(all_rows) + 1L]] <- rows

  density_store[[as.character(X)]] <- list(
    grid_mid = f_p$grid_mid,
    pi_pclm  = fit_p_i,
    pi_calib = fit_c_i,
    pi_exact = fit_e_i,
    truth    = truth_fine
  )
}

tail_df <- do.call(rbind, all_rows)
print(tail_df)
write.csv(tail_df,
          file = file.path(out_dir, "open_tail_results.csv"),
          row.names = FALSE)


# ----- 4. Figure 4: 2x2 panel grid -----------------------------------------
png(file.path(out_dir, "figure4_open_tail.png"),
    width = 10, height = 8, units = "in", res = 200)
op <- par(mfrow = c(2L, 2L), mar = c(4.0, 4.4, 2.4, 1.0))
on.exit(par(op), add = TRUE)

x_show <- 70L:104L
y_max  <- max(unlist(lapply(density_store, function(d) {
  max(c(d$pi_pclm[x_show + 1L],
        d$pi_calib[x_show + 1L],
        d$pi_exact[x_show + 1L],
        d$truth[x_show + 1L]))
}))) * 1.05

for (X in terminal_grid) {
  d <- density_store[[as.character(X)]]
  plot(x_show, d$truth[x_show + 1L], pch = 16, cex = 0.7,
       col = adjustcolor("black", alpha.f = 0.7),
       ylim = c(0, y_max), xlab = "Age (years)",
       ylab = "Count per year of age",
       main = sprintf("Terminal at %d+", X))
  lines(x_show, d$pi_pclm[x_show + 1L],  col = "steelblue",  lwd = 1.7, lty = 1)
  lines(x_show, d$pi_calib[x_show + 1L], col = "darkorange", lwd = 1.7, lty = 2)
  lines(x_show, d$pi_exact[x_show + 1L], col = "firebrick",  lwd = 1.7, lty = 3)
  abline(v = X, col = "grey50", lty = 3)
  if (X == terminal_grid[1L]) {
    legend("topright", bty = "n", cex = 0.8,
           legend = c("truth (single year)", "pclm", "calib", "exact"),
           col    = c("black", "steelblue", "darkorange", "firebrick"),
           pch    = c(16, NA, NA, NA),
           lty    = c(NA, 1, 2, 3),
           lwd    = c(NA, 1.7, 1.7, 1.7))
  }
}
dev.off()


# ----- Provenance ----------------------------------------------------------
file.create(file.path(out_dir, "tail.done"))
cat("\n=== Open-ended-tail analysis complete ===\n")
cat("Outputs: open_tail_results.csv, figure4_open_tail.png\n")
