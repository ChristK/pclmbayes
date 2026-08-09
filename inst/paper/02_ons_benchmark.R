# =============================================================================
# 02_ons_benchmark.R
# -----------------------------------------------------------------------------
# Reproduces Figure 3 and ons_benchmark.csv for the paper.
#
# Real-data benchmark on three single-year-of-age series published by the
# UK Office for National Statistics (ONS) and shipped with this package:
#
#   * deaths     -- death registrations by single year of age, England and
#                   Wales, 2022 (pclmbayes::ons_deaths_ew_2022; total
#                   577,160).
#   * population -- mid-2024 population estimates by single year of age,
#                   UK persons (pclmbayes::ons_pop_uk_2024; total
#                   69,281,437).
#   * covid      -- deaths "due to COVID-19" by single year of age, summed
#                   across all 37 monthly cohorts March 2020 to March 2023,
#                   England and Wales (pclmbayes::ons_covid_due_ew;
#                   3.1-y total 168,418).
#
# For each series we:
#   1. Treat the published single-year-of-age series as the "truth" pi*.
#   2. Aggregate truth into wide bins under three grouping schemes:
#        - five_year  -- 0, 5, 10, ..., terminal
#        - ten_year   -- 0, 10, 20, ..., terminal
#        - realistic  -- typical published-summary mix: <1, 1-4, then
#                        5-year bands up to 85, with a single open
#                        terminal bin (population: 90+ since the source
#                        data stops at 90+; deaths/COVID: 85+).
#   3. Fit pclm / calibrated pclm / pclm_exact to the wide-bin counts and
#      compare the implied single-year-of-age counts back to the truth.
#
# Per-method metrics (matching ons_benchmark.csv):
#   * max_abs_count -- max |fitted_count_age - truth_count_age| (counts,
#                      not proportions).
#   * ise           -- integrated squared error between fitted and true
#                      single-year proportions.
#   * kl            -- KL divergence KL(truth || fit) of the single-year
#                      proportions, with a tiny epsilon for stability.
#   * max_bin_disc  -- max |fitted_wide_bin - observed_wide_bin| (the
#                      only metric that distinguishes pclm from calib /
#                      exact, by construction).
#
# Running this script
# -------------------
# From the package root:
#   Rscript inst/paper/02_ons_benchmark.R [out_dir]
# Default: out_dir = "paper_results". This script has no external data
# dependencies (the ONS extracts are bundled with the package).
#
# Author: pclmbayes paper code, reconstructed 2026-05-21.
# =============================================================================

suppressPackageStartupMessages({
  library(pclmbayes)
})


# ----- 0. Output directory --------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "paper_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(out_dir, "ons.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
on.exit({ sink(type = "output"); close(log_con) }, add = TRUE)


# ----- 1. Load the three bundled datasets -----------------------------------
# Each is a (age_lower, age_upper, age_label, count) data frame. The last
# row is an open terminal bin (e.g. 105+, capped at age_upper = 110).
data(ons_deaths_ew_2022,       package = "pclmbayes")
data(ons_pop_uk_2024,          package = "pclmbayes")
data(ons_covid_due_ew,         package = "pclmbayes")

# Convert each ONS dataset to a single-year-of-age vector on the integer
# fine grid [0, b]. Each row of the input is a band [age_lower, age_upper):
#   * width-1 rows are placed directly into their single-year cell;
#   * the wider open terminal bin's count is spread uniformly across the
#     ages it covers, so wide-bin aggregations of the fine grid recover
#     the original totals to machine precision and the fine grid itself
#     is equispaced (so bin breaks at integer ages align perfectly with
#     fine breaks, which is what calibrate() requires).
expand_to_singleyear <- function(df, label) {
  b_full <- max(df$age_upper)
  counts <- numeric(b_full)
  for (i in seq_len(nrow(df))) {
    lo  <- df$age_lower[i]
    hi  <- df$age_upper[i]
    w   <- hi - lo
    counts[(lo + 1L):hi] <- df$count[i] / w
  }
  list(counts  = counts,
       a       = 0,
       b       = b_full,
       label   = label,
       wb_in   = cbind(df$age_lower, df$age_upper),
       m_in    = df$count)
}
deaths_d <- expand_to_singleyear(ons_deaths_ew_2022, "deaths")
pop_d    <- expand_to_singleyear(ons_pop_uk_2024,    "population")
covid_d  <- expand_to_singleyear(ons_covid_due_ew,   "covid")

# Provenance: write the headline totals to the log so the reader can
# sanity-check the inputs.
for (d in list(deaths_d, pop_d, covid_d)) {
  cat(sprintf("%s: ages 0 to %d  total = %s\n",
              d$label, d$b,
              format(round(sum(d$counts)), big.mark = ",")))
}


# ----- 2. Grouping schemes --------------------------------------------------
# Build a wide-bin matrix [lo, hi] for one dataset's [a, b] under a named
# scheme. The "realistic" scheme is the same shape across datasets; the
# upper bound differs because the source data ends at different ages.
make_scheme <- function(scheme_name, a, b) {
  switch(scheme_name,
    five_year = {
      brks <- unique(c(seq(a, b, by = 5L), b))
      cbind(head(brks, -1L), tail(brks, -1L))
    },
    ten_year  = {
      brks <- unique(c(seq(a, b, by = 10L), b))
      cbind(head(brks, -1L), tail(brks, -1L))
    },
    realistic = {
      brks <- c(0, 1, 5, seq(10, 85, by = 5))
      brks <- brks[brks < b]
      brks <- c(brks, b)
      cbind(head(brks, -1L), tail(brks, -1L))
    },
    stop("Unknown scheme: ", scheme_name)
  )
}


# ----- 3. Aggregation, fit, and per-scheme metrics --------------------------
aggregate_to_wide <- function(counts_fine, fine_breaks, wb) {
  C <- bin_matrix(wb, fine_breaks)
  as.numeric(C %*% counts_fine)
}

fit_and_metrics <- function(dataset, scheme_name) {

  # Equispaced fine grid: one interval per single year of age over [0, b].
  # This guarantees that any wide-bin scheme with integer breakpoints
  # partitions the support cleanly -- a requirement of calibrate().
  counts_fine <- dataset$counts
  a <- dataset$a
  b <- dataset$b
  ngrid       <- b - a
  fine_breaks <- seq(a, b, by = 1L)

  wb     <- make_scheme(scheme_name, a, b)
  m_wide <- aggregate_to_wide(counts_fine, fine_breaks, wb)

  # Fit the three methods. `ndx` is left at its default, so the basis
  # dimension adapts to the scheme: K = min(max(J + 7, 20), ngrid, 200).
  # This matters here -- the five_year and realistic schemes give J = 19 to
  # 22 wide bins, and a fixed K = 20 would sit inside the weak-identification
  # band K < J + 4. See ?pclm, "Choosing the basis dimension".
  f_p <- pclm(m = m_wide, wide_breaks = wb,
              a = a, b = b, ngrid = ngrid,
              degree = 3L, penalty_order = 3L,
              select = "BIC")
  f_c <- calibrate(f_p)
  f_e <- pclm_exact(m = m_wide, wide_breaks = wb,
                    a = a, b = b, ngrid = ngrid,
                    degree = 3L, penalty_order = 3L)

  # Implied single-year counts = m_+ * pi_i.
  N        <- sum(counts_fine)
  fit_p_i  <- N * f_p$pi
  fit_c_i  <- N * f_c$pi
  fit_e_i  <- N * f_e$pi

  # Truth proportions
  p_true <- counts_fine / N

  # Stable KL: KL(p || q) = sum p log(p / q) with p log(p / q) = 0 when
  # p = 0; we floor q at a tiny epsilon to avoid log(0) when q ≈ 0.
  eps <- 1e-15
  kl  <- function(p, q) {
    pos <- p > 0
    sum(p[pos] * log(p[pos] / pmax(q[pos], eps)))
  }
  ise <- function(p, q) sum((p - q) ^ 2)
  max_abs <- function(fit_counts) max(abs(fit_counts - counts_fine))

  rows <- data.frame(
    dataset       = dataset$label,
    scheme        = scheme_name,
    method        = c("pclm", "calib", "exact"),
    max_abs_count = c(max_abs(fit_p_i), max_abs(fit_c_i), max_abs(fit_e_i)),
    ise           = c(ise(p_true, f_p$pi),
                      ise(p_true, f_c$pi),
                      ise(p_true, f_e$pi)),
    kl            = c(kl(p_true, f_p$pi),
                      kl(p_true, f_c$pi),
                      kl(p_true, f_e$pi)),
    max_bin_disc  = c(max(abs(f_p$fitted_counts - m_wide)),
                      max(abs(f_c$fitted_counts - m_wide)),
                      max(abs(f_e$fitted_counts - m_wide))),
    stringsAsFactors = FALSE
  )

  attr(rows, "densities") <- list(
    grid_mid = f_p$grid_mid,
    grid     = f_p$grid,
    truth    = p_true,
    pi_pclm  = f_p$pi,
    pi_calib = f_c$pi,
    pi_exact = f_e$pi,
    m_wide   = m_wide,
    wb       = wb
  )
  rows
}


# ----- 4. Run the benchmark -------------------------------------------------
datasets <- list(deaths = deaths_d, population = pop_d, covid = covid_d)
schemes  <- c("five_year", "ten_year", "realistic")

all_rows  <- list()
densities <- list()

for (dn in names(datasets)) {
  for (sn in schemes) {
    cat(sprintf("\n== %s / %s ==\n", dn, sn))
    rows <- fit_and_metrics(datasets[[dn]], sn)
    densities[[paste(dn, sn, sep = "/")]] <- attr(rows, "densities")
    attr(rows, "densities") <- NULL
    for (k in seq_len(nrow(rows))) {
      cat(sprintf("  %-5s: max_abs=%.1f  ISE=%.2e  KL=%.4f  bin_disc=%.2e\n",
                  rows$method[k], rows$max_abs_count[k],
                  rows$ise[k], rows$kl[k], rows$max_bin_disc[k]))
    }
    all_rows[[length(all_rows) + 1L]] <- rows
  }
}

bench_df <- do.call(rbind, all_rows)
write.csv(bench_df,
          file = file.path(out_dir, "ons_benchmark.csv"),
          row.names = FALSE)


# ----- 5. Figure 3: density panels (3 datasets x 3 schemes) ----------------
png(file.path(out_dir, "figure3_ons_benchmark.png"),
    width = 12, height = 9, units = "in", res = 200)
op <- par(mfrow = c(3L, 3L), mar = c(3.6, 4.0, 2.4, 1.0))
on.exit(par(op), add = TRUE)

row_ymax <- vapply(names(datasets), function(dn) {
  vals <- unlist(lapply(schemes, function(sn) {
    d <- densities[[paste(dn, sn, sep = "/")]]
    N <- sum(datasets[[dn]]$counts)
    Delta <- diff(d$grid)
    max(c(d$truth / Delta, d$pi_pclm / Delta, d$pi_exact / Delta)) * N
  }))
  max(vals) * 1.05
}, numeric(1L))
names(row_ymax) <- names(datasets)

for (dn in names(datasets)) {
  N <- sum(datasets[[dn]]$counts)
  for (sn in schemes) {
    d     <- densities[[paste(dn, sn, sep = "/")]]
    Delta <- diff(d$grid)
    mids  <- d$grid_mid

    plot(NA, xlim = range(d$grid),
         ylim = c(0, row_ymax[dn]),
         xlab = if (dn == names(datasets)[3L]) "Age (years)" else "",
         ylab = if (sn == schemes[1L])         "Count per year of age" else "",
         main = sprintf("%s / %s", dn, sn))
    points(mids, N * d$truth / Delta, pch = 16, cex = 0.5,
           col = adjustcolor("black", alpha.f = 0.6))
    lines(mids, N * d$pi_pclm  / Delta, col = "steelblue", lwd = 1.6, lty = 1)
    lines(mids, N * d$pi_calib / Delta, col = "darkorange", lwd = 1.6, lty = 2)
    lines(mids, N * d$pi_exact / Delta, col = "firebrick",  lwd = 1.6, lty = 3)
    if (dn == names(datasets)[1L] && sn == schemes[length(schemes)]) {
      legend("topleft", bty = "n", cex = 0.85,
             legend = c("truth (single year)", "pclm", "calib", "exact"),
             col    = c("black", "steelblue", "darkorange", "firebrick"),
             pch    = c(16, NA, NA, NA),
             lty    = c(NA, 1, 2, 3),
             lwd    = c(NA, 1.6, 1.6, 1.6))
    }
  }
}
dev.off()


# ----- 6. Per-dataset median across schemes --------------------------------
agg <- aggregate(cbind(max_abs_count, ise, kl, max_bin_disc)
                 ~ dataset + method,
                 data = bench_df,
                 FUN = function(x) median(x, na.rm = TRUE))
agg <- agg[order(agg$dataset, agg$method), ]
agg_print <- agg
agg_print$max_abs_count <- signif(agg_print$max_abs_count, 5L)
agg_print$ise           <- signif(agg_print$ise,           4L)
agg_print$kl            <- signif(agg_print$kl,            5L)
agg_print$max_bin_disc  <- signif(agg_print$max_bin_disc,  4L)
print(agg_print)
writeLines(c("ONS real-data benchmark: 3 datasets x 3 schemes x 3 methods.",
             "",
             "Median across schemes:",
             capture.output(print(agg_print))),
           con = file.path(out_dir, "ons_benchmark_summary.txt"))


# ----- Provenance ----------------------------------------------------------
file.create(file.path(out_dir, "ons.done"))
cat("\n=== ONS benchmark complete ===\n")
cat("Outputs: ons_benchmark.csv, figure3_ons_benchmark.png, ons_benchmark_summary.txt\n")
