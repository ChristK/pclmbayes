# =============================================================================
# run_all.R
# -----------------------------------------------------------------------------
# Top-level driver that reproduces every figure, table and JSON output for
# the paper in one call. Intended to be run from the package root:
#
#   Rscript inst/paper/run_all.R [out_dir]
#
# Argument:
#   out_dir -- folder to write all outputs into. Default: "paper_results".
#
# Outputs (all under out_dir):
#   figure1_simulation.png      -- four-panel simulation summary.
#   figure2_tbdeaths.png        -- TB Netherlands 1907 applied example.
#   figure3_ons_benchmark.png   -- 3x3 panel ONS benchmark.
#   figure4_open_tail.png       -- 2x2 panel open-ended-tail sensitivity.
#   sim_results.csv             -- per-replicate simulation metrics.
#   tbdeaths_table.csv          -- observed vs fitted counts (TB example).
#   ons_benchmark.csv           -- ONS benchmark metrics (3 datasets x 3
#                                  schemes x 3 methods).
#   ons_benchmark_summary.txt   -- median across schemes per dataset.
#   open_tail_results.csv       -- open-tail experiment metrics.
#   results.json                -- numeric summaries cited in the paper.
#   sessionInfo.txt, run.log,
#   ons.log, tail.log,
#   run.done, ons.done, tail.done
#
# All three numbered scripts are now data-self-contained: they read only
# from datasets bundled with the package (`tbdeaths1907` and four ONS
# extracts; see ?ons_deaths_ew_2022 etc.).
# =============================================================================

script_dir <- (function() {
  # 1) Rscript invocation
  ca <- commandArgs(trailingOnly = FALSE)
  m  <- regmatches(ca, regexpr("(?<=^--file=).*", ca, perl = TRUE))
  if (length(m)) return(normalizePath(dirname(m[[1L]])))
  # 2) sourced from R
  if (!is.null(sys.frames()[[1L]]$ofile)) {
    return(normalizePath(dirname(sys.frames()[[1L]]$ofile)))
  }
  # 3) fallback
  getwd()
})()

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "paper_results"

message(sprintf("[run_all] script dir : %s", script_dir))
message(sprintf("[run_all] out_dir    : %s", out_dir))

sys_call_args <- commandArgs(trailingOnly = FALSE)

# 01 - simulation + applied example: needs out_dir only.
message("\n--- 01_applied_and_simulation.R ---")
local({
  commandArgs <- function(trailingOnly = FALSE) {
    if (trailingOnly) c(out_dir) else sys_call_args
  }
  source(file.path(script_dir, "01_applied_and_simulation.R"),
         local = TRUE, echo = FALSE)
})

# 02 - ONS benchmark: needs out_dir only (data are bundled).
message("\n--- 02_ons_benchmark.R ---")
local({
  commandArgs <- function(trailingOnly = FALSE) {
    if (trailingOnly) c(out_dir) else sys_call_args
  }
  source(file.path(script_dir, "02_ons_benchmark.R"),
         local = TRUE, echo = FALSE)
})

# 03 - open tail: needs out_dir only (data are bundled).
message("\n--- 03_open_tail.R ---")
local({
  commandArgs <- function(trailingOnly = FALSE) {
    if (trailingOnly) c(out_dir) else sys_call_args
  }
  source(file.path(script_dir, "03_open_tail.R"),
         local = TRUE, echo = FALSE)
})

message("\n[run_all] All scripts complete. See: ", normalizePath(out_dir))
