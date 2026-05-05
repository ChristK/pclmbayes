# =============================================================================
# data-raw/make-data.R
# -----------------------------------------------------------------------------
# Run this script ONCE (e.g. `Rscript data-raw/make-data.R`) to (re)generate
# the package datasets in `data/`. The script is not installed with the
# package; users do not need to run it.
#
# Datasets:
#   - bloodlead     : real data from Hasselblad, Stead and Galke (1980), as
#                     reproduced exactly in Lambert and Eilers (2009),
#                     Section 6.1.
#   - tbdeaths1907  : an ILLUSTRATIVE reconstruction of the wide-age-band
#                     tuberculosis-deaths dataset used in Lambert and Eilers
#                     (2009), Section 6.2. The original underlying data come
#                     from Statistics Netherlands (CBS) historical archives
#                     and were "unpublished" at the time of the paper. The
#                     reconstruction below preserves: (i) the published total
#                     of 9440 deaths; (ii) the wide age-band structure
#                     suggested by the paper's Fig. 6; (iii) the typical
#                     shape of pre-antibiotic-era TB mortality (peak in young
#                     adults aged 20-40, smaller peak in infancy). It is
#                     intended for code testing and demonstration, NOT for
#                     historical inference. To use real CBS data, replace
#                     this object before installing.
#
# References:
#   Hasselblad, V., Stead, A. G. and Galke, W. (1980). Analysis of coarsely
#     grouped data from the lognormal distribution. JASA 75, 771-778.
#   Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
#     from grouped continuous data. CSDA 53, 1388-1399.
# =============================================================================

# ---- bloodlead --------------------------------------------------------------
# Source: Lambert & Eilers (2009), Section 6.1, p. 1394 (table verbatim).
# The original "65+" upper bin was bounded above by 80 in the paper's analysis.
bloodlead <- data.frame(
  lower = c( 0, 15, 25, 35, 45, 55, 65),
  upper = c(15, 25, 35, 45, 55, 65, 80),
  count = c(27, 71, 32,  6,  3,  0,  0)
)
stopifnot(sum(bloodlead$count) == 139L)

# ---- tbdeaths1907 -----------------------------------------------------------
# ILLUSTRATIVE reconstruction (see header comment for full disclosure).
# Twelve wide bins on (0, 100) plus an extra zero-count bin (100, 120) to
# encourage smooth tapering at the upper tail (Lambert & Eilers 2009,
# Section 4). Counts sum exactly to 9440.
tbdeaths1907 <- data.frame(
  lower = c(  0,   1,   5,  10,  15,  20,  30,  40,  50,  60,  70,  80, 100),
  upper = c(  1,   5,  10,  15,  20,  30,  40,  50,  60,  70,  80, 100, 120),
  count = c(540, 980, 320, 240, 580, 1880, 1690, 1230, 850, 590, 380, 160,  0)
)
stopifnot(sum(tbdeaths1907$count) == 9440L)

# ---- Save -------------------------------------------------------------------
# Use save() with version 2 for compatibility with R >= 3.5.
out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save(bloodlead,    file = file.path(out_dir, "bloodlead.rda"),    version = 2)
save(tbdeaths1907, file = file.path(out_dir, "tbdeaths1907.rda"), version = 2)
message("Wrote: ", file.path(out_dir, "bloodlead.rda"))
message("Wrote: ", file.path(out_dir, "tbdeaths1907.rda"))
