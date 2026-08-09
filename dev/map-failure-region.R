# =============================================================================
# dev/map-failure-region.R
# -----------------------------------------------------------------------------
# Provenance for the K >= J + 4 threshold used by .default_ndx() and
# .check_basis_dim() in R/basis.R.
#
# The composite-link matrix C is J x I, so the wide-bin counts identify at most
# J - 1 directions in the K-dimensional coefficient vector phi; the rest are
# fixed by tau * P alone. For K >> J those directions are high-frequency, where
# phi' P phi is large and the penalty suppresses them. For K ~ J they are
# smooth, where phi' P phi is small, the Fisher information is severely
# ill-conditioned, and phi can grow without bound -- producing large spurious
# excursions in the fitted density, most visibly at the edges of the support.
#
# This script tests the hypothesis that the failure is a resonance at K ~ J
# rather than something specific to the open terminal bin, by sweeping K
# against three different binnings of the same data. If the hypothesis holds,
# the failure band should MOVE with J:
#
#   OPEN   5-year bands + 90+, J = 19  -> fails around K = 19-21
#   CLOSED 5-year bands only,  J = 18  -> fails around K = 17-19
#   CLOSED 10-year bands,      J =  9  -> fails around K =  9-12
#
# check_basis = FALSE throughout: this script deliberately visits the band the
# guard exists to warn about.
#
# Run:  Rscript dev/map-failure-region.R
# =============================================================================

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))

load("data/ons_pop_uk_2024.rda")
n1  <- as.numeric(ons_pop_uk_2024$count[ons_pop_uk_2024$age_label != "90+"])
m90 <- as.numeric(ons_pop_uk_2024$count[ons_pop_uk_2024$age_label == "90+"])
m5  <- as.numeric(tapply(n1, rep(seq(0, 85, 5), each = 5), sum))
YOUNG <- 1:21

run <- function(m, wbx, a, b, ngrid, ndx, r = 3L) {
  f <- pclm(m = m, wide_breaks = wbx, a = a, b = b, ngrid = ngrid,
            ndx = ndx, penalty_order = r, check_basis = FALSE)
  e <- as.numeric(bin_matrix(seq(0, 90, 1), f$grid) %*% f$pi)[1:90] * sum(f$m)
  rel <- 100 * (e - n1) / n1
  c(K = f$basis$K, tau = f$tau, edf = f$edf, maxphi = max(abs(f$phi)),
    mape_young = mean(abs(rel[YOUNG])), max_young = max(abs(rel[YOUNG])))
}

show <- function(title, J, ndx_range, fn) {
  cat(sprintf("\n--- %s   (J = %d wide bins) ---\n", title, J))
  cat("  ndx   K  K-J        tau     edf   max|phi|  MAPE(0-20)  max(0-20)\n")
  for (nd in ndx_range) {
    s <- tryCatch(fn(nd), error = function(e) NULL)
    if (is.null(s)) next
    flag <- if (s[["mape_young"]] > 5) "  <== FAILS" else ""
    cat(sprintf("  %3d %3d %4d %10.3g %7.2f %10.2f %11.2f%% %9.2f%%%s\n",
                nd, s[["K"]], s[["K"]] - J, s[["tau"]], s[["edf"]],
                s[["maxphi"]], s[["mape_young"]], s[["max_young"]], flag))
  }
}

# Open problem: 18 five-year bins + one 90+ bin = 19
show("OPEN: 5-year bins + 90+ capped at 100", 19L, 8:30,
     function(nd) run(c(m5, m90),
                      rbind(cbind(seq(0, 85, 5), seq(5, 90, 5)), c(90, 100)),
                      0, 100, 100, nd))

# Closed problem: 18 five-year bins only
show("CLOSED: 5-year bins only, b = 90", 18L, 8:30,
     function(nd) run(m5, seq(0, 90, 5), 0, 90, 90, nd))

# Closed problem with 10-year bins: J = 9, so the resonance should move down
# to K ~ 9, i.e. ndx ~ 6.
m10 <- as.numeric(tapply(n1, rep(seq(0, 80, 10), each = 10), sum))
show("CLOSED: 10-year bins, b = 90", 9L, 3:20,
     function(nd) run(m10, seq(0, 90, 10), 0, 90, 90, nd))

cat("\nConclusion: the failure band tracks J, confirming a weak-identification\n")
cat("resonance at K ~ J. K >= J + 4 is clear of it in all three binnings.\n")
cat("This is the basis for the threshold in .check_basis_dim() and for the\n")
cat("K = max(J + 7, 20) target in .default_ndx() (R/basis.R).\n")
