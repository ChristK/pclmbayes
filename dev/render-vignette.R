# Source all R files into globalenv, then knit the vignette in place.
# Quick way to validate the vignette without a full R CMD build.
script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) sub("^--file=", "", m[1L]) else NA_character_
})()
pkg_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

# Pre-source the package code
for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                     full.names = TRUE)) source(f, local = globalenv())

# Make data available
load(file.path(pkg_root, "data", "bloodlead.rda"))
load(file.path(pkg_root, "data", "tbdeaths1907.rda"))

# Use knitr::knit to compile chunks, just to confirm they execute.
# We avoid pandoc -- a plain "knit" is enough to catch code errors.
out_md <- file.path(pkg_root, "tests", "vignette-output.md")
ok <- tryCatch({
  knitr::knit(input  = file.path(pkg_root, "vignettes", "pclmbayes-intro.Rmd"),
              output = out_md,
              quiet  = TRUE)
  TRUE
}, error = function(e) {
  cat("VIGNETTE ERROR: ", conditionMessage(e), "\n")
  FALSE
})
cat(sprintf("Vignette compiled: %s\n  output -> %s (%d bytes)\n",
            ok, out_md, if (file.exists(out_md)) file.size(out_md) else 0))
