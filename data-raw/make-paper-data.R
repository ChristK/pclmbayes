# =============================================================================
# data-raw/make-paper-data.R
# -----------------------------------------------------------------------------
# Run this script once (e.g. `Rscript data-raw/make-paper-data.R [ons_dir]`)
# to (re)generate the four ONS-derived datasets bundled with the package:
#
#   ons_deaths_ew_2022        England and Wales deaths 2022 by single year
#                              of age 0..104 plus 105+ open bin.
#                              Source: deaths_dr2022.xlsx (DR series).
#   ons_pop_uk_2024           UK mid-2024 population estimates by single
#                              year of age 0..89 plus 90+ open bin.
#                              Source: pop_mye24.xlsx, sheet MYE2-Persons.
#   ons_centenarians_uk_2024  UK 2024 population at ages 90..104 plus
#                              105+ open bin.
#                              Source: centenarians_uk_2024.csv.
#   ons_covid_due_ew          E&W deaths "due to COVID-19", 37 monthly
#                              cohorts March 2020 to March 2023, by
#                              single year of age 0..104 plus 105+ open
#                              bin.
#                              Source: covid_deaths_mar2023.xlsx (ONS
#                              "Single year of age and average age of
#                              death of people whose death was due to or
#                              involved coronavirus (COVID-19), England
#                              and Wales").
#
# The script is NOT installed with the package; end users do not need to
# run it. It exists so that:
#   (i)  the extraction logic from the original ONS sources is auditable;
#   (ii) the bundled datasets can be refreshed when ONS reissues a file
#        (e.g. the 2025 mid-year estimates supersede the 2024 ones).
#
# Dependencies: readxl (for the .xlsx files). Base R is otherwise enough.
#
# Usage
# -----
# From the package root:
#   Rscript data-raw/make-paper-data.R                  # uses ons_data/
#   Rscript data-raw/make-paper-data.R /path/to/ons     # explicit folder
# =============================================================================

args     <- commandArgs(trailingOnly = TRUE)
ons_dir  <- if (length(args) >= 1L) args[[1L]] else "ons_data"
out_dir  <- "data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required: install.packages('readxl').")
}

stopifnot(dir.exists(ons_dir))
message("[make-paper-data] ons_dir = ", normalizePath(ons_dir))

# Small helper: assemble a (age_lower, age_upper, age_label, count) frame.
# `ages_lo` and `ages_hi` parametrise the bin geometry; the last bin is the
# open one and `open_cap` provides its upper finite cap.
make_age_table <- function(counts, ages_lo, ages_hi, age_label) {
  data.frame(
    age_lower = as.integer(ages_lo),
    age_upper = as.integer(ages_hi),
    age_label = as.character(age_label),
    count     = as.integer(counts),
    stringsAsFactors = FALSE
  )
}


# ----- 1. ons_deaths_ew_2022 ------------------------------------------------
# Sheets 5 (males) and 6 (females) of dr2022.xlsx. The 2022 row's single-
# year columns are 0..104 plus a "105 and above" open bin (106 columns).
read_dr_sheet <- function(sheet) {
  raw <- readxl::read_excel(file.path(ons_dir, "deaths_dr2022.xlsx"),
                            sheet = sheet, skip = 4L)
  raw <- raw[raw[[1L]] == 2022, ]
  as.numeric(unlist(raw[1L, -(1L:2L)]))     # drop Year, All ages
}
males   <- read_dr_sheet("5")
females <- read_dr_sheet("6")
stopifnot(length(males) == length(females), length(males) == 106L)
counts  <- males + females
stopifnot(sum(counts) == 577160L)

ons_deaths_ew_2022 <- make_age_table(
  counts    = counts,
  ages_lo   = 0:105,
  ages_hi   = c(1:105, 110L),
  age_label = c(as.character(0:104), "105+")
)


# ----- 2. ons_pop_uk_2024 ---------------------------------------------------
# MYE2 - Persons sheet; header on row 8. The UK row (code K02000001) has,
# after the four metadata columns, 91 single-year-of-age counts: 0..89
# plus a "90+" open bin.
mye <- readxl::read_excel(file.path(ons_dir, "pop_mye24.xlsx"),
                          sheet = "MYE2 - Persons", skip = 7L)
uk  <- which(mye[[1L]] == "K02000001")
if (length(uk) != 1L) {
  uk <- which(toupper(trimws(mye[[2L]])) == "UNITED KINGDOM")
}
stopifnot(length(uk) == 1L)
pop_vals <- as.numeric(unlist(mye[uk, -(1L:4L)]))
stopifnot(length(pop_vals) == 91L, sum(pop_vals) == 69281437L)

ons_pop_uk_2024 <- make_age_table(
  counts    = pop_vals,
  ages_lo   = 0:90,
  ages_hi   = c(1:90, 100L),
  age_label = c(as.character(0:89), "90+")
)


# ----- 3. ons_centenarians_uk_2024 -----------------------------------------
# CSV layout (header on row 4):
#   Sex, Year, "90 and over", "90 to 99", "100 and over",
#   90, 91, ..., 99, 100, 101, 102, 103, 104, "105 and over"
# Single ages 90..104 are at columns 6..20; the open 105+ bin is at
# column 21.
cen <- utils::read.csv(file.path(ons_dir, "centenarians_uk_2024.csv"),
                       skip = 3L, header = FALSE, stringsAsFactors = FALSE,
                       check.names = FALSE)
# Strip thousand-separator commas where present, then coerce numeric
# columns (Year and counts) to numeric. The first column is text ("Sex").
cen[] <- lapply(cen, function(col) {
  if (is.character(col)) {
    out <- suppressWarnings(as.numeric(gsub(",", "", col)))
    if (all(is.na(out))) col else out
  } else col
})
row24 <- which(cen[[1L]] == "Persons" & cen[[2L]] == 2024)
stopifnot(length(row24) == 1L)
cen_counts <- c(as.numeric(cen[row24, 6L:20L]),
                as.numeric(cen[row24,      21L]))

ons_centenarians_uk_2024 <- make_age_table(
  counts    = cen_counts,
  ages_lo   = 90:105,
  ages_hi   = c(91:105, 110L),
  age_label = c(as.character(90:104), "105+")
)


# ----- 4. ons_covid_due_ew --------------------------------------------------
# Sheet 1; header on row 4. Each row is (Month, Cause, Aged under 1,
# Aged 1, ..., Aged 105 and over, Total). We sum across the 37 cohorts
# whose Cause is "Due to COVID-19", dropping the trailing "Total" column.
cov <- readxl::read_excel(file.path(ons_dir, "covid_deaths_mar2023.xlsx"),
                          sheet = "1", skip = 3L)
due <- cov[cov[["Cause"]] == "Due to COVID-19", ]
stopifnot(nrow(due) == 37L)
age_block <- due[, -c(1L:2L, ncol(due))]
covid_counts <- colSums(matrix(as.numeric(unlist(age_block)),
                               nrow = nrow(age_block), byrow = FALSE))
stopifnot(length(covid_counts) == 106L, sum(covid_counts) == 168418L)

ons_covid_due_ew <- make_age_table(
  counts    = covid_counts,
  ages_lo   = 0:105,
  ages_hi   = c(1:105, 110L),
  age_label = c(as.character(0:104), "105+")
)


# ----- Save .rda files ------------------------------------------------------
save(ons_deaths_ew_2022,
     file = file.path(out_dir, "ons_deaths_ew_2022.rda"), version = 2L)
save(ons_pop_uk_2024,
     file = file.path(out_dir, "ons_pop_uk_2024.rda"), version = 2L)
save(ons_centenarians_uk_2024,
     file = file.path(out_dir, "ons_centenarians_uk_2024.rda"), version = 2L)
save(ons_covid_due_ew,
     file = file.path(out_dir, "ons_covid_due_ew.rda"), version = 2L)

message("[make-paper-data] wrote 4 .rda files under '", out_dir, "'.")
