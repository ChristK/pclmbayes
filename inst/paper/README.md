# Paper reproduction scripts

This directory contains the R scripts that reproduce every figure, table
and JSON summary referenced in the manuscript

> *pclmbayes: exact-preservation extensions of the penalised composite
> link model for density estimation from grouped continuous data*
> (Christodoulou et al., draft AJE 2026).

The scripts were reconstructed on 2026-05-21 from the output artefacts
that were originally written into `paper_results/` (the run logs,
`results.json`, the CSVs, and the four PNG figures). They reproduce the
intent of the original runs; exact numerical values for the
random-seed-dependent steps (the simulation study) will match the
originals when both the seed (20250507) and the package version are
the same as those used originally.

## Layout

```
inst/paper/
  README.md                       (this file)
  run_all.R                       end-to-end driver
  01_applied_and_simulation.R     -> figure1, figure2, sim_results.csv,
                                     tbdeaths_table.csv, results.json
  02_ons_benchmark.R              -> figure3, ons_benchmark.csv,
                                     ons_benchmark_summary.txt
  03_open_tail.R                  -> figure4, open_tail_results.csv
```

## Dependencies

- R >= 4.1
- The `pclmbayes` package, installed locally
  (e.g. `devtools::install_local(".")` from the package root).
- `splines` (Recommended; ships with R).
- **No CRAN dependencies for the scripts themselves.** The four ONS
  datasets the benchmark needs are bundled with the package (see
  `?ons_deaths_ew_2022`, `?ons_pop_uk_2024`,
  `?ons_centenarians_uk_2024`, `?ons_covid_due_ew`). Package `readxl`
  is in `Suggests` only because `data-raw/make-paper-data.R` uses it
  to re-extract from the original ONS spreadsheets when refreshing the
  bundled datasets.

## Reproducing the outputs

From the package root (where `DESCRIPTION` lives):

```sh
Rscript inst/paper/run_all.R                       # writes to ./paper_results
```

Or run the three scripts individually:

```sh
Rscript inst/paper/01_applied_and_simulation.R  paper_results
Rscript inst/paper/02_ons_benchmark.R           paper_results
Rscript inst/paper/03_open_tail.R               paper_results
```

Outputs are written to `paper_results/` (created if absent).

## Bundled datasets

| Dataset                     | Rows | Notes                                              |
|-----------------------------|------|----------------------------------------------------|
| `tbdeaths1907`              | 13   | Illustrative reconstruction (see `?tbdeaths1907`). |
| `ons_deaths_ew_2022`        | 106  | Deaths registered E&W, 2022 (total 577,160).        |
| `ons_pop_uk_2024`           | 91   | UK mid-2024 population (total 69,281,437).          |
| `ons_centenarians_uk_2024`  | 16   | UK 2024, ages 90+ broken down.                      |
| `ons_covid_due_ew`          | 106  | E&W "due to COVID-19" deaths, Mar 2020 - Mar 2023. |

All four ONS extracts share a common schema:
`(age_lower, age_upper, age_label, count)`. The last row of each is an
open terminal bin with `age_label` ending in `+`.

## Refreshing the ONS bundled extracts

When ONS reissues the underlying spreadsheets (e.g. the 2025 mid-year
estimates supersede the 2024 ones), regenerate the `.rda` files in
`data/` by running:

```sh
Rscript data-raw/make-paper-data.R [path/to/ons_data]
```

`data-raw/make-paper-data.R` is also the auditable record of how the
bundled extracts were derived from the public ONS files; it has no
runtime relationship to the paper scripts.

## Notes

- The simulation study (script 01) runs `nrep = 200` replicates per
  cell of a 4-target x 3-scheme grid, fitting three methods per
  replicate. Wall-clock time on a 2024-era laptop: a few minutes.
- Scripts 02 and 03 do single fits per (dataset, scheme) cell; both
  complete in seconds.
- All three scripts write a per-script log file
  (`run.log`, `ons.log`, `tail.log`) inside `out_dir`.
- The seed for the simulation study is fixed at `20250507`.
