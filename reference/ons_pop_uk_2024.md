# UK mid-2024 population estimates by single year of age

Mid-2024 population estimates for the United Kingdom (persons),
classified by single year of age from 0 to 89 plus an open "90+" bin.
Extracted from the `MYE2 - Persons` sheet of the ONS mid-year population
estimates, row for the United Kingdom (code `K02000001`).

**Cross-release rounding caveat.** The "90+" count in this dataset
(625,236) is the published aggregate from the `MYE2` table. The sum of
the single-year and 105+ counts in the companion dataset
[`ons_centenarians_uk_2024`](https://christk.github.io/pclmbayes/reference/ons_centenarians_uk_2024.md)
for the same year is 625,250, differing by 14 persons (about 0.002%).
This is not a bug: the two ONS releases are rounded independently and
use slightly different bases. Reproductions of the open-tail analysis
(`inst/paper/03_open_tail.R`) should mention the discrepancy in a
footnote where the "truth" for ages 90+ is described, because anyone
re-running the figure will notice the small mismatch.

## Usage

``` r
data(ons_pop_uk_2024)
```

## Format

A data frame with 91 rows and 4 columns:

- age_lower:

  Lower limit of the age band (years), integer.

- age_upper:

  Upper limit of the age band (years), integer. For the open "90+" bin
  this is the finite cap 100.

- age_label:

  Display label ("0", "1", ..., "89", "90+").

- count:

  Observed count in the band. Total counts equal 69,281,437.

## Source

Office for National Statistics. *Estimates of the population for the UK,
England, Wales, Scotland, and Northern Ireland*, mid-2024 edition. Sheet
`MYE2 - Persons`. Released under the Open Government Licence v3.0.
<https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland>

## See also

[`ons_centenarians_uk_2024`](https://christk.github.io/pclmbayes/reference/ons_centenarians_uk_2024.md)
for the complementary single-year breakdown at ages 90 and over.

## Examples

``` r
data(ons_pop_uk_2024)
head(ons_pop_uk_2024)
#>   age_lower age_upper age_label  count
#> 1         0         1         0 667994
#> 2         1         2         1 690113
#> 3         2         3         2 732532
#> 4         3         4         3 728400
#> 5         4         5         4 755117
#> 6         5         6         5 769993
tail(ons_pop_uk_2024)    # the "90+" open bin in the last row
#>    age_lower age_upper age_label  count
#> 86        85        86        85 281337
#> 87        86        87        86 256103
#> 88        87        88        87 224742
#> 89        88        89        88 195365
#> 90        89        90        89 167231
#> 91        90       100       90+ 625236
```
