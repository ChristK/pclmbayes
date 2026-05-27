# UK 2024 population at very old ages

Mid-2024 population estimates for the United Kingdom (persons),
classified by single year of age from 90 to 104 plus an open "105+" bin.
ONS publishes this separately from the main `MYE2` table because `MYE2`
closes the age distribution at "90+". The 105+ count (610) is the input
the open-tail analysis relies on (see `inst/paper/03_open_tail.R`).

**Cross-release rounding caveat.** The total of this dataset
(single-year counts at ages 90 to 104 plus the 105+ count) is 625,250,
which differs from the `MYE2` "90+" aggregate of 625,236 in
[`ons_pop_uk_2024`](https://christk.github.io/pclmbayes/reference/ons_pop_uk_2024.md)
by 14 persons (about 0.002%). This is not a bug: the two ONS releases
are rounded independently. Manuscripts reproducing the open-tail figure
should mention the discrepancy in a footnote where the "truth" at ages
90+ is described.

## Usage

``` r
data(ons_centenarians_uk_2024)
```

## Format

A data frame with 16 rows and 4 columns:

- age_lower:

  Lower limit of the age band (years), integer.

- age_upper:

  Upper limit of the age band (years), integer. For the open "105+" bin
  this is the finite cap 110.

- age_label:

  Display label ("90", "91", ..., "104", "105+").

- count:

  Observed count in the band.

## Source

Office for National Statistics. *Mid-year population estimates of the
very old, including centenarians: UK*, edition covering 2002 to 2024.
Row for Persons, 2024. Released under the Open Government Licence v3.0.
<https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/ageing/datasets/midyearpopulationestimatesoftheveryoldincludingcentenariansunitedkingdom>

## See also

[`ons_pop_uk_2024`](https://christk.github.io/pclmbayes/reference/ons_pop_uk_2024.md)
for the complementary single-year breakdown at ages 0 to 89.

## Examples

``` r
data(ons_centenarians_uk_2024)
ons_centenarians_uk_2024
#>    age_lower age_upper age_label  count
#> 1         90        91        90 136960
#> 2         91        92        91 114880
#> 3         92        93        92  95850
#> 4         93        94        93  76860
#> 5         94        95        94  59550
#> 6         95        96        95  44020
#> 7         96        97        96  31760
#> 8         97        98        97  22670
#> 9         98        99        98  15650
#> 10        99       100        99  10390
#> 11       100       101       100   6720
#> 12       101       102       101   4260
#> 13       102       103       102   2620
#> 14       103       104       103   1560
#> 15       104       105       104    890
#> 16       105       110      105+    610
```
