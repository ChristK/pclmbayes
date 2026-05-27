# Deaths registered in England and Wales, 2022, by single year of age

Counts of deaths registered in England and Wales during the calendar
year 2022, classified by single year of age (0, 1, ..., 104) plus an
open "105+" bin. Persons (males and females summed). Used by the paper's
ONS benchmark in `inst/paper/02_ons_benchmark.R`.

## Usage

``` r
data(ons_deaths_ew_2022)
```

## Format

A data frame with 106 rows and 4 columns:

- age_lower:

  Lower limit of the age band (years), integer.

- age_upper:

  Upper limit of the age band (years), integer. For the open "105+" bin
  this is the finite cap 110.

- age_label:

  Display label ("0", "1", ..., "104", "105+").

- count:

  Observed count in the band. Total counts equal 577,160.

## Source

Office for National Statistics. *Deaths registered in England and Wales
(series DR), 2022*. Sheets 5 (males) and 6 (females), row for 2022.
<https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsregisteredinenglandandwalesseriesdrreferencetables>

## Examples

``` r
data(ons_deaths_ew_2022)
head(ons_deaths_ew_2022)
#>   age_lower age_upper age_label count
#> 1         0         1         0  2411
#> 2         1         2         1   152
#> 3         2         3         2    95
#> 4         3         4         3    71
#> 5         4         5         4    62
#> 6         5         6         5    65
sum(ons_deaths_ew_2022$count)   # 577160
#> [1] 577160
```
