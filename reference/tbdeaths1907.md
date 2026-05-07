# Tuberculosis deaths by age, The Netherlands, 1907

Counts of deaths attributed to tuberculosis in The Netherlands in 1907,
classified by wide age bands. Used as Example 6.2 by Lambert and Eilers
(2009) to illustrate density estimation from heavily grouped mortality
data. Total deaths in the dataset: 9440.

## Usage

``` r
tbdeaths1907
```

## Format

A data frame with 12 rows and 3 columns:

- lower:

  Lower limit of the age band (years).

- upper:

  Upper limit of the age band (years).

- count:

  Number of deaths attributed to tuberculosis in the band.

## Details

Yearly population numbers in single-year age intervals were available,
but TB deaths were only recorded in wide, irregular age bands (a common
feature of historical mortality data). An extra interval \[100, 120)
with zero count is added at the upper tail to force the estimated
density to taper smoothly to zero, as suggested by Lambert and Eilers
(2009, Section 4).

## References

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *Computational Statistics and Data
Analysis*, 53(4), 1388–1399.

## Examples

``` r
data(tbdeaths1907)
tbdeaths1907
#>    lower upper count
#> 1      0     1   540
#> 2      1     5   980
#> 3      5    10   320
#> 4     10    15   240
#> 5     15    20   580
#> 6     20    30  1880
#> 7     30    40  1690
#> 8     40    50  1230
#> 9     50    60   850
#> 10    60    70   590
#> 11    70    80   380
#> 12    80   100   160
#> 13   100   120     0
```
