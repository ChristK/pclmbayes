# Tuberculosis deaths by age, The Netherlands, 1907 (illustrative)

Counts of deaths attributed to tuberculosis in The Netherlands in 1907
by wide age bands. Used as Example 6.2 by Lambert and Eilers (2009) to
illustrate density estimation from heavily grouped mortality data. Total
deaths in the dataset: 9440.

**Note:** the original CBS records used in the paper were unpublished.
The values shipped here are an *illustrative reconstruction* that
preserves the published total and a plausible early-20th-century
age-mortality shape, intended for code testing and demonstration of the
methodology, not for historical inference. See `data-raw/make-data.R` in
the package source for substitution instructions if you have access to
real data.

## Usage

``` r
data(tbdeaths1907)
```

## Format

A data frame with 13 rows (12 age bands plus an upper-tail zero) and 3
columns: `lower`, `upper` (age band, years) and `count` (number of
deaths).

## References

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *CSDA*, 53(4), 1388–1399.

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
