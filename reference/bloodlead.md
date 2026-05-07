# Lead concentration in the blood of New York children, 1974

Interval-censored measurements of lead concentration (in \\\mu\\g/dl) in
the blood of young Puerto Ricans aged 1-12 years living in New York in
1974, as reported by Hasselblad et al. (1980) and re-used as Example 6.1
by Lambert and Eilers (2009).

## Usage

``` r
bloodlead
```

## Format

A data frame with 7 rows and 3 columns:

- lower:

  Lower limit of the wide bin (\\\mu\\g/dl).

- upper:

  Upper limit of the wide bin (\\\mu\\g/dl).

- count:

  Observed frequency \\m_j\\ in the bin.

## Details

The original instrument's resolution was limited; observations were
recorded only as falling into one of seven wide intervals. The last
interval, *65+*, was bounded above by 80 in the analysis of Lambert and
Eilers (2009). Total \\n = 139\\.

## References

Hasselblad, V., Stead, A. G. and Galke, W. (1980). Analysis of coarsely
grouped data from the lognormal distribution. *Journal of the American
Statistical Association*, 75, 771–778.

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *Computational Statistics and Data
Analysis*, 53(4), 1388–1399.

## Examples

``` r
data(bloodlead)
bloodlead
#>   lower upper count
#> 1     0    15    27
#> 2    15    25    71
#> 3    25    35    32
#> 4    35    45     6
#> 5    45    55     3
#> 6    55    65     0
#> 7    65    80     0
```
