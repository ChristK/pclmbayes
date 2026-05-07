# Shape-constraint indicators

Returns `TRUE` if `pi` satisfies the requested shape property. Used by
[`bpclm`](https://christk.github.io/pclmbayes/reference/bpclm.md) to
enforce the rejection prior in Eq. (7) of Lambert and Eilers (2009).

## Usage

``` r
is_unimodal(pi, tol = 1e-12)
is_logconcave(pi, tol = 1e-10)
is_monotonic(pi,
             direction = c("either", "increasing", "decreasing"),
             tol = 1e-12)
```

## Arguments

- pi:

  Numeric vector of fine-grid probabilities (a discrete density).

- tol:

  Numerical tolerance.

- direction:

  For monotonicity: `"either"`, `"increasing"` or `"decreasing"`.

## Value

Logical scalar.

## Examples

``` r
is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1))
#> [1] TRUE
is_logconcave(dnorm(seq(-3, 3, length.out = 31)) /
              sum(dnorm(seq(-3, 3, length.out = 31))))
#> [1] TRUE
is_monotonic(c(0.5, 0.3, 0.1), direction = "decreasing")
#> [1] TRUE
```
