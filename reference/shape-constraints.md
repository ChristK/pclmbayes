# Shape predicates for discrete densities

Predicates that check whether a discrete latent density `pi` (typically
a fine-grid probability vector from a PCLM fit) satisfies a particular
shape constraint:

## Usage

``` r
is_unimodal(pi, tol = 1e-12)

is_logconcave(pi, tol = 1e-10)

is_monotonic(
  pi,
  direction = c("either", "increasing", "decreasing"),
  tol = 1e-12
)
```

## Arguments

- pi:

  Numeric vector of fine-grid probabilities (length \\I\\).

- tol:

  Numerical tolerance for the inequality comparisons.

- direction:

  For `is_monotonic`: either `"either"` (the default), `"increasing"`,
  or `"decreasing"`.

## Value

Logical scalar.

## Details

- `is_unimodal`:

  Non-decreasing up to some index and non-increasing thereafter (a flat
  plateau is permitted).

- `is_logconcave`:

  \\\log \pi\_{i-1} + \log \pi\_{i+1} \le 2 \log \pi_i\\ at every
  interior index. Bins with zero probability are allowed provided they
  are not interior to a strictly positive segment.

- `is_monotonic`:

  Non-decreasing or non-increasing across the whole support, possibly
  restricted to one direction.

Each predicate uses a small numerical tolerance to absorb floating-point
noise. They are used internally by
[`bpclm`](https://christk.github.io/pclmbayes/reference/bpclm.md) to
enforce shape priors via rejection.

## Examples

``` r
is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1))   # TRUE
#> [1] TRUE
is_unimodal(c(0.1, 0.4, 0.2, 0.4, 0.1))   # FALSE (bimodal)
#> [1] FALSE
is_logconcave(dnorm(seq(-3, 3, length.out = 31)))            # TRUE
#> [1] TRUE
is_monotonic(c(0.5, 0.3, 0.1), direction = "decreasing")     # TRUE
#> [1] TRUE
```
