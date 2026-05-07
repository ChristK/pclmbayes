# Calibrate a fitted PCLM to preserve wide-bin totals exactly

Post-hoc renormalisation of a fitted `pclm` or `bpclm` object so that
the implied wide-bin counts exactly equal the observed counts `m`. For
each wide bin \\j\\, the fine-grid probabilities \\\pi_i\\ (\\i\\ in bin
\\j\\) are scaled by \\(m_j / m\_+) / \gamma_j\\.

Calibration assumes the wide bins partition the support without overlap
so that each fine interval belongs to exactly one wide bin. For
overlapping or misaligned bins, use
[`pclm_exact`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
instead.

## Usage

``` r
calibrate(fit, ...)

# S3 method for class 'pclm'
calibrate(fit, ...)

# S3 method for class 'bpclm'
calibrate(fit, ...)
```

## Arguments

- fit:

  A fitted `"pclm"` or `"bpclm"` object.

- ...:

  Currently unused.

## Value

The same object with calibrated `pi` (and, for `bpclm`, calibrated
`pi_chain`, `pi_lower`, `pi_upper`). The component `fitted_counts` is
set to equal `m` exactly, and `calibrated = TRUE` is recorded.

## See also

[`pclm_exact`](https://christk.github.io/pclmbayes/reference/pclm_exact.md).

## Examples

``` r
data(bloodlead)
fit   <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
              penalty_order = 3)
fit_c <- calibrate(fit)
max(abs(fit_c$fitted_counts - fit_c$m))
#> [1] 1.421085e-14
```
