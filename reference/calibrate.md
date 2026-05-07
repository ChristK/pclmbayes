# Calibrate a fitted PCLM to preserve wide-bin totals exactly

Post-hoc adjustment of a fitted `pclm` or `bpclm` object so that the
implied wide-bin counts match the observed counts `m` exactly. For each
wide bin \\j\\, the fine-grid probabilities \\\pi_i\\ (\\i\\ in bin
\\j\\) are multiplied by \\(m_j / m\_+) / \gamma_j\\, where \\\gamma_j =
(C\pi)\_j\\ is the fitted bin probability before calibration. The
within-bin shape from the smooth fit is preserved; small jumps at the
bin boundaries are the price paid for exact preservation.

## Usage

``` r
calibrate(fit, ...)
```

## Arguments

- fit:

  A fitted `"pclm"` or `"bpclm"` object.

- ...:

  Passed to methods (currently unused).

## Value

The same object with calibrated `pi` (and, for `bpclm`, calibrated
`pi_chain`, `pi_lower`, `pi_upper`). The element `fitted_counts` is set
equal to the input counts `m`, and `calibrated = TRUE` is added.

## Details

Calibration is meaningful only when the wide bins partition the support
without overlap, so that each fine interval belongs to exactly one wide
bin. The function checks this and errors if the bin assignment is
ambiguous (e.g. overlapping bins, or fine-grid breakpoints misaligned
with the wide-bin boundaries).

For `bpclm` objects, calibration is applied separately to every
posterior draw and the posterior summaries (`pi`, `pi_lower`,
`pi_upper`, `phi`) are recomputed. Posterior credible intervals for
derived quantities (means, quantiles) recomputed from the calibrated
chain will respect the constraint.

## See also

[`pclm_exact`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
for a constrained MAP fit that enforces the same property without
introducing kinks at bin boundaries.

## Examples

``` r
data(bloodlead)
fit <- pclm(m = bloodlead$count,
            wide_breaks = with(bloodlead, cbind(lower, upper)),
            a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
            penalty_order = 3)
max(abs(fit$fitted_counts - fit$m))            # 7.4
#> [1] 7.371174
fit_c <- calibrate(fit)
max(abs(fit_c$fitted_counts - fit_c$m))        # ~ 0
#> [1] 1.421085e-14
```
