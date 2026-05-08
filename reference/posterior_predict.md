# Posterior or parametric predictive distribution of fine-cell counts

Generates draws from the posterior or parametric predictive distribution
of the latent fine-cell (typically single-year) counts. Two flavours are
supported:

## Usage

``` r
posterior_predict(fit, ...)

# S3 method for class 'bpclm'
posterior_predict(
  fit,
  type = c("predictive", "rate"),
  level = 0.9,
  seed = NULL,
  ...
)

# S3 method for class 'pclm'
posterior_predict(
  fit,
  type = c("predictive", "rate"),
  level = 0.9,
  n_draws = 2000L,
  seed = NULL,
  ...
)

# S3 method for class 'pclm_posterior_predict'
print(x, n = 6L, ...)

# S3 method for class 'pclm_posterior_predict'
summary(object, ...)

# S3 method for class 'pclm_posterior_predict'
plot(
  x,
  show_bins = TRUE,
  xlab = "Fine-grid value",
  ylab = "Counts per cell",
  main = NULL,
  xlim = NULL,
  ylim = NULL,
  lwd = 2,
  ...
)
```

## Arguments

- fit:

  A fitted `"pclm"`, `"pclm_exact"` or `"bpclm"` object.

- ...:

  Currently unused; for future extension.

- type:

  Either `"predictive"` (the default) or `"rate"`.

- level:

  Credible/prediction level for the returned interval (default 0.9).

- seed:

  Optional integer for reproducibility.

- n_draws:

  For frequentist input objects (`"pclm"`, `"pclm_exact"`), the number
  of parametric multinomial draws to generate. Ignored for `bpclm`
  input. Default 2000.

- x:

  A `"pclm_posterior_predict"` object (in `print` and `plot` methods).

- n:

  Number of fine-grid rows shown by `print`. Default 6.

- object:

  A `"pclm_posterior_predict"` object (in `summary`).

- show_bins:

  Logical: overlay the observed wide-bin histogram on the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

- xlab, ylab, main, xlim, ylim, lwd:

  Standard graphical parameters passed to `plot`.

## Value

An object of class `"pclm_posterior_predict"`: a list with components

- `draws`:

  An `nsim x ngrid` matrix of posterior (or parametric) draws of the
  per-cell counts.

- `mean, median, lower, upper`:

  Per-cell summaries.

- `level, type`:

  The arguments used.

- `grid, grid_mid`:

  Fine-grid breakpoints and midpoints, copied from the input fit.

- `wide_breaks, m`:

  The wide-bin definitions and observed counts.

## Details

- `type = "rate"`:

  The smooth latent rate \\m\_+ \pi\\. For `bpclm` fits this returns a
  chain of rate draws (one per posterior draw of \\\pi\\). For `pclm`
  and `pclm_exact` fits this returns the point estimate \\m\_+ \pi\\ (no
  uncertainty). Mean = posterior mean rate; the credible band is narrow
  at large \\N\\.

- `type = "predictive"` (default):

  The actual realised fine-cell counts. Within each wide bin \\j\\, the
  conditional distribution of the fine-cell counts given the band total
  \\m_j\\ and the within-band probabilities \\\pi_y/\gamma_j\\ is a
  multinomial with size \\m_j\\. For `bpclm` this is the full posterior
  predictive (sample within each posterior draw of \\\pi\\); for `pclm`
  and `pclm_exact` it is a parametric multinomial bootstrap around the
  point estimate (`n_draws` samples).

Both flavours *exactly* preserve the wide-bin totals on every draw under
`type = "predictive"` (the multinomial draws sum to \\m_j\\ within each
band by construction). Under `type = "rate"` the totals are preserved
exactly only when the input fit has been calibrated
([`calibrate`](https://christk.github.io/pclmbayes/reference/calibrate.md))
or is a
[`pclm_exact`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
fit.

For an ungrouping use case (band counts \\m_j\\ treated as the data,
fine-cell counts treated as the unknown), the recommended workflow is


      fit <- bpclm(m, wide_breaks, ...)
      fit <- calibrate(fit)                              # exact band totals
      pp  <- posterior_predict(fit, type = "predictive") # plausible counts

which gives credible-interval coverage close to nominal for both the
fine-cell counts and any cumulative sum of them.

## See also

[`calibrate`](https://christk.github.io/pclmbayes/reference/calibrate.md),
[`pclm_exact`](https://christk.github.io/pclmbayes/reference/pclm_exact.md).

## Examples

``` r
# \donttest{
# Bayesian workflow with exact preservation and uncertainty:
data(bloodlead)
fit <- bpclm(m = bloodlead$count,
             wide_breaks = with(bloodlead, cbind(lower, upper)),
             a = 0, b = 80, ngrid = 80, ndx = 17,
             niter = 2000, burnin = 500, adapt = 300, seed = 1)
fit <- calibrate(fit)                                # exact band totals
pp  <- posterior_predict(fit, type = "predictive")
plot(pp)

# }
```
