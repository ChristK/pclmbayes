# Package index

## Package overview

- [`pclmbayes-package`](https://christk.github.io/pclmbayes/reference/pclmbayes-package.md)
  [`pclmbayes`](https://christk.github.io/pclmbayes/reference/pclmbayes-package.md)
  : Bayesian density estimation from grouped continuous data

## Frequentist fit

- [`pclm()`](https://christk.github.io/pclmbayes/reference/pclm.md) :
  Frequentist penalised composite link model for grouped data
- [`pclm_exact()`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
  : Smoothest density that exactly preserves wide-bin counts

## Bayesian fit

- [`bpclm()`](https://christk.github.io/pclmbayes/reference/bpclm.md) :
  Bayesian penalised composite link model for grouped data
- [`posterior_predict()`](https://christk.github.io/pclmbayes/reference/posterior_predict.md)
  : Posterior or parametric predictive distribution of fine-cell counts
- [`calibrate()`](https://christk.github.io/pclmbayes/reference/calibrate.md)
  : Calibrate a fitted PCLM to preserve wide-bin totals exactly

## Building blocks

- [`bspline_basis()`](https://christk.github.io/pclmbayes/reference/bspline_basis.md)
  : Equally-spaced B-spline basis on a fine grid
- [`bin_matrix()`](https://christk.github.io/pclmbayes/reference/bin_matrix.md)
  : Build the bin (composite-link) matrix C
- [`diff_penalty()`](https://christk.github.io/pclmbayes/reference/diff_penalty.md)
  : Difference penalty matrix
- [`latent_density()`](https://christk.github.io/pclmbayes/reference/latent_density.md)
  : Latent density on the fine grid

## Shape constraints

- [`is_unimodal()`](https://christk.github.io/pclmbayes/reference/shape-constraints.md)
  [`is_logconcave()`](https://christk.github.io/pclmbayes/reference/shape-constraints.md)
  [`is_monotonic()`](https://christk.github.io/pclmbayes/reference/shape-constraints.md)
  : Shape-constraint indicators

## Datasets

- [`bloodlead`](https://christk.github.io/pclmbayes/reference/bloodlead.md)
  : Lead concentration in the blood of New York children, 1974
- [`tbdeaths1907`](https://christk.github.io/pclmbayes/reference/tbdeaths1907.md)
  : Tuberculosis deaths by age, The Netherlands, 1907 (illustrative)
