# Latent density on the fine grid

Returns the discrete latent density values \\\pi_i = f_Y(u_i) \Delta\\
at the fine-grid midpoints given a vector of B-spline coefficients.

## Usage

``` r
latent_density(phi, basis)
```

## Arguments

- phi:

  Numeric vector of B-spline coefficients (length `basis$K`).

- basis:

  A `"pclm_basis"` object from
  [`bspline_basis`](https://christk.github.io/pclmbayes/reference/bspline_basis.md).

## Value

Numeric vector of length \\I\\ (number of fine grid intervals), summing
to 1.
