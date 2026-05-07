# Evaluate the latent density on the fine grid

Returns the discrete latent density values \\\pi_i = f_Y(u_i) \Delta\\
at the fine-grid midpoints \\u_i\\, given a vector of B-spline
coefficients.

## Usage

``` r
latent_density(phi, basis)
```

## Arguments

- phi:

  Numeric vector of B-spline coefficients of length \\K\\.

- basis:

  A `"pclm_basis"` object, as returned by
  [`bspline_basis`](https://christk.github.io/pclmbayes/reference/bspline_basis.md).

## Value

Numeric vector of length \\I\\ (number of fine grid intervals), summing
to 1.

## Examples

``` r
grid <- seq(0, 10, length.out = 51)
mids <- (head(grid, -1) + tail(grid, -1)) / 2
basis <- bspline_basis(mids, a = 0, b = 10, ndx = 17)
phi <- rep(0, basis$K)
pi_grid <- latent_density(phi, basis)
sum(pi_grid)  # 1
#> [1] 1
```
