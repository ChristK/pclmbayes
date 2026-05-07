# Equally-spaced B-spline basis on a fine grid

Builds a B-spline basis evaluated at `x`, with equally-spaced knots
extending beyond the support to give partition-of-unity at all interior
points (the standard P-spline construction of Eilers and Marx, 1996).
The number of basis functions is `K = ndx + degree`.

## Usage

``` r
bspline_basis(x, a, b, ndx = 17L, degree = 3L)
```

## Arguments

- x:

  Numeric vector of evaluation points.

- a, b:

  Lower and upper limits of the support \\(a, b)\\.

- ndx:

  Number of equally-spaced knot intervals on \\(a, b)\\.

- degree:

  Spline degree (default 3 for cubic B-splines).

## Value

A list of class `"pclm_basis"` with components `B` (the basis matrix),
`knots`, `a`, `b`, `ndx`, `degree`, `K`.

## References

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89–121.

## Examples

``` r
mids <- seq(0.05, 9.95, by = 0.1)
bs   <- bspline_basis(mids, a = 0, b = 10, ndx = 17, degree = 3)
dim(bs$B)
#> [1] 100  20
```
