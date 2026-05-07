# Build an equally-spaced B-spline basis on a fine grid

Constructs a B-spline basis \\B\\ of dimension \\I \times K\\ evaluated
at the midpoints \\u_i\\ of the fine-grid intervals partitioning \\(a,
b)\\. The knots are equally spaced and extend beyond the support so that
the basis is well-defined at the boundaries (the standard P-spline
construction of Eilers and Marx, 1996). The number of basis functions is
\\K = ndx + degree\\, where `ndx` is the number of knot intervals inside
\\(a, b)\\ and `degree` is the spline degree (cubic by default).

## Usage

``` r
bspline_basis(x, a, b, ndx = 17L, degree = 3L)
```

## Arguments

- x:

  Numeric vector of evaluation points (typically the fine-grid
  midpoints).

- a, b:

  Lower and upper limits of the support \\(a, b)\\.

- ndx:

  Number of equally-spaced knot intervals on \\(a, b)\\. The number of
  B-splines is `K = ndx + degree`.

- degree:

  Spline degree (default 3 for cubic B-splines).

## Value

An object of class `"pclm_basis"`, a list with components

- B:

  Numeric matrix of size `length(x) x K` with the basis evaluated at
  `x`.

- knots:

  The augmented knot sequence used internally.

- a, b, ndx, degree, K:

  The arguments and derived dimension.

## References

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89–121.

## Examples

``` r
grid <- seq(0, 10, length.out = 100)
basis <- bspline_basis(grid, a = 0, b = 10, ndx = 17, degree = 3)
dim(basis$B)  # 100 x 20
#> [1] 100  20
```
