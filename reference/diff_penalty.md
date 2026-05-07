# Difference penalty matrix

Returns the matrix \\P = D'D\\ that gives the sum of squared
\\r\\th-order differences of a vector of length \\K\\, i.e. \\\sum_k
(\Delta^r \phi_k)^2 = \phi' P \phi\\. The matrix \\D\\ itself can be
returned via `return_D = TRUE`.

## Usage

``` r
diff_penalty(K, r = 3L, return_D = FALSE)
```

## Arguments

- K:

  Length of the coefficient vector (number of B-splines).

- r:

  Penalty order (default 3, as recommended by Lambert and Eilers 2009
  because it tends to a normal density at the limit of strong
  smoothing).

- return_D:

  Logical. If `TRUE`, the function returns a list with components `D`
  (the \\(K-r) \times K\\ difference matrix) and `P` (the \\K \times K\\
  penalty matrix). If `FALSE` (the default), only `P` is returned as a
  matrix.

## Value

Either a matrix `P` or a list (see `return_D`).

## Details

The rank of \\P\\ is \\K - r\\; its null space is spanned by polynomials
of degree up to \\r-1\\ (in the sequence index). This rank-deficiency is
the reason for the identifiability constraint \\\sum_k \phi_k = 0\\ in
the model.

## Examples

``` r
diff_penalty(K = 6, r = 2)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1   -2    1    0    0    0
#> [2,]   -2    5   -4    1    0    0
#> [3,]    1   -4    6   -4    1    0
#> [4,]    0    1   -4    6   -4    1
#> [5,]    0    0    1   -4    5   -2
#> [6,]    0    0    0    1   -2    1
```
