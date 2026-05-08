# Constrained MAP fit: smoothest density that exactly preserves wide-bin counts

Solves the equality-constrained optimisation problem \$\$\min\_\phi
\tfrac{1}{2} \phi' P \phi \quad \text{subject to} \quad C \pi(\phi) = m
/ m\_+,\$\$ where \\\pi(\phi)\\ is the softmax of \\B \phi\\, \\P\\ is
the \\r\\th-order difference penalty, \\C\\ is the bin matrix, and \\m /
m\_+\\ is the vector of observed wide-bin proportions. The result is the
smoothest density on the fine grid that exactly reproduces the wide-bin
counts when summed.

## Usage

``` r
pclm_exact(
  m,
  wide_breaks,
  a = NULL,
  b = NULL,
  ngrid = 100L,
  ndx = 17L,
  degree = 3L,
  penalty_order = 3L,
  max_iter = 200L,
  tol = 1e-09,
  eps_ridge = 1e-08,
  phi_start = NULL,
  verbose = FALSE
)
```

## Arguments

- m:

  Numeric vector of non-negative wide-bin counts (length \\J\\).

- wide_breaks:

  Wide-bin boundaries; either a \\J + 1\\ vector (contiguous bins
  partitioning \\(a, b)\\) or a \\J \times 2\\ matrix/data frame of
  (lower, upper) limits.

- a, b:

  Lower and upper limits of the support \\(a, b)\\ on which to estimate
  the density. Defaults to the smallest lower limit and largest upper
  limit of `wide_breaks`.

- ngrid:

  Number of fine-grid intervals \\I\\. Defaults to 100.

- ndx:

  Number of equally-spaced knot intervals on \\(a, b)\\. The number of
  B-splines is `ndx + degree`. Default 17 (so \\K = 20\\ cubic
  B-splines, matching the paper's examples).

- degree:

  B-spline degree (default 3, cubic).

- penalty_order:

  Order \\r\\ of the difference penalty (default 3, as in the paper's
  examples).

- max_iter, tol:

  Convergence controls.

- eps_ridge:

  Small ridge added to the (singular) penalty matrix for numerical
  stability. Default `1e-8`.

- phi_start:

  Optional starting value for \\\phi\\ (length `ndx + degree`). Defaults
  to the zero vector (uniform density).

- verbose:

  Logical: if `TRUE`, print one line per scoring iteration.

## Value

An object of class `c("pclm_exact", "pclm")` with the same components as
a [`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md)
object, plus `lambda` (the converged Lagrange multipliers) and
`constraint_residual` (max \|γ - m/m\_+\|).

## Details

This is the principled alternative to
[`calibrate`](https://christk.github.io/pclmbayes/reference/calibrate.md):
no kinks are introduced at the bin boundaries because the smoothness
penalty is minimised across the whole fine grid subject to the
constraint. Solved by sequential quadratic programming (Newton steps on
the Lagrangian, with the constraint linearised at each step).

## References

Eilers, P. H. C. (2007). Ill-posed problems with counts, the composite
link model and penalized likelihood. *Statistical Modelling*, 7(3),
239–254. (Underlying penalised composite-link model on a fine grid.)

Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd ed.
Springer. Chapter 18 (Sequential Quadratic Programming, used here to
enforce the bin-total constraint exactly).

## See also

[`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md),
[`calibrate`](https://christk.github.io/pclmbayes/reference/calibrate.md).

## Examples

``` r
data(bloodlead)
fit <- pclm_exact(m = bloodlead$count,
                   wide_breaks = with(bloodlead, cbind(lower, upper)),
                   a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
                   penalty_order = 3)
max(abs(fit$fitted_counts - fit$m))    # < 1e-8
#> [1] 3.552714e-15
```
