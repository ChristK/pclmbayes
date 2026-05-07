# Smoothest density that exactly preserves wide-bin counts

Solves the equality-constrained problem \$\$\min\_\phi \tfrac{1}{2}
\phi' P \phi \quad \text{s.t.} \quad C \pi(\phi) = m / m\_+,\$\$
returning the smoothest density on the fine grid that exactly reproduces
the wide-bin counts when summed. Solved by sequential quadratic
programming with the constraint linearised at each step.

## Usage

``` r
pclm_exact(m, wide_breaks,
           a = NULL, b = NULL,
           ngrid = 100L,
           ndx = 17L, degree = 3L,
           penalty_order = 3L,
           max_iter = 200L, tol = 1e-9,
           eps_ridge = 1e-8,
           phi_start = NULL,
           verbose = FALSE)
```

## Arguments

- m, wide_breaks, a, b, ngrid, ndx, degree, penalty_order:

  See [`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md).

- max_iter, tol, eps_ridge:

  Convergence and numerical stability controls for the SQP iteration.

- phi_start:

  Optional starting coefficient vector. Default uses a mildly-smoothed
  [`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md) fit.

- verbose:

  Print one line per SQP iteration.

## Value

An object of class `c("pclm_exact", "pclm")`, with components as in
[`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md), plus
`lambda` (the converged Lagrange multipliers) and `constraint_residual`
(max \\\|\gamma - m / m\_+\|\\).

## See also

[`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md),
[`calibrate`](https://christk.github.io/pclmbayes/reference/calibrate.md).

## Examples

``` r
data(bloodlead)
fit <- pclm_exact(m = bloodlead$count,
                  wide_breaks = with(bloodlead, cbind(lower, upper)),
                  a = 0, b = 80, ngrid = 80, ndx = 17,
                  degree = 3, penalty_order = 3)
fit$constraint_residual
#> [1] 2.775558e-17
max(abs(fit$fitted_counts - fit$m))
#> [1] 3.552714e-15
```
