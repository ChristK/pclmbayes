# Frequentist penalised composite link model for grouped data

Fits the model of Lambert and Eilers (2009) by penalised maximum
likelihood. The log of the latent density is modelled as a linear
combination of B-splines with an \\r\\th-order difference penalty on the
coefficients. Wide-bin counts `m` are assumed to follow a multinomial
distribution with probabilities \\\gamma = C \pi(\phi)\\.

## Usage

``` r
# S3 method for class 'pclm'
print(x, digits = 4L, ...)

# S3 method for class 'pclm'
summary(object, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), ...)

# S3 method for class 'pclm'
coef(object, ...)

# S3 method for class 'pclm'
fitted(object, ...)

# S3 method for class 'pclm'
logLik(object, ...)

# S3 method for class 'pclm'
plot(
  x,
  add = FALSE,
  density_col = "black",
  hist_col = grDevices::adjustcolor("grey70", alpha.f = 0.6),
  hist_border = "grey40",
  xlab = "y",
  ylab = "Density",
  main = "Frequentist PCLM fit",
  xlim = NULL,
  ylim = NULL,
  lwd = 2,
  ...
)

# S3 method for class 'pclm'
predict(object, newdata, ...)

# S3 method for class 'pclm'
quantile(x, probs = c(0.25, 0.5, 0.75), ...)

pclm(
  m,
  wide_breaks,
  a = NULL,
  b = NULL,
  ngrid = 100L,
  ndx = 17L,
  degree = 3L,
  penalty_order = 3L,
  tau = NULL,
  select = c("BIC", "AIC"),
  max_iter = 100L,
  tol = 1e-07,
  phi_start = NULL,
  verbose = FALSE
)
```

## Arguments

- probs:

  Probabilities at which to compute quantiles.

- add:

  Logical: if `TRUE`, add to an existing plot.

- density_col:

  Colour for the fitted density line.

- hist_col:

  Fill colour for the histogram rectangles.

- hist_border:

  Border colour for histogram rectangles.

- xlab, ylab, main, ylim, xlim, lwd:

  Standard graphical parameters (passed to plotting methods of `pclm` /
  `bpclm` fits).

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

- tau:

  Either a positive scalar smoothing parameter, or a numeric vector of
  candidate values to evaluate. If `NULL` (the default), a
  logarithmically-spaced grid from `1e-2` to `1e6` (length 25) is used.

- select:

  Information criterion used to pick \\\tau\\ when a grid is provided:
  either `"BIC"` (default) or `"AIC"`.

- max_iter:

  Maximum number of scoring iterations.

- tol:

  Convergence tolerance on the largest absolute change in \\\phi\\.

- phi_start:

  Optional starting value for \\\phi\\ (length `ndx + degree`). Defaults
  to the zero vector (uniform density).

- verbose:

  Logical: if `TRUE`, print one line per scoring iteration.

## Value

An object of class `"pclm"`, a list with components:

- phi:

  Estimated B-spline coefficients (length \\K\\, summing to 0).

- tau:

  Selected (or supplied) smoothing parameter.

- tau_grid:

  Numeric vector: the full grid of \\\tau\\ values evaluated. `NULL` if
  a single \\\tau\\ was passed.

- ic:

  Information-criterion value at the selected \\\tau\\.

- ic_grid:

  Vector of IC values at the candidate \\\tau\\ values (`NULL` if a
  single \\\tau\\).

- select:

  Which IC was used (`"BIC"` or `"AIC"`).

- logL:

  Log-likelihood at convergence.

- edf:

  Effective degrees of freedom, \\\mathrm{tr}((I + \tau P)^{-1} I)\\.

- vcov:

  Estimated variance-covariance matrix of \\\phi\\ (Eq. 3 of the paper).

- pi:

  Latent grid probabilities \\\pi_i\\.

- gamma:

  Fitted wide-bin probabilities \\\gamma_j\\.

- fitted_counts:

  \\m\_+ \gamma\\, the multinomial expectations.

- grid:

  Fine-grid breakpoints (length \\I + 1\\).

- grid_mid:

  Fine-grid midpoints (length \\I\\).

- basis:

  The `"pclm_basis"` object used.

- C:

  The bin matrix \\C\\.

- m, wide_breaks:

  The supplied data.

- call:

  The matched call.

## Details

For a given smoothing parameter \\\tau\\, the fit is obtained by the
penalised scoring algorithm (Eq. 2 of the paper). When a grid of
candidate \\\tau\\ values is supplied via `tau` (a numeric vector of
length \\\> 1\\), the value minimising the chosen information criterion
(`select = "BIC"`, the default, or `"AIC"`) is selected. The
variance-covariance matrix of the coefficients is the inverse of the
penalised Fisher information at convergence (Eq. 3 of the paper).

## References

Eilers, P. H. C. (2007). Ill-posed problems with counts, the composite
link model and penalized likelihood. *Statistical Modelling*, 7(3),
239–254. (Original frequentist penalised PCLM.)

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *Computational Statistics and Data
Analysis*, 53(4), 1388–1399. (Section 2, Eq. (2): the specific
penalised-scoring form implemented here.)

## See also

[`bpclm`](https://christk.github.io/pclmbayes/reference/bpclm.md) for
the Bayesian variant.

## Examples

``` r
data(bloodlead)
fit <- pclm(m = bloodlead$count,
            wide_breaks = with(bloodlead, cbind(lower, upper)),
            a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
            penalty_order = 3)
summary(fit)
#> Penalised composite link model (frequentist)
#> Call: pclm(m = bloodlead$count, wide_breaks = with(bloodlead, cbind(lower, 
#>     upper)), a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3, 
#>     penalty_order = 3)
#> 
#> Number of wide bins: 7  | total counts: 139 
#> Fine grid:80intervals on (0, 80)
#> B-spline basis: K = 20 (degree = 3 )
#> Penalty order r = 3 
#> Selected tau = 100 (BIC = 355.79, edf = 2.56)
#> Converged = TRUE in 13 iterations.  log-likelihood = -171.582
#> 
#> Fitted summary statistics of the latent density:
#>   mean = 21.7395, sd = 8.5228
#>      5%     10%     25%     50%     75%     90%     95% 
#>  8.6085 11.2664 15.8611 21.2147 26.9966 32.8015 36.6231 
#> 
#> Goodness of fit (observed vs fitted counts):
#>  lower upper obs   exp
#>      0    15  27 29.96
#>     15    25  71 63.63
#>     25    35  32 36.04
#>     35    45   6  8.12
#>     45    55   3  1.13
#>     55    65   0  0.11
#>     65    80   0  0.01
```
