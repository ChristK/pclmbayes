# Bayesian penalised composite link model for grouped data

Fits the model of Lambert and Eilers (2009) by Markov chain Monte Carlo.
The latent log-density is modelled by P-splines on a fine grid; wide-bin
counts have a multinomial likelihood; the spline coefficients have a
smoothness prior controlled by a precision parameter \\\tau\\ with a
vague \\\Gamma(a, b)\\ hyperprior.

## Usage

``` r
bpclm(
  m,
  wide_breaks,
  a = NULL,
  b = NULL,
  ngrid = 100L,
  ndx = 17L,
  degree = 3L,
  penalty_order = 3L,
  niter = 5000L,
  burnin = NULL,
  thin = 1L,
  adapt = 500L,
  tau_init = NULL,
  a_tau = 1e-04,
  b_tau = 1e-04,
  delta_init = NULL,
  target_accept = 0.57,
  Sigma = NULL,
  shape = NULL,
  shape_args = list(),
  phi_init = NULL,
  seed = NULL,
  verbose = FALSE
)

# S3 method for class 'bpclm'
print(x, digits = 4L, ...)

# S3 method for class 'bpclm'
summary(
  object,
  probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
  cred = 0.9,
  ...
)

# S3 method for class 'bpclm'
coef(object, ...)

# S3 method for class 'bpclm'
fitted(object, ...)

# S3 method for class 'bpclm'
plot(
  x,
  add = FALSE,
  cred = NULL,
  density_col = "black",
  band_col = grDevices::adjustcolor("steelblue", alpha.f = 0.25),
  hist_col = grDevices::adjustcolor("grey70", alpha.f = 0.6),
  hist_border = "grey40",
  xlab = "y",
  ylab = "Density",
  main = "Bayesian PCLM fit (posterior mean + credible band)",
  xlim = NULL,
  ylim = NULL,
  lwd = 2,
  ...
)

# S3 method for class 'bpclm'
predict(
  object,
  newdata,
  summary = c("mean", "median", "sample"),
  cred = NULL,
  ...
)

# S3 method for class 'bpclm'
quantile(
  x,
  probs = c(0.25, 0.5, 0.75),
  summary = c("mean", "median", "sample"),
  cred = NULL,
  ...
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

- niter:

  Total number of MCMC iterations (including burn-in). Default 5000.

- burnin:

  Number of initial iterations to discard. Default `floor(niter / 5)`.

- thin:

  Thinning interval (default 1). Only every `thin`-th post-burn-in draw
  is retained.

- adapt:

  Number of additional iterations of adaptive tuning of the step size
  `delta` performed before the main run. Default 500.

- tau_init:

  Initial value of the precision \\\tau\\. If `NULL` (the default), the
  BIC-selected value from a frequentist warm-start fit is used.

- a_tau, b_tau:

  Hyperparameters of the \\\Gamma(a, b)\\ hyperprior on \\\tau\\.
  Defaults `1e-4` each give a near- flat improper prior on \\\log
  \tau\\.

- delta_init:

  Initial step size for the Langevin proposal. If `NULL`, the value
  \\1.65^2 / K^{1/3}\\ from Roberts and Rosenthal (1998) is used.

- target_accept:

  Target acceptance rate during adaptation (default 0.57).

- Sigma:

  Optional positive-definite \\K \times K\\ proposal covariance. If
  `NULL`, the frequentist warm-start variance-covariance matrix is used.

- shape:

  Optional character vector of shape constraints. Any subset of
  `"unimodal"`, `"logconcave"` (or `"log-concave"`) and `"monotonic"`.
  Default `NULL` (unconstrained).

- shape_args:

  A list of additional arguments for the shape indicators. Currently
  used only for monotonicity, where
  `shape_args = list(direction = "decreasing")` (or `"increasing"`,
  default `"either"`) restricts the sense.

- phi_init:

  Optional starting value for the chain on \\\phi\\. Defaults to the
  warm-start frequentist MLE.

- seed:

  Optional integer for reproducibility.

- verbose:

  Logical: print progress every 10% of iterations.

- probs:

  Probabilities at which to compute quantiles.

- cred:

  Credible level for the band (default `x$cred_level`). For
  `quantile.bpclm`, also reports lower/upper credible limits when
  `summary = "mean"`.

- add:

  Logical: if `TRUE`, add to an existing plot.

- density_col:

  Colour for the fitted density line.

- band_col:

  Fill colour for the credible band (alpha-blended).

- hist_col:

  Fill colour for the histogram rectangles.

- hist_border:

  Border colour for histogram rectangles.

- xlab, ylab, main, ylim, xlim, lwd:

  Standard graphical parameters (passed to plotting methods of `pclm` /
  `bpclm` fits).

- summary:

  One of `"mean"` (posterior mean of the quantile or predicted density),
  `"median"`, or `"sample"` (returns the full posterior sample as a
  matrix).

## Value

An object of class `"bpclm"`, a list with components:

- phi:

  Posterior mean of \\\phi\\ (length \\K\\).

- phi_chain:

  `nsim x K` matrix of post-burn-in, post-thin draws of \\\phi\\.

- tau_chain:

  Numeric vector of post-burn-in draws of \\\tau\\.

- pi_chain:

  `nsim x I` matrix of fine-grid probabilities (each row a draw).

- pi:

  Posterior mean of \\\pi\\ (latent grid density, summing to 1).

- pi_lower, pi_upper:

  Pointwise 90% (or other, see `cred_level`) credible-interval limits
  for \\\pi\\.

- cred_level:

  Credible level used for the bands (default 0.90).

- accept:

  Overall acceptance rate of the \\\phi\\ step across post-adaptation
  iterations.

- delta:

  Final tuned step size.

- adapt_path:

  Numeric vector tracking \\\delta\\ over the adaptation phase (for
  diagnostics).

- warmstart:

  The frequentist `pclm` fit used for warm starting (or `NULL` if all
  relevant arguments were supplied directly).

- grid, grid_mid, basis, C, m, wide_breaks, penalty_order, a_tau, b_tau,
  shape, call:

  Bookkeeping.

## Details

The sampler is a Metropolis-within-Gibbs scheme: a modified
Langevin-Hastings (MALA) update for \\\phi\\ followed by a Gibbs draw
for \\\tau\\. The proposal covariance for \\\phi\\ is a fixed matrix
`Sigma` approximating the posterior covariance, obtained from a
frequentist warm-start fit
([`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md)) at a
user-chosen \\\tau\\.

Optional shape constraints (`"unimodal"`, `"logconcave"` and/or
`"monotonic"`) are imposed by rejecting any proposal whose induced
density violates them, in line with Eq. (7) of the paper.

## References

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *Computational Statistics and Data
Analysis*, 53(4), 1388–1399.

Roberts, G. O. and Rosenthal, J. S. (1998). Optimal scaling of discrete
approximations to Langevin diffusions. *Journal of the Royal Statistical
Society, Series B*, 60(1), 255–268.

Atchade, Y. F. and Rosenthal, J. S. (2005). On adaptive Markov chain
Monte Carlo algorithms. *Bernoulli*, 11(5), 815–828.

## Examples

``` r
# \donttest{
data(bloodlead)
fit <- bpclm(m = bloodlead$count,
             wide_breaks = with(bloodlead, cbind(lower, upper)),
             a = 0, b = 80,
             ngrid = 80, ndx = 17,
             niter = 2000, burnin = 500,
             shape = "unimodal", seed = 1)
summary(fit)
#> Bayesian penalised composite link model
#> Call: bpclm(m = bloodlead$count, wide_breaks = with(bloodlead, cbind(lower, 
#>     upper)), a = 0, b = 80, ngrid = 80, ndx = 17, niter = 2000, 
#>     burnin = 500, shape = "unimodal", seed = 1)
#> 
#> Number of wide bins: 7  | total counts: 139 
#> Fine grid:80intervals on (0, 80)
#> B-spline basis: K = 20 (degree = 3 )
#> Penalty order r = 3 
#> MCMC: niter = 2000, burnin = 500, thin = 1, kept = 1500
#> Final delta = 1.09  |  acceptance rate = 0.63
#> Shape constraint(s): unimodal 
#> Posterior mean tau = 19.5 (sd 46.7)
#> 
#> Posterior of mean(Y):  21.7779  (90% CI: 20.5296, 23.0597)
#> Posterior of sd(Y):    8.2889  (90% CI: 7.1929, 9.3911)
#> 
#> Posterior summaries of quantiles (mean and 90% CI):
#>     p    mean      lo      hi
#>  0.05  9.6024  7.0299 11.6007
#>  0.10 12.0372 10.1232 13.6278
#>  0.25 16.1954 14.9382 17.4892
#>  0.50 20.9888 19.7886 22.4045
#>  0.75 26.4068 24.8010 28.1725
#>  0.90 32.5078 30.0854 35.0970
#>  0.95 36.9349 33.7977 40.6505
# }
```
