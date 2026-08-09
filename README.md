# pclmbayes

<!-- badges: start -->
[![R-CMD-check](https://github.com/ChristK/pclmbayes/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ChristK/pclmbayes/actions/workflows/R-CMD-check.yml)
[![test-coverage](https://github.com/ChristK/pclmbayes/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/ChristK/pclmbayes/actions/workflows/test-coverage.yml)
[![Codecov test coverage](https://codecov.io/gh/ChristK/pclmbayes/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ChristK/pclmbayes)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Bayesian density estimation and exact ungrouping from grouped
continuous data, via the penalised composite link model**

`pclmbayes` covers two related problems:

1. **Density estimation** — given wide-bin counts that are a
   multinomial sample of an underlying continuous variable, recover
   the smooth latent density (Lambert and Eilers, 2009).
2. **Exact ungrouping** — given wide-bin counts that *are* the data
   (e.g. deaths in 5-year age bands), distribute them across a fine
   grid (e.g. single years of age) such that the grid totals sum back
   exactly to the band totals.

The package provides:

- `pclm()` / `bpclm()` — frequentist (BIC) and Bayesian fits of the
  Lambert–Eilers density model, with optional shape constraints
  (`unimodal`, `logconcave`, `monotonic`).
- `pclm_exact()` — constrained MAP that gives the smoothest density
  consistent with the wide-bin counts *exactly*.
- `calibrate()` — post-hoc projection of any fit (frequentist or
  Bayesian) onto the constraint manifold; for Bayesian fits, every
  posterior draw is calibrated.
- `posterior_predict()` — per-cell credible/prediction intervals for
  the ungrouped fine-cell counts, with multinomial sampling noise.
- Datasets `bloodlead` (Hasselblad et al. 1980) and `tbdeaths1907`
  (illustrative reconstruction) used in Section 6 of the paper.
- A full vignette walking through all of the above.

## Installation

### From GitHub (recommended)

```r
# install.packages("remotes")
remotes::install_github("ChristK/pclmbayes")

# To also build the vignettes locally:
remotes::install_github("ChristK/pclmbayes", build_vignettes = TRUE)
```

A specific branch, tag, or commit can be installed with the `ref`
argument, e.g. `remotes::install_github("ChristK/pclmbayes", ref = "main")`.

### From a local clone

The package source is self-contained: the `.rda` data files and the
`man/` documentation files are bundled, so no roxygen2 round-trip is
required to install. From the directory containing the package:

```r
install.packages("pclmbayes", repos = NULL, type = "source")
# or, with devtools:
# devtools::install("/path/to/pclmbayes")
```

To replace `tbdeaths1907` with real CBS data (or otherwise rebuild the
data objects), edit and re-run `data-raw/make-data.R`:

```r
source("data-raw/make-data.R")
```

## Quick start

### Density estimation (Lambert–Eilers paradigm)

```r
library(pclmbayes)
data(bloodlead)

fit_b <- bpclm(
  m           = bloodlead$count,
  wide_breaks = with(bloodlead, cbind(lower, upper)),
  a = 0, b = 80, ngrid = 80, degree = 3, penalty_order = 3,
  niter = 5000, burnin = 1000, adapt = 500,
  shape = "unimodal", seed = 1
)
plot(fit_b)
summary(fit_b)
```

### Exact ungrouping with uncertainty

```r
library(pclmbayes)
data(tbdeaths1907)   # wide-bin (uneven) death counts + band edges

m           <- tbdeaths1907$count
wide_breaks <- cbind(tbdeaths1907$lower, tbdeaths1907$upper)

# Bayesian ungrouping to single years of age on (0, 120).
fit <- bpclm(
  m, wide_breaks,
  a = 0, b = 120, ngrid = 120,
  niter = 5000, burnin = 1000, seed = 1
)
fit <- calibrate(fit)                                 # exact band totals
pp  <- posterior_predict(fit, type = "predictive")    # 90% PI for single-year counts
plot(pp)

# or, point-estimate only:
fit_e <- pclm_exact(m, wide_breaks, a = 0, b = 120, ngrid = 120)
```

See the vignette (`vignette("pclmbayes-intro")`) for a full
walkthrough.

### A note on the basis dimension

The wide-bin counts identify at most $J - 1$ directions in the
coefficient vector $\phi$, where $J$ is the number of wide bins; the
rest are set by the penalty alone. When the number of B-splines $K$ is
close to $J$ those directions are *smooth*, a difference penalty barely
resists them, and the fit can develop large spurious excursions at the
edges of the support. A larger basis is therefore more stable, not less.

Leaving `ndx` at its default is the safe choice: $K$ is then derived
from the problem as `min(max(J + 7, 20), ngrid, 200)`. If you set `ndx`
yourself, keep $K =$ `ndx + degree` at or above $J + 4$; `pclm()`,
`bpclm()` and `pclm_exact()` warn otherwise. See
`?pclm`, section *Choosing the basis dimension*.

## Method summary

The latent random variable $Y$ is supported on $(a, b)$, partitioned
into $I$ fine-grid intervals of equal width $\Delta$ with midpoints
$u_i$. The log of $\pi_i \approx f_Y(u_i)\Delta$ is modelled as a
linear combination of $K$ B-splines:

$$
\pi_i = \frac{e^{\eta_i}}{\sum_{l=1}^I e^{\eta_l}}, \qquad \eta = B\phi
$$

with the identifiability constraint $\sum_k \phi_k = 0$. An $r$-th order
discrete difference penalty on the spline coefficients $\phi$ provides
smoothing. The wide-bin probabilities $\gamma = C\pi$ are linked to
the observed counts $m$ through a multinomial likelihood. The
Bayesian variant places a Gaussian prior on $\Delta^r\phi$ with
precision $\tau$, with a vague Gamma hyperprior on $\tau$.

For full methodological details see the package vignette and the
original paper.

## Dataset notes

* `bloodlead` is the verbatim seven-bin frequency table from
  Hasselblad, Stead and Galke (1980) reproduced in Lambert and Eilers
  (2009, Section 6.1). The "65+" bin is bounded above by 80 µg/dl as
  in the paper.
* `tbdeaths1907` is an **illustrative reconstruction** of the dataset
  used in Section 6.2 of the paper. The original CBS records were
  unpublished; the reconstruction matches the published total of 9440
  deaths and a plausible early-20th-century age-mortality shape, and
  is intended for demonstration of the methodology rather than for
  historical inference. See `data-raw/make-data.R` for substitution
  instructions.

## Acknowledgements

This package implements a method developed by **Philippe Lambert** and
**Paul H. C. Eilers** (Lambert and Eilers, 2009). We gratefully
acknowledge their foundational work, on which `pclmbayes` builds. Any
errors in this implementation are our own.

## References

- Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
  from grouped continuous data. *Computational Statistics and Data
  Analysis*, 53(4), 1388–1399.
  doi:[10.1016/j.csda.2008.11.022](https://doi.org/10.1016/j.csda.2008.11.022)
- Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
  B-splines and penalties. *Statistical Science*, 11(2), 89–121.
- Roberts, G. O. and Rosenthal, J. S. (1998). Optimal scaling of
  discrete approximations to Langevin diffusions. *JRSS B*, 60(1),
  255–268.
- Atchadé, Y. F. and Rosenthal, J. S. (2005). On adaptive Markov chain
  Monte Carlo algorithms. *Bernoulli*, 11(5), 815–828.

## License

GPL (>= 3).
