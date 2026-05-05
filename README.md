# pclmbayes

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
  a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3, penalty_order = 3,
  niter = 5000, burnin = 1000, adapt = 500,
  shape = "unimodal", seed = 1
)
plot(fit_b)
summary(fit_b)
```

### Exact ungrouping with uncertainty

```r
# wide-bin death counts (e.g. 5-year age bands) and band edges
fit  <- bpclm(m, wide_breaks, ...)
fit  <- calibrate(fit)                                # exact band totals
pp   <- posterior_predict(fit, type = "predictive")   # 90% PI for single-year counts
plot(pp)

# or, point-estimate only:
fit_e <- pclm_exact(m, wide_breaks, ...)
```

See the vignette (`vignette("pclmbayes-intro")`) for a full
walkthrough.

## Method summary

The latent random variable $Y$ is supported on $(a, b)$, partitioned
into $I$ fine-grid intervals of equal width $\Delta$ with midpoints
$u_i$. The log of $\pi_i \approx f_Y(u_i)\Delta$ is modelled as a
linear combination of $K$ B-splines:

$$
\pi_i = \frac{e^{\eta_i}}{\sum_{l=1}^I e^{\eta_l}}, \qquad \eta = B\phi
$$

with the identifiability constraint $\sum_k \phi_k = 0$. An $r$th-order
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
