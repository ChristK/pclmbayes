# Changelog

## pclmbayes 0.2.0

### Bug fix: erratic fits when the basis dimension is close to the bin count

- [`pclm()`](https://christk.github.io/pclmbayes/reference/pclm.md),
  [`bpclm()`](https://christk.github.io/pclmbayes/reference/bpclm.md)
  and
  [`pclm_exact()`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
  could return wildly erratic densities — most visibly large spurious
  excursions at the edges of the support — whenever the number of
  B-splines happened to fall near the number of wide bins .

  The cause is structural, not numerical. The composite-link matrix is ,
  so the wide-bin counts identify at most directions in the coefficient
  vector ; the rest are fixed by alone. For those directions are
  high-frequency, where is large and even a small suppresses them. For
  they are smooth, where is small, the Fisher information is severely
  ill-conditioned, and can grow without bound. A *larger* basis is
  therefore more stable, not less.

  Concretely: ungrouping the ONS mid-2024 UK population from 5-year
  bands plus an open 90+ bin () back to single years of age gave a
  fitted 2,365,717 people at age 0 against a true 667,994 — a mean
  absolute error of 30% across ages 0–20, peaking at 254%, with . With
  the new default the same call gives 1.0% mean / 2.5% maximum and .
  [`bpclm()`](https://christk.github.io/pclmbayes/reference/bpclm.md)
  was affected identically: the Gibbs step on does not protect against
  this, because the posterior collapses towards zero for the same reason
  the BIC-selected value does.

- **Breaking change.** The default `ndx` in
  [`pclm()`](https://christk.github.io/pclmbayes/reference/pclm.md),
  [`bpclm()`](https://christk.github.io/pclmbayes/reference/bpclm.md)
  and
  [`pclm_exact()`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
  is now `NULL`, meaning the basis dimension is derived from the problem
  as `K = min(max(J + 7, 20), ngrid, 200)` rather than being fixed at
  `K = 20`. This keeps `K = 20` — the value used in the examples of
  Lambert and Eilers (2009) — for every problem with , so fits of
  `bloodlead` () and `tbdeaths1907` () are bit-for-bit unchanged.
  Problems with more wide bins now get a larger basis. Pass `ndx`
  explicitly to restore the old behaviour exactly.

- New argument `check_basis` (default `TRUE`) in
  [`pclm()`](https://christk.github.io/pclmbayes/reference/pclm.md),
  [`bpclm()`](https://christk.github.io/pclmbayes/reference/bpclm.md)
  and
  [`pclm_exact()`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
  warns when the basis dimension falls inside the weak-identification
  band , naming the constraint that binds (`ndx` or `ngrid`). Set
  `check_basis = FALSE` to silence it.

- New section **Choosing the basis dimension** in
  [`?pclm`](https://christk.github.io/pclmbayes/reference/pclm.md) and
  [`?pclm_exact`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
  explains the mechanism, the safe threshold, and the related caution
  that at very large the multinomial log-likelihood dwarfs the
  effective-degrees-of-freedom term, so `BIC` and `AIC` can be monotone
  in across the whole candidate grid and the selected is pinned at an
  endpoint. Check with `fit$tau %in% range(fit$tau_grid)`;
  [`pclm_exact()`](https://christk.github.io/pclmbayes/reference/pclm_exact.md)
  avoids the issue entirely, having no to select.

## pclmbayes 0.1.0

- Initial release. Implements the method of Lambert and Eilers (2009)
  “Bayesian density estimation from grouped continuous data”
  (*Computational Statistics and Data Analysis* 53(4), 1388-1399).
- Frequentist fit
  [`pclm()`](https://christk.github.io/pclmbayes/reference/pclm.md):
  penalised scoring with BIC/AIC selection of the smoothing parameter.
- Bayesian fit
  [`bpclm()`](https://christk.github.io/pclmbayes/reference/bpclm.md):
  modified Langevin-Hastings sampler with rotation by the Cholesky of a
  frequentist warm-start variance- covariance, adaptive tuning of the
  step size, and a Gibbs step on the smoothing precision.
- Optional shape-constraint priors: unimodality, log-concavity and
  monotonicity (Eq. 7 of the paper).
- Example datasets `bloodlead` (verbatim from Hasselblad et al. 1980 /
  the paper) and `tbdeaths1907` (illustrative reconstruction; see help
  page for caveats).
- Vignette reproducing the worked examples of the paper.
