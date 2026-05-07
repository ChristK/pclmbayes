# pclmbayes: Bayesian density estimation from grouped continuous data

Implements the penalised composite link model (PCLM) of Lambert and
Eilers (2009) for estimating a smooth latent continuous density from
grouped (binned) observations. The package provides:

## Details

- [`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md): a
  frequentist penalised-scoring fit, with selection of the smoothing
  parameter by AIC or BIC.

- [`bpclm`](https://christk.github.io/pclmbayes/reference/bpclm.md): a
  fully Bayesian fit using a modified Langevin-Hastings sampler with an
  adaptive step size and a Gibbs step on the penalty parameter.

- Optional shape-constraint priors (unimodality, log-concavity,
  monotonicity), enforced through rejection in the MCMC.

- Two example datasets used in the original paper:
  [`bloodlead`](https://christk.github.io/pclmbayes/reference/bloodlead.md)
  (Hasselblad, Stead and Galke, 1980) and
  [`tbdeaths1907`](https://christk.github.io/pclmbayes/reference/tbdeaths1907.md)
  (deaths by tuberculosis in The Netherlands in 1907).

The model assumes that the support \\(a, b)\\ of the latent random
variable \\Y\\ is partitioned into a fine grid of \\I\\ intervals of
equal width \\\Delta\\. The log-density on the fine grid is modelled
with a B-spline basis \\B\\ with \\K\\ basis functions and a discrete
\\r\\th-order difference penalty on the spline coefficients \\\phi\\
(Eilers and Marx, 1996). Letting \\\eta = B \phi\\, the latent grid
probabilities are \$\$\pi_i = \exp(\eta_i) / \sum_j \exp(\eta_j),\$\$
the wide-bin probabilities are \\\gamma = C \pi\\ where \\C\\ is a
binning matrix, and the observed wide-bin counts \\m\\ follow a
multinomial distribution with probabilities \\\gamma\\.

For full methodological details, see the package vignette.

## References

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *Computational Statistics and Data
Analysis*, 53(4), 1388–1399.
[doi:10.1016/j.csda.2008.11.022](https://doi.org/10.1016/j.csda.2008.11.022)

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89–121.

Eilers, P. H. C. (2007). Ill-posed problems with counts, the composite
link model, and penalized likelihood. *Statistical Modelling*, 7(3),
239–254.

Roberts, G. O. and Rosenthal, J. S. (1998). Optimal scaling of discrete
approximations to Langevin diffusions. *Journal of the Royal Statistical
Society, Series B*, 60(1), 255–268.

Atchade, Y. F. and Rosenthal, J. S. (2005). On adaptive Markov chain
Monte Carlo algorithms. *Bernoulli*, 11(5), 815–828.

## See also

Useful links:

- <https://github.com/ChristK/pclmbayes>

- <https://christk.github.io/pclmbayes/>

- Report bugs at <https://github.com/ChristK/pclmbayes/issues>

## Author

**Maintainer**: Chris Kypridemos <christodoulosk@gmail.com>

Authors:

- Chris Kypridemos <christodoulosk@gmail.com>

Other contributors:

- Philippe Lambert (Author of the underlying method (Lambert & Eilers,
  2009)) \[contributor\]

- Paul H. C. Eilers (Author of the underlying method (Lambert & Eilers,
  2009)) \[contributor\]
