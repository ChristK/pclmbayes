# Bayesian density estimation from grouped continuous data

Implements the penalised composite link model (PCLM) of Lambert and
Eilers (2009) for estimating a smooth latent continuous density from
grouped (binned) observations. Provides both a frequentist
penalised-scoring fit
([`pclm`](https://christk.github.io/pclmbayes/reference/pclm.md)) and a
fully Bayesian fit
([`bpclm`](https://christk.github.io/pclmbayes/reference/bpclm.md))
using a modified Langevin-Hastings sampler with an adaptive step size
and a Gibbs step on the smoothing parameter. Optional shape-constraint
priors (unimodality, log-concavity, monotonicity) are supported.

## References

Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. *Computational Statistics and Data
Analysis*, 53(4), 1388–1399.

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89–121.

Roberts, G. O. and Rosenthal, J. S. (1998). Optimal scaling of discrete
approximations to Langevin diffusions. *Journal of the Royal Statistical
Society, Series B*, 60(1), 255–268.

Atchade, Y. F. and Rosenthal, J. S. (2005). On adaptive Markov chain
Monte Carlo algorithms. *Bernoulli*, 11(5), 815–828.
