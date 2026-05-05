#' pclmbayes: Bayesian density estimation from grouped continuous data
#'
#' Implements the penalised composite link model (PCLM) of Lambert and
#' Eilers (2009) for estimating a smooth latent continuous density from
#' grouped (binned) observations. The package provides:
#'
#' \itemize{
#'   \item \code{\link{pclm}}: a frequentist penalised-scoring fit, with
#'     selection of the smoothing parameter by AIC or BIC.
#'   \item \code{\link{bpclm}}: a fully Bayesian fit using a modified
#'     Langevin-Hastings sampler with an adaptive step size and a Gibbs
#'     step on the penalty parameter.
#'   \item Optional shape-constraint priors (unimodality, log-concavity,
#'     monotonicity), enforced through rejection in the MCMC.
#'   \item Two example datasets used in the original paper:
#'     \code{\link{bloodlead}} (Hasselblad, Stead and Galke, 1980) and
#'     \code{\link{tbdeaths1907}} (deaths by tuberculosis in The
#'     Netherlands in 1907).
#' }
#'
#' The model assumes that the support \eqn{(a, b)} of the latent random
#' variable \eqn{Y} is partitioned into a fine grid of \eqn{I} intervals
#' of equal width \eqn{\Delta}. The log-density on the fine grid is
#' modelled with a B-spline basis \eqn{B} with \eqn{K} basis functions
#' and a discrete \eqn{r}th-order difference penalty on the spline
#' coefficients \eqn{\phi} (Eilers and Marx, 1996). Letting
#' \eqn{\eta = B \phi}, the latent grid probabilities are
#' \deqn{\pi_i = \exp(\eta_i) / \sum_j \exp(\eta_j),}
#' the wide-bin probabilities are \eqn{\gamma = C \pi} where \eqn{C} is
#' a binning matrix, and the observed wide-bin counts \eqn{m} follow a
#' multinomial distribution with probabilities \eqn{\gamma}.
#'
#' For full methodological details, see the package vignette.
#'
#' @references
#' Lambert, P. and Eilers, P. H. C. (2009).
#' Bayesian density estimation from grouped continuous data.
#' \emph{Computational Statistics and Data Analysis}, 53(4),
#' 1388--1399. \doi{10.1016/j.csda.2008.11.022}
#'
#' Eilers, P. H. C. and Marx, B. D. (1996).
#' Flexible smoothing with B-splines and penalties.
#' \emph{Statistical Science}, 11(2), 89--121.
#'
#' Eilers, P. H. C. (2007).
#' Ill-posed problems with counts, the composite link model, and
#' penalized likelihood. \emph{Statistical Modelling}, 7(3), 239--254.
#'
#' Roberts, G. O. and Rosenthal, J. S. (1998).
#' Optimal scaling of discrete approximations to Langevin diffusions.
#' \emph{Journal of the Royal Statistical Society, Series B}, 60(1),
#' 255--268.
#'
#' Atchade, Y. F. and Rosenthal, J. S. (2005).
#' On adaptive Markov chain Monte Carlo algorithms.
#' \emph{Bernoulli}, 11(5), 815--828.
#'
#' @keywords internal
"_PACKAGE"
