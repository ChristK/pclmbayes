\name{bpclm}
\alias{bpclm}
\alias{print.bpclm}
\alias{summary.bpclm}
\alias{plot.bpclm}
\alias{coef.bpclm}
\alias{fitted.bpclm}
\alias{predict.bpclm}
\alias{quantile.bpclm}
\title{Bayesian penalised composite link model for grouped data}
\usage{
bpclm(m, wide_breaks,
      a = NULL, b = NULL,
      ngrid = 100L,
      ndx = 17L, degree = 3L,
      penalty_order = 3L,
      niter = 5000L,
      burnin = NULL,
      thin = 1L,
      adapt = 500L,
      tau_init = NULL,
      a_tau = 1e-4, b_tau = 1e-4,
      delta_init = NULL,
      target_accept = 0.57,
      Sigma = NULL,
      shape = NULL,
      shape_args = list(),
      phi_init = NULL,
      seed = NULL,
      verbose = FALSE)
}
\arguments{
\item{m, wide_breaks, a, b, ngrid, ndx, degree, penalty_order}{See
\code{\link{pclm}}.}
\item{niter}{Number of MCMC iterations (including burn-in). Default
5000.}
\item{burnin}{Number of initial iterations to discard. Default
\code{floor(niter / 5)}.}
\item{thin}{Thinning interval (default 1).}
\item{adapt}{Number of additional adaptive-tuning iterations
performed before the main run. Default 500.}
\item{tau_init}{Initial value of the precision \eqn{\tau}.
\code{NULL} (default) uses the BIC-selected value of the
frequentist warm-start fit.}
\item{a_tau, b_tau}{Hyperparameters of the \eqn{\Gamma(a, b)}
hyperprior on \eqn{\tau}. Defaults \code{1e-4} each.}
\item{delta_init}{Initial Langevin step size; default
\eqn{1.65^2 / K^{1/3}}.}
\item{target_accept}{Target acceptance rate during adaptation
(default 0.57).}
\item{Sigma}{Optional positive-definite \eqn{K \times K} proposal
covariance. \code{NULL} (default) uses the frequentist warm-start
variance-covariance.}
\item{shape}{Optional shape constraint(s). Any subset of
\code{"unimodal"}, \code{"logconcave"}, \code{"monotonic"}.}
\item{shape_args}{Additional arguments for the shape indicators
(currently only \code{direction} for monotonicity).}
\item{phi_init}{Optional starting value for \eqn{\phi}; defaults to
the frequentist MLE.}
\item{seed}{Optional integer for reproducibility.}
\item{verbose}{Print progress every 10\% of iterations.}
}
\value{
An object of class \code{"bpclm"}; see \code{\link{pclmbayes-package}}.
}
\description{
Fits the model of Lambert and Eilers (2009) by Markov chain Monte
Carlo. The sampler is Metropolis-within-Gibbs: a modified
Langevin-Hastings update for the spline coefficients followed by a
Gibbs draw for the smoothing precision \eqn{\tau}. Optional shape
constraints (\code{"unimodal"}, \code{"logconcave"}, \code{"monotonic"})
are imposed by rejecting proposals whose induced density violates
them.
}
\examples{
\donttest{
data(bloodlead)
fit <- bpclm(m = bloodlead$count,
             wide_breaks = with(bloodlead, cbind(lower, upper)),
             a = 0, b = 80, ngrid = 80, ndx = 17,
             niter = 2000, burnin = 500, shape = "unimodal", seed = 1)
summary(fit)
plot(fit)
}
}
\references{
Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. \emph{CSDA}, 53(4), 1388--1399.

Roberts, G. O. and Rosenthal, J. S. (1998). Optimal scaling of
discrete approximations to Langevin diffusions. \emph{JRSS B},
60(1), 255--268.
}
\seealso{
\code{\link{pclm}}.
}
