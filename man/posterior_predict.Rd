\name{posterior_predict}
\alias{posterior_predict}
\alias{posterior_predict.pclm}
\alias{posterior_predict.bpclm}
\alias{print.pclm_posterior_predict}
\alias{summary.pclm_posterior_predict}
\alias{plot.pclm_posterior_predict}
\title{Posterior or parametric predictive distribution of fine-cell counts}
\usage{
posterior_predict(fit, ...)

\method{posterior_predict}{bpclm}(fit,
                  type = c("predictive", "rate"),
                  level = 0.9, seed = NULL, ...)

\method{posterior_predict}{pclm}(fit,
                  type = c("predictive", "rate"),
                  level = 0.9, n_draws = 2000L, seed = NULL, ...)
}
\arguments{
\item{fit}{A fitted \code{"pclm"}, \code{"pclm_exact"} or
\code{"bpclm"} object.}
\item{type}{Either \code{"predictive"} (the default) or \code{"rate"}.
See Details.}
\item{level}{Credible/prediction level for the returned interval
(default 0.9).}
\item{n_draws}{For frequentist input objects, the number of
parametric multinomial draws to generate (default 2000).  Ignored
for \code{bpclm} input.}
\item{seed}{Optional integer for reproducibility.}
\item{...}{Currently unused.}
}
\value{
An object of class \code{"pclm_posterior_predict"} -- a list with
components \code{draws}, \code{mean}, \code{median}, \code{lower},
\code{upper}, \code{level}, \code{type}, \code{grid},
\code{grid_mid}, \code{wide_breaks}, \code{m}.
}
\description{
Generates draws from the posterior or parametric predictive
distribution of the latent fine-cell (typically single-year) counts.
With \code{type = "predictive"} (default), within each wide bin
\eqn{j} the conditional distribution of the fine-cell counts given
the band total \eqn{m_j} is multinomial of size \eqn{m_j} with
within-band probabilities \eqn{\pi_y / \gamma_j}.  With
\code{type = "rate"} the function returns draws of the latent
smooth rate \eqn{m_+ \pi}.
}
\details{
Two distinct sources of uncertainty exist:
\itemize{
  \item \emph{Rate uncertainty}: epistemic uncertainty about the
    smooth underlying rate \eqn{m_+ \pi}.  Captured by the
    \code{bpclm} posterior on \eqn{\pi}; absent from \code{pclm} /
    \code{pclm_exact} point estimates.
  \item \emph{Predictive uncertainty}: the multinomial sampling
    variability of the realised fine-cell counts within each wide
    bin, conditional on the band total.
}
For ungrouping use cases (band counts treated as data, fine-cell
counts as the unknown), \code{type = "predictive"} is usually the
right object, since it gives credible intervals that cover the
realised counts at the nominal level.

The wide-bin totals are reproduced exactly on every draw under
\code{type = "predictive"} (multinomial sampling within each band).
Under \code{type = "rate"} they are reproduced exactly only when
the input \code{fit} has been calibrated (see
\code{\link{calibrate}}) or is a \code{\link{pclm_exact}} fit.
}
\examples{
\donttest{
data(bloodlead)
fit <- bpclm(m = bloodlead$count,
             wide_breaks = with(bloodlead, cbind(lower, upper)),
             a = 0, b = 80, ngrid = 80, ndx = 17,
             niter = 2000, burnin = 500, adapt = 300, seed = 1)
fit <- calibrate(fit)
pp  <- posterior_predict(fit, type = "predictive")
plot(pp)
print(pp)
}
}
\seealso{
\code{\link{bpclm}}, \code{\link{calibrate}},
\code{\link{pclm_exact}}.
}
