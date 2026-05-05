\name{calibrate}
\alias{calibrate}
\alias{calibrate.pclm}
\alias{calibrate.bpclm}
\title{Calibrate a fitted PCLM to preserve wide-bin totals exactly}
\usage{
calibrate(fit, ...)

\method{calibrate}{pclm}(fit, ...)

\method{calibrate}{bpclm}(fit, ...)
}
\arguments{
\item{fit}{A fitted \code{"pclm"} or \code{"bpclm"} object.}
\item{...}{Currently unused.}
}
\value{
The same object with calibrated \code{pi} (and, for \code{bpclm},
calibrated \code{pi_chain}, \code{pi_lower}, \code{pi_upper}). The
component \code{fitted_counts} is set to equal \code{m} exactly, and
\code{calibrated = TRUE} is recorded.
}
\description{
Post-hoc renormalisation of a fitted \code{pclm} or \code{bpclm}
object so that the implied wide-bin counts exactly equal the
observed counts \code{m}. For each wide bin \eqn{j}, the fine-grid
probabilities \eqn{\pi_i} (\eqn{i} in bin \eqn{j}) are scaled by
\eqn{(m_j / m_+) / \gamma_j}.

Calibration assumes the wide bins partition the support without
overlap so that each fine interval belongs to exactly one wide bin.
For overlapping or misaligned bins, use \code{\link{pclm_exact}}
instead.
}
\examples{
data(bloodlead)
fit   <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
              penalty_order = 3)
fit_c <- calibrate(fit)
max(abs(fit_c$fitted_counts - fit_c$m))
}
\seealso{
\code{\link{pclm_exact}}.
}
