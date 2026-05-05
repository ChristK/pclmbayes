\name{pclm}
\alias{pclm}
\alias{print.pclm}
\alias{summary.pclm}
\alias{plot.pclm}
\alias{coef.pclm}
\alias{fitted.pclm}
\alias{predict.pclm}
\alias{quantile.pclm}
\alias{logLik.pclm}
\title{Frequentist penalised composite link model for grouped data}
\usage{
pclm(m, wide_breaks,
     a = NULL, b = NULL,
     ngrid = 100L,
     ndx = 17L, degree = 3L,
     penalty_order = 3L,
     tau = NULL,
     select = c("BIC", "AIC"),
     max_iter = 100L,
     tol = 1e-7,
     phi_start = NULL,
     verbose = FALSE)
}
\arguments{
\item{m}{Numeric vector of non-negative wide-bin counts.}
\item{wide_breaks}{Either a sorted vector of \eqn{J + 1} contiguous
wide-bin breakpoints or a \eqn{J \times 2} matrix/data frame of
(lower, upper) limits.}
\item{a, b}{Lower and upper limits of the support \eqn{(a, b)}.
Defaults to the smallest and largest wide-bin limits.}
\item{ngrid}{Number of fine-grid intervals \eqn{I}. Default 100.}
\item{ndx}{Number of equally-spaced knot intervals on \eqn{(a, b)}.
The number of B-splines is \code{ndx + degree}. Default 17.}
\item{degree}{B-spline degree. Default 3.}
\item{penalty_order}{Order of the difference penalty. Default 3.}
\item{tau}{Smoothing parameter; either a scalar or a vector of
candidate values. \code{NULL} (default) uses an internal log-spaced
grid of length 25.}
\item{select}{Information criterion: \code{"BIC"} (default) or
\code{"AIC"}.}
\item{max_iter, tol}{Convergence controls for the penalised scoring
algorithm.}
\item{phi_start}{Optional starting coefficient vector.}
\item{verbose}{Logical; print iteration progress.}
}
\value{
An object of class \code{"pclm"}; see \code{\link{pclmbayes-package}}.
}
\description{
Fits the model of Lambert and Eilers (2009) by penalised maximum
likelihood. The smoothing parameter \eqn{\tau} is selected by the
chosen information criterion when a grid of candidates is supplied.
}
\examples{
data(bloodlead)
fit <- pclm(m = bloodlead$count,
            wide_breaks = with(bloodlead, cbind(lower, upper)),
            a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
            penalty_order = 3)
summary(fit)
plot(fit)
}
\references{
Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. \emph{CSDA}, 53(4), 1388--1399.
}
\seealso{
\code{\link{bpclm}}.
}
