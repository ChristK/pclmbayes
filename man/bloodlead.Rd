\name{bloodlead}
\alias{bloodlead}
\docType{data}
\title{Lead concentration in the blood of New York children, 1974}
\description{
Interval-censored measurements of lead concentration (in
\eqn{\mu}g/dl) in the blood of young Puerto Ricans aged 1-12 years
living in New York in 1974, as reported by Hasselblad et al. (1980)
and re-used as Example 6.1 by Lambert and Eilers (2009). Total
\eqn{n = 139}.
}
\usage{data(bloodlead)}
\format{
A data frame with 7 rows and 3 columns: \code{lower}, \code{upper}
(bin limits in \eqn{\mu}g/dl) and \code{count} (frequency).
}
\source{
Hasselblad, V., Stead, A. G. and Galke, W. (1980). Analysis of
coarsely grouped data from the lognormal distribution.
\emph{JASA}, 75, 771--778.
}
\references{
Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. \emph{CSDA}, 53(4), 1388--1399.
}
\examples{
data(bloodlead)
bloodlead
}
\keyword{datasets}
