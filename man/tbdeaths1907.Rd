\name{tbdeaths1907}
\alias{tbdeaths1907}
\docType{data}
\title{Tuberculosis deaths by age, The Netherlands, 1907 (illustrative)}
\description{
Counts of deaths attributed to tuberculosis in The Netherlands in
1907 by wide age bands. Used as Example 6.2 by Lambert and Eilers
(2009) to illustrate density estimation from heavily grouped
mortality data. Total deaths in the dataset: 9440.

\strong{Note:} the original CBS records used in the paper were
unpublished. The values shipped here are an \emph{illustrative
reconstruction} that preserves the published total and a plausible
early-20th-century age-mortality shape, intended for code testing
and demonstration of the methodology, not for historical inference.
See \code{data-raw/make-data.R} in the package source for
substitution instructions if you have access to real data.
}
\usage{data(tbdeaths1907)}
\format{
A data frame with 13 rows (12 age bands plus an upper-tail zero) and
3 columns: \code{lower}, \code{upper} (age band, years) and
\code{count} (number of deaths).
}
\references{
Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
from grouped continuous data. \emph{CSDA}, 53(4), 1388--1399.
}
\examples{
data(tbdeaths1907)
tbdeaths1907
}
\keyword{datasets}
