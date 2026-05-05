\name{is_unimodal}
\alias{is_unimodal}
\alias{is_logconcave}
\alias{is_monotonic}
\title{Shape-constraint indicators}
\usage{
is_unimodal(pi, tol = 1e-12)
is_logconcave(pi, tol = 1e-10)
is_monotonic(pi,
             direction = c("either", "increasing", "decreasing"),
             tol = 1e-12)
}
\arguments{
\item{pi}{Numeric vector of fine-grid probabilities (a discrete
density).}
\item{tol}{Numerical tolerance.}
\item{direction}{For monotonicity: \code{"either"}, \code{"increasing"}
or \code{"decreasing"}.}
}
\value{
Logical scalar.
}
\description{
Returns \code{TRUE} if \code{pi} satisfies the requested shape
property. Used by \code{\link{bpclm}} to enforce the rejection prior
in Eq. (7) of Lambert and Eilers (2009).
}
\examples{
is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1))
is_logconcave(dnorm(seq(-3, 3, length.out = 31)) /
              sum(dnorm(seq(-3, 3, length.out = 31))))
is_monotonic(c(0.5, 0.3, 0.1), direction = "decreasing")
}
