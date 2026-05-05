\name{bspline_basis}
\alias{bspline_basis}
\title{Equally-spaced B-spline basis on a fine grid}
\usage{
bspline_basis(x, a, b, ndx = 17L, degree = 3L)
}
\arguments{
\item{x}{Numeric vector of evaluation points.}
\item{a, b}{Lower and upper limits of the support \eqn{(a, b)}.}
\item{ndx}{Number of equally-spaced knot intervals on \eqn{(a, b)}.}
\item{degree}{Spline degree (default 3 for cubic B-splines).}
}
\value{
A list of class \code{"pclm_basis"} with components \code{B}
(the basis matrix), \code{knots}, \code{a}, \code{b}, \code{ndx},
\code{degree}, \code{K}.
}
\description{
Builds a B-spline basis evaluated at \code{x}, with equally-spaced
knots extending beyond the support to give partition-of-unity at all
interior points (the standard P-spline construction of Eilers and
Marx, 1996). The number of basis functions is
\code{K = ndx + degree}.
}
\examples{
mids <- seq(0.05, 9.95, by = 0.1)
bs   <- bspline_basis(mids, a = 0, b = 10, ndx = 17, degree = 3)
dim(bs$B)
}
\references{
Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. \emph{Statistical Science}, 11(2), 89--121.
}
