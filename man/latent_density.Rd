\name{latent_density}
\alias{latent_density}
\title{Latent density on the fine grid}
\usage{
latent_density(phi, basis)
}
\arguments{
\item{phi}{Numeric vector of B-spline coefficients (length
\code{basis$K}).}
\item{basis}{A \code{"pclm_basis"} object from \code{\link{bspline_basis}}.}
}
\value{
Numeric vector of length \eqn{I} (number of fine grid intervals),
summing to 1.
}
\description{
Returns the discrete latent density values
\eqn{\pi_i = f_Y(u_i) \Delta} at the fine-grid midpoints given a
vector of B-spline coefficients.
}
