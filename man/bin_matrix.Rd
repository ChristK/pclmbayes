\name{bin_matrix}
\alias{bin_matrix}
\title{Build the bin (composite-link) matrix C}
\usage{
bin_matrix(wide_breaks, fine_breaks)
}
\arguments{
\item{wide_breaks}{Numeric matrix or data frame with two columns giving
the lower and upper limits of each wide bin, or a sorted numeric
vector of length \eqn{J + 1} when the wide bins partition
\eqn{(a, b)} contiguously.}
\item{fine_breaks}{Numeric vector of length \eqn{I + 1} giving the
fine-grid breakpoints.}
}
\value{
A \eqn{J \times I} matrix \code{C}, where \code{C[j, i]} is the
proportion of fine interval \eqn{i} that lies in wide bin \eqn{j}.
}
\description{
Builds the composite-link matrix \eqn{C} such that the wide-bin
probabilities are \eqn{\gamma = C \pi}, where \eqn{\pi} is the
vector of fine-grid probabilities. Handles both aligned and
misaligned wide bins via the rectangle method.
}
\examples{
fine  <- seq(0, 10, by = 0.5)
wide  <- seq(0, 10, by = 2)
C     <- bin_matrix(wide, fine)
dim(C)
}
