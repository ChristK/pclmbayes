\name{diff_penalty}
\alias{diff_penalty}
\title{Difference penalty matrix}
\usage{
diff_penalty(K, r = 3L, return_D = FALSE)
}
\arguments{
\item{K}{Length of the coefficient vector.}
\item{r}{Penalty order (default 3).}
\item{return_D}{If \code{TRUE}, return both \code{D} (the
\eqn{(K-r) \times K} difference matrix) and \code{P} (the
\eqn{K \times K} penalty matrix).}
}
\value{
Either a matrix \code{P} or a list with components \code{D} and
\code{P} (see \code{return_D}).
}
\description{
Returns the matrix \eqn{P = D'D} that gives the sum of squared
\eqn{r}th-order differences of a vector of length \code{K}, i.e.
\deqn{\sum_k (\Delta^r \phi_k)^2 = \phi' P \phi.}
The rank of \eqn{P} is \eqn{K - r}; its null space is spanned by
polynomials in the index of degree up to \eqn{r-1}.
}
\examples{
diff_penalty(K = 6, r = 2)
}
