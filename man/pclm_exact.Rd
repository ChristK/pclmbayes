\name{pclm_exact}
\alias{pclm_exact}
\title{Smoothest density that exactly preserves wide-bin counts}
\usage{
pclm_exact(m, wide_breaks,
           a = NULL, b = NULL,
           ngrid = 100L,
           ndx = 17L, degree = 3L,
           penalty_order = 3L,
           max_iter = 200L, tol = 1e-9,
           eps_ridge = 1e-8,
           phi_start = NULL,
           verbose = FALSE)
}
\arguments{
\item{m, wide_breaks, a, b, ngrid, ndx, degree, penalty_order}{See
\code{\link{pclm}}.}
\item{max_iter, tol, eps_ridge}{Convergence and numerical stability
controls for the SQP iteration.}
\item{phi_start}{Optional starting coefficient vector.  Default uses
a mildly-smoothed \code{\link{pclm}} fit.}
\item{verbose}{Print one line per SQP iteration.}
}
\value{
An object of class \code{c("pclm_exact", "pclm")}, with components
as in \code{\link{pclm}}, plus \code{lambda} (the converged Lagrange
multipliers) and \code{constraint_residual} (max
\eqn{|\gamma - m / m_+|}).
}
\description{
Solves the equality-constrained problem
\deqn{\min_\phi \tfrac{1}{2} \phi' P \phi
  \quad \text{s.t.} \quad
  C \pi(\phi) = m / m_+,}
returning the smoothest density on the fine grid that exactly
reproduces the wide-bin counts when summed.  Solved by sequential
quadratic programming with the constraint linearised at each step.
}
\examples{
data(bloodlead)
fit <- pclm_exact(m = bloodlead$count,
                  wide_breaks = with(bloodlead, cbind(lower, upper)),
                  a = 0, b = 80, ngrid = 80, ndx = 17,
                  degree = 3, penalty_order = 3)
fit$constraint_residual
max(abs(fit$fitted_counts - fit$m))
}
\seealso{
\code{\link{pclm}}, \code{\link{calibrate}}.
}
