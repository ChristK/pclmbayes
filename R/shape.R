# =============================================================================
# Shape-constraint indicators (Section 3.3, Eq. 7 of Lambert & Eilers 2009)
# -----------------------------------------------------------------------------
# These functions return TRUE if the latent density on the fine grid satisfies
# the requested shape property, FALSE otherwise. They are used in the Bayesian
# sampler to reject proposals that fall outside the indicator set Phi
# associated with the desired prior.
#
# The constraints operate on the latent grid probabilities pi, not directly on
# the spline coefficients phi. This means a proposal is accepted only if its
# induced fine-grid density has the requested property; the constraint is in
# that sense "soft" with respect to phi.
# =============================================================================

#' Shape predicates for discrete densities
#'
#' Predicates that check whether a discrete latent density \code{pi}
#' (typically a fine-grid probability vector from a PCLM fit) satisfies
#' a particular shape constraint:
#'
#' \describe{
#'   \item{\code{is_unimodal}}{Non-decreasing up to some index and
#'     non-increasing thereafter (a flat plateau is permitted).}
#'   \item{\code{is_logconcave}}{
#'     \eqn{\log \pi_{i-1} + \log \pi_{i+1} \le 2 \log \pi_i} at every
#'     interior index. Bins with zero probability are allowed provided
#'     they are not interior to a strictly positive segment.}
#'   \item{\code{is_monotonic}}{Non-decreasing or non-increasing across
#'     the whole support, possibly restricted to one direction.}
#' }
#'
#' Each predicate uses a small numerical tolerance to absorb
#' floating-point noise. They are used internally by \code{\link{bpclm}}
#' to enforce shape priors via rejection.
#'
#' @param pi Numeric vector of fine-grid probabilities (length \eqn{I}).
#' @param tol Numerical tolerance for the inequality comparisons.
#' @param direction For \code{is_monotonic}: either \code{"either"} (the
#'   default), \code{"increasing"}, or \code{"decreasing"}.
#'
#' @return Logical scalar.
#'
#' @examples
#' is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1))   # TRUE
#' is_unimodal(c(0.1, 0.4, 0.2, 0.4, 0.1))   # FALSE (bimodal)
#' is_logconcave(dnorm(seq(-3, 3, length.out = 31)))            # TRUE
#' is_monotonic(c(0.5, 0.3, 0.1), direction = "decreasing")     # TRUE
#'
#' @name shape-constraints
#' @rdname shape-constraints
#' @export
is_unimodal <- function(pi, tol = 1e-12) {
  pi <- as.numeric(pi)
  if (length(pi) < 3L) return(TRUE)
  d <- diff(pi)
  # Find sign changes ignoring tiny ones
  s <- ifelse(d >  tol,  1L,
       ifelse(d < -tol, -1L, 0L))
  s <- s[s != 0L]
  if (length(s) <= 1L) return(TRUE)
  # Allowed pattern: zero or more +1 followed by zero or more -1
  n_changes <- sum(diff(s) != 0L)
  if (n_changes == 0L)              return(TRUE)
  if (n_changes == 1L && s[1L] == 1L) return(TRUE)
  FALSE
}

#' @rdname shape-constraints
#' @export
is_logconcave <- function(pi, tol = 1e-10) {
  pi <- as.numeric(pi)
  if (length(pi) < 3L) return(TRUE)
  if (any(pi < 0))     return(FALSE)
  # Replace zeros with a tiny number for log-difference; alternatively, log on
  # zero gives -Inf and the inequality may still hold or fail consistently.
  lp <- log(pmax(pi, .Machine$double.xmin))
  d2 <- diff(lp, differences = 2L)
  all(d2 <= tol)
}

#' @rdname shape-constraints
#' @export
is_monotonic <- function(pi,
                         direction = c("either", "increasing", "decreasing"),
                         tol = 1e-12) {
  direction <- match.arg(direction)
  pi <- as.numeric(pi)
  if (length(pi) < 2L) return(TRUE)
  d <- diff(pi)
  if (direction == "increasing")  return(all(d >= -tol))
  if (direction == "decreasing")  return(all(d <=  tol))
  all(d >= -tol) || all(d <= tol)
}

# Internal: combine multiple shape constraints into a single indicator.
# `shape` is either NULL (no constraint), a single string, or a character
# vector of constraint names. The "monotonic" entry can carry a direction via
# the named list element `direction` (passed through bpclm via `shape_args`).
.shape_ok <- function(pi, shape, shape_args = list()) {
  if (is.null(shape) || length(shape) == 0L) return(TRUE)
  shape <- tolower(shape)
  for (s in shape) {
    ok <- switch(s,
      "unimodal"    = is_unimodal(pi),
      "logconcave"  = is_logconcave(pi),
      "log-concave" = is_logconcave(pi),
      "monotonic"   = is_monotonic(
                        pi,
                        direction = if (!is.null(shape_args$direction))
                                       shape_args$direction
                                    else "either"),
      stop("Unknown shape constraint: '", s, "'.")
    )
    if (!isTRUE(ok)) return(FALSE)
  }
  TRUE
}
