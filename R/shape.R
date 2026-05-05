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

#' Check whether a discrete density is unimodal
#'
#' Returns \code{TRUE} when the vector \code{pi} (a discrete latent
#' density on a fine grid) is unimodal, i.e. it is non-decreasing up to
#' some index and non-increasing thereafter (a flat plateau is
#' permitted). A small numerical tolerance is allowed.
#'
#' @param pi Numeric vector of fine-grid probabilities (length \eqn{I}).
#' @param tol Tolerance for "non-decrease" / "non-increase" comparisons
#'   (default \code{1e-12}). Allows floating-point wobble.
#'
#' @return Logical scalar.
#'
#' @examples
#' is_unimodal(c(0.1, 0.2, 0.4, 0.2, 0.1))   # TRUE
#' is_unimodal(c(0.1, 0.4, 0.2, 0.4, 0.1))   # FALSE (bimodal)
#'
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

#' Check whether a discrete density is log-concave
#'
#' A discrete density \code{pi} is log-concave when
#' \eqn{\log \pi_{i-1} + \log \pi_{i+1} \le 2 \log \pi_i} for all
#' interior indices \eqn{i}. Bins with zero probability are allowed
#' provided they are not interior to a strictly positive segment (which
#' would violate log-concavity by sending the second log-difference to
#' \eqn{-\infty} on one side and \eqn{+\infty} on the other).
#'
#' @param pi Numeric vector of fine-grid probabilities.
#' @param tol Numerical tolerance.
#'
#' @return Logical scalar.
#'
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

#' Check whether a discrete density is monotonic
#'
#' Returns \code{TRUE} if \code{pi} is non-decreasing or non-increasing
#' (within \code{tol}). Useful to enforce, e.g., a strictly decreasing
#' density on a positive support.
#'
#' @param pi Numeric vector.
#' @param direction Either \code{"either"} (default), \code{"increasing"}
#'   or \code{"decreasing"}.
#' @param tol Numerical tolerance.
#'
#' @return Logical scalar.
#'
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
