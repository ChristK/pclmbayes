# =============================================================================
# B-spline basis, difference penalty, and bin matrix
# -----------------------------------------------------------------------------
# These low-level helpers build the design objects used by both pclm() and
# bpclm(). The conventions follow Eilers & Marx (1996) and Lambert & Eilers
# (2009).
# =============================================================================

#' Build an equally-spaced B-spline basis on a fine grid
#'
#' Constructs a B-spline basis \eqn{B} of dimension \eqn{I \times K}
#' evaluated at the midpoints \eqn{u_i} of the fine-grid intervals
#' partitioning \eqn{(a, b)}. The knots are equally spaced and extend
#' beyond the support so that the basis is well-defined at the
#' boundaries (the standard P-spline construction of Eilers and Marx,
#' 1996). The number of basis functions is \eqn{K = ndx + degree}, where
#' \code{ndx} is the number of knot intervals inside \eqn{(a, b)} and
#' \code{degree} is the spline degree (cubic by default).
#'
#' @param x Numeric vector of evaluation points (typically the fine-grid
#'   midpoints).
#' @param a,b Lower and upper limits of the support \eqn{(a, b)}.
#' @param ndx Number of equally-spaced knot intervals on \eqn{(a, b)}.
#'   The number of B-splines is \code{K = ndx + degree}.
#' @param degree Spline degree (default 3 for cubic B-splines).
#'
#' @return An object of class \code{"pclm_basis"}, a list with components
#'   \describe{
#'     \item{B}{Numeric matrix of size \code{length(x) x K} with the
#'       basis evaluated at \code{x}.}
#'     \item{knots}{The augmented knot sequence used internally.}
#'     \item{a, b, ndx, degree, K}{The arguments and derived dimension.}
#'   }
#'
#' @references
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
#' B-splines and penalties. \emph{Statistical Science}, 11(2), 89--121.
#'
#' @examples
#' grid <- seq(0, 10, length.out = 100)
#' basis <- bspline_basis(grid, a = 0, b = 10, ndx = 17, degree = 3)
#' dim(basis$B)  # 100 x 20
#'
#' @export
bspline_basis <- function(x, a, b, ndx = 17L, degree = 3L) {

  # nocov start
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop("`x` must be a numeric vector with finite values.")
  }
  if (!is.numeric(a) || !is.numeric(b) || a >= b) {
    stop("Require `a < b`.")
  }
  if (any(x < a - 1e-8) || any(x > b + 1e-8)) {
    stop("All points in `x` must lie inside [a, b].")
  }
  # nocov end
  ndx    <- as.integer(ndx)
  degree <- as.integer(degree)
  # nocov start
  if (ndx < 1L)    stop("`ndx` must be at least 1.")
  if (degree < 0L) stop("`degree` must be non-negative.")
  # nocov end

  # equally-spaced augmented knot sequence
  h <- (b - a) / ndx
  knots <- seq(from = a - degree * h,
               to   = b + degree * h,
               by   = h)

  # number of B-splines
  K <- ndx + degree

  # Evaluate basis using base-R splines::splineDesign (no extra dependency)
  B <- splines::splineDesign(knots = knots,
                             x      = x,
                             ord    = degree + 1L,
                             outer.ok = TRUE)
  # B has length(x) rows and K columns
  # nocov start
  if (ncol(B) != K) {
    # Should not happen with a correct knot sequence, but guard anyway
    stop("Internal error: basis matrix has ", ncol(B),
         " columns, expected ", K, ".")
  }
  # nocov end

  structure(list(B      = B,
                 knots  = knots,
                 a      = a,
                 b      = b,
                 ndx    = ndx,
                 degree = degree,
                 K      = K),
            class = "pclm_basis")
}

#' Difference penalty matrix
#'
#' Returns the matrix \eqn{P = D'D} that gives the sum of squared
#' \eqn{r}th-order differences of a vector of length \eqn{K}, i.e.
#' \eqn{\sum_k (\Delta^r \phi_k)^2 = \phi' P \phi}. The matrix \eqn{D}
#' itself can be returned via \code{return_D = TRUE}.
#'
#' The rank of \eqn{P} is \eqn{K - r}; its null space is spanned by
#' polynomials of degree up to \eqn{r-1} (in the sequence index). This
#' rank-deficiency is the reason for the identifiability constraint
#' \eqn{\sum_k \phi_k = 0} in the model.
#'
#' @param K Length of the coefficient vector (number of B-splines).
#' @param r Penalty order (default 3, as recommended by Lambert and
#'   Eilers 2009 because it tends to a normal density at the limit of
#'   strong smoothing).
#' @param return_D Logical. If \code{TRUE}, the function returns a list
#'   with components \code{D} (the \eqn{(K-r) \times K} difference
#'   matrix) and \code{P} (the \eqn{K \times K} penalty matrix). If
#'   \code{FALSE} (the default), only \code{P} is returned as a matrix.
#'
#' @return Either a matrix \code{P} or a list (see \code{return_D}).
#'
#' @examples
#' diff_penalty(K = 6, r = 2)
#'
#' @export
diff_penalty <- function(K, r = 3L, return_D = FALSE) {
  K <- as.integer(K)
  r <- as.integer(r)
  # nocov start
  if (K < 2L)        stop("`K` must be at least 2.")
  if (r < 1L)        stop("`r` must be at least 1.")
  if (r >= K)        stop("`r` must be smaller than `K`.")
  # nocov end

  D <- diff(diag(K), differences = r)  # (K - r) x K
  P <- crossprod(D)                    # K x K
  if (return_D) list(D = D, P = P) else P
}

#' Build the bin (composite-link) matrix C
#'
#' Given the wide bin boundaries \eqn{(L_j, U_j)} (\eqn{j = 1,\dots,J})
#' and the fine-grid breakpoints, this function constructs the
#' \eqn{J \times I} composite-link matrix \eqn{C} such that the
#' wide-bin probabilities are \eqn{\gamma = C \pi}, where \eqn{\pi} is
#' the vector of fine-grid probabilities.
#'
#' Element \eqn{c_{ji}} is the proportion of the \eqn{i}th fine
#' interval that falls inside the \eqn{j}th wide bin. When wide-bin
#' boundaries align with the fine grid, this gives 0/1 entries (as in
#' Section 2 of Lambert and Eilers, 2009). Misaligned or overlapping
#' wide bins are handled by the rectangle method (Section 4 of the same
#' paper).
#'
#' @param wide_breaks Numeric matrix or data frame with two columns
#'   giving the lower and upper limit of each wide bin, or a sorted
#'   numeric vector of length \eqn{J + 1} when the wide bins partition
#'   \eqn{(a, b)} contiguously.
#' @param fine_breaks Numeric vector of length \eqn{I + 1} giving the
#'   fine-grid breakpoints.
#'
#' @return A \eqn{J \times I} numeric matrix \code{C}.
#'
#' @examples
#' fine_breaks <- seq(0, 10, by = 0.5)
#' # five wide bins of width 2
#' wide_breaks <- seq(0, 10, by = 2)
#' C <- bin_matrix(wide_breaks, fine_breaks)
#' dim(C)            # 5 x 20
#' rowSums(C)        # all 1: every wide bin fully covered
#'
#' @export
bin_matrix <- function(wide_breaks, fine_breaks) {

  fine_breaks <- as.numeric(fine_breaks)
  # nocov start
  if (any(diff(fine_breaks) <= 0)) {
    stop("`fine_breaks` must be strictly increasing.")
  }
  # nocov end
  I <- length(fine_breaks) - 1L

  # Convert wide_breaks to a J x 2 matrix [L_j, U_j]
  if (is.matrix(wide_breaks) || is.data.frame(wide_breaks)) {
    wb <- as.matrix(wide_breaks)
    # nocov start
    if (ncol(wb) != 2L) {
      stop("`wide_breaks` matrix/data frame must have two columns ",
           "(lower, upper limits).")
    }
    # nocov end
    L <- wb[, 1L]
    U <- wb[, 2L]
  } else {
    wb <- sort(as.numeric(wide_breaks))
    # nocov start
    if (length(wb) < 2L) {
      stop("`wide_breaks` must contain at least two breakpoints.")
    }
    # nocov end
    L <- head(wb, -1L)
    U <- tail(wb, -1L)
  }
  # nocov start
  if (any(U <= L)) {
    stop("Each wide bin must have upper > lower limit.")
  }
  # nocov end
  J <- length(L)

  # Pre-compute fine interval bounds and widths
  fL <- head(fine_breaks, -1L)            # lower bounds of fine intervals
  fU <- tail(fine_breaks, -1L)            # upper bounds
  fW <- fU - fL                           # widths

  C <- matrix(0, nrow = J, ncol = I)
  for (j in seq_len(J)) {
    lo <- pmax(L[j], fL)
    hi <- pmin(U[j], fU)
    overlap <- pmax(0, hi - lo)
    # weight is the proportion of the fine interval that lies inside [L_j, U_j]
    C[j, ] <- overlap / fW
  }
  C
}

# ---- Internal helpers -------------------------------------------------------

# Numerically stable softmax: pi_i = exp(eta_i) / sum_l exp(eta_l).
# Always returns a numeric vector (drops any matrix dim that comes from
# `B %*% phi` etc.), so downstream element-wise broadcasts (e.g. `pi * B`)
# behave correctly.
.softmax <- function(eta) {
  eta <- as.numeric(eta)
  m  <- max(eta)
  ex <- exp(eta - m)
  ex / sum(ex)
}

# Centre a vector to sum to zero (used to enforce the identifiability
# constraint \sum_k phi_k = 0)
.centre <- function(phi) phi - mean(phi)

# Latent grid probabilities pi(phi) given a basis B
.pi_phi <- function(phi, B) .softmax(B %*% phi)

# Wide-bin probabilities gamma(phi)
.gamma_phi <- function(phi, B, C) as.numeric(C %*% .pi_phi(phi, B))

#' Evaluate the latent density on the fine grid
#'
#' Returns the discrete latent density values
#' \eqn{\pi_i = f_Y(u_i) \Delta} at the fine-grid midpoints
#' \eqn{u_i}, given a vector of B-spline coefficients.
#'
#' @param phi Numeric vector of B-spline coefficients of length \eqn{K}.
#' @param basis A \code{"pclm_basis"} object, as returned by
#'   \code{\link{bspline_basis}}.
#'
#' @return Numeric vector of length \eqn{I} (number of fine grid
#'   intervals), summing to 1.
#'
#' @examples
#' grid <- seq(0, 10, length.out = 51)
#' mids <- (head(grid, -1) + tail(grid, -1)) / 2
#' basis <- bspline_basis(mids, a = 0, b = 10, ndx = 17)
#' phi <- rep(0, basis$K)
#' pi_grid <- latent_density(phi, basis)
#' sum(pi_grid)  # 1
#'
#' @export
latent_density <- function(phi, basis) {
  if (!inherits(basis, "pclm_basis")) {
    stop("`basis` must be a 'pclm_basis' object (see bspline_basis()).")
  }
  if (length(phi) != basis$K) {
    stop("`phi` must have length basis$K = ", basis$K, ".")
  }
  if (any(!is.finite(phi))) {
    stop("`phi` contains non-finite values.")
  }
  .pi_phi(phi, basis$B)
}

# -----------------------------------------------------------------------------
# Multinomial composite-link likelihood and its gradient
# -----------------------------------------------------------------------------
#
# Let pi_i(phi) = exp(eta_i) / sum_l exp(eta_l) with eta = B phi,
#     gamma = C pi.
# Log-likelihood:   ell(phi) = sum_j m_j log gamma_j
# Gradient w.r.t. phi:
#   d pi / d phi  = diag(pi) %*% B - pi %*% (t(pi) %*% B)
#   d gamma / d phi = C %*% (d pi / d phi)
#   d ell / d phi   = (m / gamma)' %*% (d gamma / d phi)
#
# Equivalently (more efficient): with omega_j = m_j / gamma_j,
# alpha = t(C) %*% omega, then
#   d ell / d phi = t(B) %*% (alpha * pi) - (sum(alpha * pi)) * (t(B) %*% pi)
#
# Fisher information (negative expected Hessian):
#   I(phi) = m_+ * t(d gamma / d phi) %*% diag(1 / gamma) %*% (d gamma / d phi)
# where m_+ = sum_j m_j.
#
# The implementation below returns log-likelihood, gradient, and the Fisher
# information in one place to amortise shared work.
# -----------------------------------------------------------------------------

# Internal: shared computation of pi, gamma, log-likelihood, gradient, and
# Fisher information given phi, design matrices, and observed counts m.
#
# Returns a list with:
#   pi     : length I
#   gamma  : length J
#   logL   : numeric (uses 0 * log(0) = 0 convention)
#   grad   : length K
#   FI     : K x K Fisher information (only if `compute_FI = TRUE`)
#
# Note on R broadcasting: `pi * B` (pi length I, B is I x K) recycles pi
# in column-major order. Because B has I rows, this is equivalent to
# diag(pi) %*% B (each ROW of B is multiplied by the corresponding pi[i]).
# Similarly, `inv_g * dgamma` (inv_g length J, dgamma is J x K) is
# equivalent to diag(inv_g) %*% dgamma.
.pclm_eval <- function(phi, B, C, m, compute_FI = TRUE) {

  m_plus <- sum(m)

  pi  <- .softmax(B %*% phi)             # length I
  ga  <- as.numeric(C %*% pi)            # length J

  # Log-likelihood (multinomial up to a constant in m)
  pos  <- m > 0
  if (any(ga[pos] <= 0)) {
    logL <- -Inf
  } else {
    logL <- sum(m[pos] * log(ga[pos]))
  }

  # Gradient: t(B) %*% (alpha * pi) - (alpha %*% pi) * (t(B) %*% pi)
  # where alpha = t(C) %*% (m / gamma). Use safe division.
  omega <- numeric(length(m))
  ok    <- ga > 0
  omega[ok] <- m[ok] / ga[ok]
  alpha     <- as.numeric(crossprod(C, omega))  # length I
  ap        <- sum(alpha * pi)
  grad      <- as.numeric(crossprod(B, alpha * pi) - ap * crossprod(B, pi))

  out <- list(pi = pi, gamma = ga, logL = logL, grad = grad)

  if (compute_FI) {
    # d gamma / d phi  =  C %*% (diag(pi) %*% B - pi %*% (pi' B))
    pB     <- as.numeric(crossprod(B, pi))           # length K, t(B) %*% pi
    dpidphi <- (pi * B) - tcrossprod(pi, pB)         # I x K
    dgamma  <- C %*% dpidphi                          # J x K
    # Inverse-variance weights: 1 / gamma_j (with safe handling of zeros)
    inv_g <- numeric(length(ga))
    inv_g[ok] <- 1 / ga[ok]
    # Fisher info = m_+ * t(dgamma) %*% diag(inv_g) %*% dgamma
    out$FI <- m_plus * crossprod(dgamma, inv_g * dgamma)
  }

  out
}
