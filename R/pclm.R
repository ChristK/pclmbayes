# =============================================================================
# Frequentist penalised composite link model fit
# =============================================================================
#
# Implements the penalised scoring (penalised IRLS / penalised Newton-Raphson)
# algorithm of Lambert and Eilers (2009), Eq. (2). The same paper's Eq. (3)
# gives the variance-covariance matrix of the spline coefficients at
# convergence, which is also used as the proposal covariance Sigma in the
# Bayesian sampler (bpclm).
# =============================================================================

#' Frequentist penalised composite link model for grouped data
#'
#' Fits the model of Lambert and Eilers (2009) by penalised maximum
#' likelihood. The log of the latent density is modelled as a linear
#' combination of B-splines with an \eqn{r}th-order difference penalty
#' on the coefficients. Wide-bin counts \code{m} are assumed to follow a
#' multinomial distribution with probabilities
#' \eqn{\gamma = C \pi(\phi)}.
#'
#' For a given smoothing parameter \eqn{\tau}, the fit is obtained by
#' the penalised scoring algorithm (Eq. 2 of the paper). When a grid of
#' candidate \eqn{\tau} values is supplied via \code{tau} (a numeric
#' vector of length \eqn{> 1}), the value minimising the chosen
#' information criterion (\code{select = "BIC"}, the default, or
#' \code{"AIC"}) is selected. The variance-covariance matrix of the
#' coefficients is the inverse of the penalised Fisher information at
#' convergence (Eq. 3 of the paper).
#'
#' @param m Numeric vector of non-negative wide-bin counts (length
#'   \eqn{J}).
#' @param wide_breaks Wide-bin boundaries; either a \eqn{J + 1} vector
#'   (contiguous bins partitioning \eqn{(a, b)}) or a \eqn{J \times 2}
#'   matrix/data frame of (lower, upper) limits.
#' @param a,b Lower and upper limits of the support \eqn{(a, b)} on
#'   which to estimate the density. Defaults to the smallest lower limit
#'   and largest upper limit of \code{wide_breaks}.
#' @param ngrid Number of fine-grid intervals \eqn{I}. Defaults to 100.
#' @param ndx Number of equally-spaced knot intervals on \eqn{(a, b)}.
#'   The number of B-splines is \code{ndx + degree}. Default 17 (so
#'   \eqn{K = 20} cubic B-splines, matching the paper's examples).
#' @param degree B-spline degree (default 3, cubic).
#' @param penalty_order Order \eqn{r} of the difference penalty
#'   (default 3, as in the paper's examples).
#' @param tau Either a positive scalar smoothing parameter, or a numeric
#'   vector of candidate values to evaluate. If \code{NULL} (the
#'   default), a logarithmically-spaced grid from \code{1e-2} to
#'   \code{1e6} (length 25) is used.
#' @param select Information criterion used to pick \eqn{\tau} when a
#'   grid is provided: either \code{"BIC"} (default) or \code{"AIC"}.
#' @param max_iter Maximum number of scoring iterations.
#' @param tol Convergence tolerance on the largest absolute change in
#'   \eqn{\phi}.
#' @param phi_start Optional starting value for \eqn{\phi} (length
#'   \code{ndx + degree}). Defaults to the zero vector (uniform
#'   density).
#' @param verbose Logical: if \code{TRUE}, print one line per scoring
#'   iteration.
#'
#' @return An object of class \code{"pclm"}, a list with components:
#'   \describe{
#'     \item{phi}{Estimated B-spline coefficients (length \eqn{K},
#'       summing to 0).}
#'     \item{tau}{Selected (or supplied) smoothing parameter.}
#'     \item{tau_grid}{Numeric vector: the full grid of \eqn{\tau}
#'       values evaluated. \code{NULL} if a single \eqn{\tau} was
#'       passed.}
#'     \item{ic}{Information-criterion value at the selected
#'       \eqn{\tau}.}
#'     \item{ic_grid}{Vector of IC values at the candidate \eqn{\tau}
#'       values (\code{NULL} if a single \eqn{\tau}).}
#'     \item{select}{Which IC was used (\code{"BIC"} or \code{"AIC"}).}
#'     \item{logL}{Log-likelihood at convergence.}
#'     \item{edf}{Effective degrees of freedom,
#'       \eqn{\mathrm{tr}((I + \tau P)^{-1} I)}.}
#'     \item{vcov}{Estimated variance-covariance matrix of \eqn{\phi}
#'       (Eq. 3 of the paper).}
#'     \item{pi}{Latent grid probabilities \eqn{\pi_i}.}
#'     \item{gamma}{Fitted wide-bin probabilities \eqn{\gamma_j}.}
#'     \item{fitted_counts}{\eqn{m_+ \gamma}, the multinomial
#'       expectations.}
#'     \item{grid}{Fine-grid breakpoints (length \eqn{I + 1}).}
#'     \item{grid_mid}{Fine-grid midpoints (length \eqn{I}).}
#'     \item{basis}{The \code{"pclm_basis"} object used.}
#'     \item{C}{The bin matrix \eqn{C}.}
#'     \item{m, wide_breaks}{The supplied data.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @references
#' Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
#' from grouped continuous data. \emph{Computational Statistics and
#' Data Analysis}, 53(4), 1388--1399.
#'
#' @seealso \code{\link{bpclm}} for the Bayesian variant.
#'
#' @examples
#' data(bloodlead)
#' fit <- pclm(m = bloodlead$count,
#'             wide_breaks = with(bloodlead, cbind(lower, upper)),
#'             a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
#'             penalty_order = 3)
#' summary(fit)
#'
#' @export
pclm <- function(m, wide_breaks,
                 a = NULL, b = NULL,
                 ngrid = 100L,
                 ndx = 17L, degree = 3L,
                 penalty_order = 3L,
                 tau = NULL,
                 select = c("BIC", "AIC"),
                 max_iter = 100L,
                 tol = 1e-7,
                 phi_start = NULL,
                 verbose = FALSE) {

  cl <- match.call()
  select <- match.arg(select)

  # ---- Validate inputs ---------------------------------------------------
  if (!is.numeric(m) || any(m < 0) || any(!is.finite(m))) {
    stop("`m` must be a non-negative numeric vector.")
  }
  if (sum(m) <= 0) stop("`sum(m)` must be positive.")

  # Standardise wide_breaks to J x 2
  if (is.matrix(wide_breaks) || is.data.frame(wide_breaks)) {
    wb <- as.matrix(wide_breaks)
    if (ncol(wb) != 2L) {
      stop("If `wide_breaks` is a matrix, it must have two columns.")
    }
  } else {
    wb_vec <- sort(as.numeric(wide_breaks))
    wb <- cbind(head(wb_vec, -1), tail(wb_vec, -1))
  }
  if (nrow(wb) != length(m)) {
    stop("Number of wide bins (rows of `wide_breaks`) must equal length(m).")
  }

  if (is.null(a)) a <- min(wb[, 1L])
  if (is.null(b)) b <- max(wb[, 2L])
  if (a >= b) stop("Require `a < b`.")

  # Fine grid
  ngrid       <- as.integer(ngrid)
  fine_breaks <- seq(a, b, length.out = ngrid + 1L)
  mids        <- (head(fine_breaks, -1L) + tail(fine_breaks, -1L)) / 2

  # Basis and penalty
  basis <- bspline_basis(mids, a = a, b = b, ndx = ndx, degree = degree)
  K     <- basis$K
  P     <- diff_penalty(K = K, r = penalty_order)
  rankP <- K - penalty_order

  # Bin matrix
  C <- bin_matrix(wb, fine_breaks)

  # Starting value
  if (is.null(phi_start)) {
    phi0 <- rep(0, K)
  } else {
    if (length(phi_start) != K) stop("`phi_start` must have length K = ", K)
    phi0 <- .centre(as.numeric(phi_start))
  }

  # ---- Tau grid setup ----------------------------------------------------
  if (is.null(tau)) {
    tau_grid <- 10 ^ seq(-2, 6, length.out = 25L)
  } else {
    tau_grid <- sort(unique(as.numeric(tau)))
    if (any(tau_grid <= 0) || any(!is.finite(tau_grid))) {
      stop("`tau` must contain positive finite values.")
    }
  }

  # ---- Fit at each tau ---------------------------------------------------
  fits     <- vector("list", length(tau_grid))
  ic_vec   <- numeric(length(tau_grid))
  for (s in seq_along(tau_grid)) {
    fits[[s]] <- .pclm_score(phi0 = phi0, B = basis$B, C = C, m = m,
                             P = P, tau = tau_grid[s],
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)
    ic_vec[s] <- if (select == "BIC")
                   -2 * fits[[s]]$logL + log(sum(m)) * fits[[s]]$edf
                 else
                   -2 * fits[[s]]$logL + 2 * fits[[s]]$edf
  }

  # ---- Pick best ---------------------------------------------------------
  if (length(tau_grid) > 1L) {
    best_s   <- which.min(ic_vec)
    fit_best <- fits[[best_s]]
    tau_best <- tau_grid[best_s]
  } else {
    best_s   <- 1L
    fit_best <- fits[[1L]]
    tau_best <- tau_grid[1L]
  }

  out <- list(
    phi          = fit_best$phi,
    tau          = tau_best,
    tau_grid     = if (length(tau_grid) > 1L) tau_grid else NULL,
    ic           = ic_vec[best_s],
    ic_grid      = if (length(tau_grid) > 1L) ic_vec  else NULL,
    select       = select,
    logL         = fit_best$logL,
    edf          = fit_best$edf,
    vcov         = fit_best$vcov,
    pi           = fit_best$pi,
    gamma        = fit_best$gamma,
    fitted_counts = sum(m) * fit_best$gamma,
    grid         = fine_breaks,
    grid_mid     = mids,
    basis        = basis,
    C            = C,
    m            = m,
    wide_breaks  = wb,
    penalty_order = penalty_order,
    iter         = fit_best$iter,
    converged    = fit_best$converged,
    call         = cl
  )
  class(out) <- "pclm"
  out
}

# -----------------------------------------------------------------------------
# Internal: penalised scoring algorithm for a given tau
# -----------------------------------------------------------------------------
# Iteration: phi <- phi + (FI(phi) + tau P + ridge)^{-1} (grad - tau P phi)
# After each step, re-centre phi to satisfy the identifiability constraint
# sum(phi) = 0. The small ridge `eps_ridge` is added for numerical stability:
# both FI and P have the constant vector in their null spaces, so the system
# would be singular without it; the ridge shrinks the solution towards
# Σphi = 0 which we then enforce exactly by re-centring.
.pclm_score <- function(phi0, B, C, m, P, tau,
                        max_iter = 100L, tol = 1e-7,
                        eps_ridge = 1e-8,
                        verbose = FALSE) {

  K <- length(phi0)
  phi <- .centre(phi0)
  converged <- FALSE

  for (it in seq_len(max_iter)) {
    ev <- .pclm_eval(phi, B = B, C = C, m = m, compute_FI = TRUE)
    if (!is.finite(ev$logL)) {
      stop("Non-finite log-likelihood at iteration ", it,
           "; consider increasing tau or refining the grid.")
    }
    H   <- ev$FI + tau * P + eps_ridge * diag(K)
    rhs <- ev$grad - tau * (P %*% phi)
    step <- tryCatch(as.numeric(solve(H, rhs)),
                     error = function(e) {
                       # Fall back to a tiny extra ridge if needed
                       as.numeric(solve(H + 1e-4 * diag(K), rhs))
                     })
    # Damp step if it is huge (early iterations)
    step_norm <- max(abs(step))
    if (step_norm > 5) step <- step * (5 / step_norm)

    phi_new <- .centre(phi + step)
    delta   <- max(abs(phi_new - phi))
    phi     <- phi_new

    if (verbose) {
      message(sprintf("[pclm] iter %3d  logL = %.4f  ||step|| = %.2e",
                      it, ev$logL, delta))
    }
    if (delta < tol) { converged <- TRUE; break }
  }

  # Final evaluation (use latest phi)
  ev <- .pclm_eval(phi, B = B, C = C, m = m, compute_FI = TRUE)

  H        <- ev$FI + tau * P + eps_ridge * diag(K)
  vcov_mat <- tryCatch(solve(H),
                       error = function(e) solve(H + 1e-4 * diag(K)))
  # effective degrees of freedom = tr( H^{-1} FI )
  edf <- sum(diag(vcov_mat %*% ev$FI))

  list(
    phi       = phi,
    pi        = ev$pi,
    gamma     = ev$gamma,
    logL      = ev$logL,
    grad      = ev$grad,
    FI        = ev$FI,
    vcov      = vcov_mat,
    edf       = edf,
    iter      = it,
    converged = converged
  )
}
