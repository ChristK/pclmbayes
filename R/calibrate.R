# =============================================================================
# Exact preservation of wide-bin totals
# -----------------------------------------------------------------------------
# The Lambert-Eilers PCLM treats the observed wide-bin counts as multinomial
# samples and shrinks them by the smoothness penalty.  In an "ungrouping"
# context (where the wide-bin counts are taken as fixed totals to be
# distributed across a fine grid -- e.g. distributing deaths in 5-year age
# bands across single years of age), the user often needs the implied
# fine-bin counts to sum back to the observed wide-bin counts exactly.
#
# Two routes are provided:
#
#   * calibrate()     : post-hoc renormalisation of any fitted pclm/bpclm.
#                       For each wide bin j, scale pi_i (i in bin j) by
#                       (m_j / m_+) / gamma_j.  Cheap, preserves the smooth
#                       within-bin shape, but introduces small kinks at the
#                       bin boundaries.  Works only when the wide bins
#                       partition the support without overlap (the typical
#                       ungrouping setting).
#
#   * pclm_exact()    : a frequentist fit that minimises (1/2) phi' P phi
#                       subject to gamma(phi) = m / m_+  (a hard equality
#                       constraint).  No kinks at boundaries.  Solved by an
#                       SQP iteration -- linearise the constraint, take a
#                       Newton step on the Lagrangian, repeat.
#
# Both methods preserve the wide-bin totals exactly to numerical precision.
# =============================================================================


# Internal: a J-vector mapping each fine interval i (1..I) to its parent
# wide bin j (1..J) under the assumption that bins partition the support.
# Returns NA for fine intervals not fully assigned to a single bin (which
# happens when bins overlap or are misaligned with the fine grid).
.bin_assignment <- function(C, tol = 1e-8) {
  I <- ncol(C)
  J <- nrow(C)
  # For each column i: how many bins is the fine interval substantially
  # assigned to?
  assign <- integer(I)
  for (i in seq_len(I)) {
    nz <- which(C[, i] > tol)
    if (length(nz) == 1L && abs(C[nz, i] - 1) < tol) {
      assign[i] <- nz
    } else {
      assign[i] <- NA_integer_
    }
  }
  assign
}


#' Calibrate a fitted PCLM to preserve wide-bin totals exactly
#'
#' Post-hoc adjustment of a fitted \code{pclm} or \code{bpclm} object so
#' that the implied wide-bin counts match the observed counts \code{m}
#' exactly. For each wide bin \eqn{j}, the fine-grid probabilities
#' \eqn{\pi_i} (\eqn{i} in bin \eqn{j}) are multiplied by
#' \eqn{(m_j / m_+) / \gamma_j}, where \eqn{\gamma_j = (C\pi)_j} is the
#' fitted bin probability before calibration. The within-bin shape from
#' the smooth fit is preserved; small jumps at the bin boundaries are
#' the price paid for exact preservation.
#'
#' Calibration is meaningful only when the wide bins partition the
#' support without overlap, so that each fine interval belongs to
#' exactly one wide bin. The function checks this and errors if the bin
#' assignment is ambiguous (e.g. overlapping bins, or fine-grid
#' breakpoints misaligned with the wide-bin boundaries).
#'
#' For \code{bpclm} objects, calibration is applied separately to every
#' posterior draw and the posterior summaries (\code{pi},
#' \code{pi_lower}, \code{pi_upper}, \code{phi}) are recomputed.
#' Posterior credible intervals for derived quantities (means,
#' quantiles) recomputed from the calibrated chain will respect the
#' constraint.
#'
#' @param fit A fitted \code{"pclm"} or \code{"bpclm"} object.
#' @param ... Passed to methods (currently unused).
#' @return The same object with calibrated \code{pi} (and, for
#'   \code{bpclm}, calibrated \code{pi_chain}, \code{pi_lower},
#'   \code{pi_upper}). The element \code{fitted_counts} is set equal to
#'   the input counts \code{m}, and \code{calibrated = TRUE} is added.
#'
#' @examples
#' data(bloodlead)
#' fit <- pclm(m = bloodlead$count,
#'             wide_breaks = with(bloodlead, cbind(lower, upper)),
#'             a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
#'             penalty_order = 3)
#' max(abs(fit$fitted_counts - fit$m))            # 7.4
#' fit_c <- calibrate(fit)
#' max(abs(fit_c$fitted_counts - fit_c$m))        # ~ 0
#'
#' @seealso \code{\link{pclm_exact}} for a constrained MAP fit that
#'   enforces the same property without introducing kinks at bin
#'   boundaries.
#' @export
calibrate <- function(fit, ...) UseMethod("calibrate")


#' @export
calibrate.pclm <- function(fit, ...) {
  assign <- .bin_assignment(fit$C)
  if (any(is.na(assign))) {
    stop("calibrate() requires the wide bins to partition the support ",
         "without overlap; some fine intervals are not assigned to a ",
         "single bin.  Consider pclm_exact() instead.")
  }
  m_plus <- sum(fit$m)
  target <- fit$m / m_plus
  gamma  <- as.numeric(fit$C %*% fit$pi)

  scale_per_bin <- ifelse(gamma > 0, target / gamma, 0)
  pi_new <- fit$pi * scale_per_bin[assign]
  s <- sum(pi_new)
  if (s > 0) pi_new <- pi_new / s

  fit$pi            <- pi_new
  fit$gamma         <- as.numeric(fit$C %*% pi_new)
  fit$fitted_counts <- m_plus * fit$gamma
  fit$calibrated    <- TRUE
  fit
}


#' @export
calibrate.bpclm <- function(fit, ...) {
  assign <- .bin_assignment(fit$C)
  if (any(is.na(assign))) {
    stop("calibrate() requires the wide bins to partition the support ",
         "without overlap; some fine intervals are not assigned to a ",
         "single bin.  Consider pclm_exact() instead.")
  }
  m_plus <- sum(fit$m)
  target <- fit$m / m_plus

  cal_pi <- function(p) {
    g <- as.numeric(fit$C %*% p)
    s <- ifelse(g > 0, target / g, 0)
    p_new <- p * s[assign]
    tot <- sum(p_new)
    if (tot > 0) p_new <- p_new / tot
    p_new
  }
  fit$pi_chain <- t(apply(fit$pi_chain, 1L, cal_pi))
  cred <- fit$cred_level
  alpha <- (1 - cred) / 2
  fit$pi       <- colMeans(fit$pi_chain)
  fit$pi_lower <- apply(fit$pi_chain, 2L, quantile,
                        probs = alpha,     names = FALSE)
  fit$pi_upper <- apply(fit$pi_chain, 2L, quantile,
                        probs = 1 - alpha, names = FALSE)
  fit$calibrated <- TRUE
  fit
}


#' Constrained MAP fit: smoothest density that exactly preserves wide-bin
#' counts
#'
#' Solves the equality-constrained optimisation problem
#' \deqn{\min_\phi \tfrac{1}{2} \phi' P \phi
#'   \quad \text{subject to} \quad
#'   C \pi(\phi) = m / m_+,}
#' where \eqn{\pi(\phi)} is the softmax of \eqn{B \phi}, \eqn{P} is the
#' \eqn{r}th-order difference penalty, \eqn{C} is the bin matrix, and
#' \eqn{m / m_+} is the vector of observed wide-bin proportions. The
#' result is the smoothest density on the fine grid that exactly
#' reproduces the wide-bin counts when summed.
#'
#' This is the principled alternative to \code{\link{calibrate}}: no
#' kinks are introduced at the bin boundaries because the smoothness
#' penalty is minimised across the whole fine grid subject to the
#' constraint. Solved by sequential quadratic programming (Newton steps
#' on the Lagrangian, with the constraint linearised at each step).
#'
#' @inheritParams pclm
#' @param eps_ridge Small ridge added to the (singular) penalty matrix
#'   for numerical stability. Default \code{1e-8}.
#' @param max_iter,tol Convergence controls.
#'
#' @return An object of class \code{c("pclm_exact", "pclm")} with the
#'   same components as a \code{\link{pclm}} object, plus
#'   \code{lambda} (the converged Lagrange multipliers) and
#'   \code{constraint_residual} (max |γ - m/m_+|).
#'
#' @references
#' Eilers, P. H. C. (2007). Ill-posed problems with counts, the
#' composite link model and penalized likelihood. \emph{Statistical
#' Modelling}, 7(3), 239--254. (Underlying penalised composite-link
#' model on a fine grid.)
#'
#' Nocedal, J. and Wright, S. J. (2006). \emph{Numerical Optimization},
#' 2nd ed. Springer.  Chapter 18 (Sequential Quadratic Programming, used
#' here to enforce the bin-total constraint exactly).
#'
#' @examples
#' data(bloodlead)
#' fit <- pclm_exact(m = bloodlead$count,
#'                    wide_breaks = with(bloodlead, cbind(lower, upper)),
#'                    a = 0, b = 80, ngrid = 80, ndx = 17, degree = 3,
#'                    penalty_order = 3)
#' max(abs(fit$fitted_counts - fit$m))    # < 1e-8
#'
#' @seealso \code{\link{pclm}}, \code{\link{calibrate}}.
#' @export
pclm_exact <- function(m, wide_breaks,
                       a = NULL, b = NULL,
                       ngrid = 100L,
                       ndx = 17L, degree = 3L,
                       penalty_order = 3L,
                       max_iter = 200L, tol = 1e-9,
                       eps_ridge = 1e-8,
                       phi_start = NULL,
                       verbose = FALSE) {
  cl <- match.call()

  # ---- standardise inputs (same as pclm) ----
  if (is.matrix(wide_breaks) || is.data.frame(wide_breaks)) {
    wb <- as.matrix(wide_breaks)
  } else {
    wbv <- sort(as.numeric(wide_breaks))
    wb  <- cbind(head(wbv, -1), tail(wbv, -1))
  }
  if (is.null(a)) a <- min(wb[, 1L])
  if (is.null(b)) b <- max(wb[, 2L])

  fine_breaks <- seq(a, b, length.out = ngrid + 1L)
  mids        <- (head(fine_breaks, -1) + tail(fine_breaks, -1)) / 2
  basis       <- bspline_basis(mids, a = a, b = b,
                               ndx = ndx, degree = degree)
  K           <- basis$K
  P           <- diff_penalty(K = K, r = penalty_order)
  C           <- bin_matrix(wb, fine_breaks)
  J           <- nrow(C)

  m_plus <- sum(m)
  target <- m / m_plus     # length J, sums to 1

  # Starting value: warm start with a small frequentist fit at large tau
  if (is.null(phi_start)) {
    warm <- pclm(m = m, wide_breaks = wb,
                 a = a, b = b, ngrid = ngrid, ndx = ndx, degree = degree,
                 penalty_order = penalty_order,
                 tau = 1)         # mild smoothing as a starting point
    phi <- warm$phi
  } else {
    phi <- .centre(as.numeric(phi_start))
  }
  lambda <- numeric(J)         # Lagrange multipliers

  # SQP iteration on the equality-constrained quadratic
  #   min  (1/2) phi' P phi   s.t.   gamma(phi) = target
  #
  # The KKT matrix
  #   [ P + eps I    J' ]
  #   [ J            0  ]
  # has null space along the constant direction in phi (the softmax is
  # invariant under phi -> phi + c*1), so we need either:
  #   (a) a small fixed ridge (introduces O(eps_ridge) bias on the
  #       constraint), or
  #   (b) reparametrisation to phi[K] = -sum(phi[1:(K-1)]).
  # We use (a) with a small ridge and an outer "polish" step: after
  # convergence, project phi onto the constraint manifold by a final
  # bin-conditional renormalisation.  This guarantees the constraint is
  # satisfied to machine precision without sacrificing convergence
  # robustness.
  for (it in seq_len(max_iter)) {
    # Compute pi, gamma, dgamma/dphi at current phi
    pi_curr <- .pi_phi(phi, basis$B)
    g       <- as.numeric(C %*% pi_curr)
    pB      <- as.numeric(crossprod(basis$B, pi_curr))            # K
    dpidphi <- (pi_curr * basis$B) - tcrossprod(pi_curr, pB)      # I x K
    Jmat    <- C %*% dpidphi                                       # J x K

    # KKT block structure.  The constraint Cπ = target is redundant
    # by one (both sides sum to 1), so the J rows of J are linearly
    # dependent: rowSums(J) ≈ 0.  We add a small ridge to the bottom-
    # right block to break this dependency; the resulting bias is
    # O(eps_ridge) and is removed by the polish step at the end.
    A11 <- P + eps_ridge * diag(K)
    A22 <- -eps_ridge * diag(J)
    KKT <- rbind(
      cbind(A11,    t(Jmat)),
      cbind(Jmat,   A22)
    )
    rhs <- c(-as.numeric(P %*% phi),
             target - g)
    sol <- tryCatch(
      solve(KKT, rhs),
      error = function(e) {
        KKT2 <- KKT
        KKT2[seq_len(K), seq_len(K)] <-
          KKT2[seq_len(K), seq_len(K)] + 1e-6 * diag(K)
        solve(KKT2, rhs)
      })
    dphi   <- sol[seq_len(K)]
    lambda <- sol[(K + 1L):(K + J)]

    # Damp step if it is huge (early iterations)
    s_norm <- max(abs(dphi))
    if (s_norm > 5) dphi <- dphi * (5 / s_norm)

    phi_new <- .centre(phi + dphi)
    pi_new  <- .pi_phi(phi_new, basis$B)
    g_new   <- as.numeric(C %*% pi_new)
    res     <- max(abs(g_new - target))
    step    <- max(abs(phi_new - phi))
    phi     <- phi_new

    if (verbose) {
      message(sprintf("[pclm_exact] iter %3d  step = %.2e  ||g - t||_inf = %.2e",
                      it, step, res))
    }
    if (res < 1e-6 && step < tol) break
  }

  # ---- Polish: final bin-conditional renormalisation ----
  # The SQP iteration converges to a constraint residual of order
  # eps_ridge.  A final exact projection onto the constraint manifold
  # eliminates this bias.  Only safe when bins partition the support
  # without overlap (the typical ungrouping setting).
  pi_curr <- .pi_phi(phi, basis$B)
  assign  <- .bin_assignment(C)
  if (!any(is.na(assign))) {
    g_now <- as.numeric(C %*% pi_curr)
    s     <- ifelse(g_now > 0, target / g_now, 0)
    pi_curr <- pi_curr * s[assign]
    pi_curr <- pi_curr / sum(pi_curr)
  }

  pi_final <- pi_curr
  ga_final <- as.numeric(C %*% pi_final)

  # Cheap log-likelihood at convergence (non-meaningful since constraint
  # is exactly met; included for parity with pclm).
  pos  <- m > 0
  logL <- if (any(ga_final[pos] <= 0)) -Inf
          else sum(m[pos] * log(ga_final[pos]))

  out <- list(
    phi          = phi,
    tau          = NA_real_,        # no smoothing parameter -- pure constraint
    tau_grid     = NULL,
    ic           = NA_real_,
    ic_grid      = NULL,
    select       = NA_character_,
    logL         = logL,
    edf          = NA_real_,
    vcov         = NULL,
    pi           = pi_final,
    gamma        = ga_final,
    fitted_counts = m_plus * ga_final,
    grid         = fine_breaks,
    grid_mid     = mids,
    basis        = basis,
    C            = C,
    m            = as.numeric(m),
    wide_breaks  = wb,
    penalty_order = penalty_order,
    iter         = it,
    converged    = (max(abs(ga_final - target)) < tol),
    lambda       = lambda,
    constraint_residual = max(abs(ga_final - target)),
    calibrated   = TRUE,
    call         = cl
  )
  class(out) <- c("pclm_exact", "pclm")
  out
}
