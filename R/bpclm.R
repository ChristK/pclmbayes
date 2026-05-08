# =============================================================================
# Bayesian penalised composite link model fit
# -----------------------------------------------------------------------------
# Implements the Metropolis-within-Gibbs sampler of Lambert and Eilers (2009),
# Section 3.3:
#
#   * Phi step: modified Langevin-Hastings (MALA) proposal
#       phi* ~ N(phi + 0.5 delta Sigma G, delta Sigma)
#     where G = grad log p(phi | tau, D), and Sigma is a fixed positive
#     definite matrix chosen to approximate the posterior covariance of
#     phi (typically the inverse of (FI + tau P) at the frequentist MLE).
#   * Tau step: Gibbs draw tau | phi ~ Gamma(a + 0.5 R(P), b + 0.5 phi'P phi)
#
#   The acceptance log-ratio uses
#     q(phi*, phi) / q(phi, phi*)
#       = exp(- 0.5 (G* + G)' (phi* - phi + (delta/4) Sigma (G* - G)))
#   (Lambert & Eilers 2009, Section 3.3.1).
#
#   Adaptive tuning of delta during burn-in targets a 0.57 acceptance rate
#   (Roberts and Rosenthal 1998; Atchade and Rosenthal 2005).
#
#   Optional shape constraints are enforced as in Eq. (7) of the paper, by
#   rejecting any proposal whose induced fine-grid density violates the
#   indicator I_Phi.
# =============================================================================

#' Bayesian penalised composite link model for grouped data
#'
#' Fits the model of Lambert and Eilers (2009) by Markov chain Monte
#' Carlo. The latent log-density is modelled by P-splines on a fine
#' grid; wide-bin counts have a multinomial likelihood; the spline
#' coefficients have a smoothness prior controlled by a precision
#' parameter \eqn{\tau} with a vague \eqn{\Gamma(a, b)} hyperprior.
#'
#' The sampler is a Metropolis-within-Gibbs scheme: a modified
#' Langevin-Hastings (MALA) update for \eqn{\phi} followed by a Gibbs
#' draw for \eqn{\tau}. The proposal covariance for \eqn{\phi} is a
#' fixed matrix \code{Sigma} approximating the posterior covariance,
#' obtained from a frequentist warm-start fit (\code{\link{pclm}}) at a
#' user-chosen \eqn{\tau}.
#'
#' Optional shape constraints (\code{"unimodal"}, \code{"logconcave"}
#' and/or \code{"monotonic"}) are imposed by rejecting any proposal
#' whose induced density violates them, in line with Eq. (7) of the
#' paper.
#'
#' @inheritParams pclm
#' @param niter Total number of MCMC iterations (including burn-in).
#'   Default 5000.
#' @param burnin Number of initial iterations to discard. Default
#'   \code{floor(niter / 5)}.
#' @param thin Thinning interval (default 1). Only every \code{thin}-th
#'   post-burn-in draw is retained.
#' @param adapt Number of additional iterations of adaptive tuning of
#'   the step size \code{delta} performed before the main run. Default
#'   500.
#' @param tau_init Initial value of the precision \eqn{\tau}. If
#'   \code{NULL} (the default), the BIC-selected value from a
#'   frequentist warm-start fit is used.
#' @param a_tau,b_tau Hyperparameters of the \eqn{\Gamma(a, b)}
#'   hyperprior on \eqn{\tau}. Defaults \code{1e-4} each give a near-
#'   flat improper prior on \eqn{\log \tau}.
#' @param delta_init Initial step size for the Langevin proposal. If
#'   \code{NULL}, the value \eqn{1.65^2 / K^{1/3}} from Roberts and
#'   Rosenthal (1998) is used.
#' @param target_accept Target acceptance rate during adaptation
#'   (default 0.57).
#' @param Sigma Optional positive-definite \eqn{K \times K} proposal
#'   covariance. If \code{NULL}, the frequentist warm-start
#'   variance-covariance matrix is used.
#' @param shape Optional character vector of shape constraints. Any
#'   subset of \code{"unimodal"}, \code{"logconcave"} (or
#'   \code{"log-concave"}) and \code{"monotonic"}. Default \code{NULL}
#'   (unconstrained).
#' @param shape_args A list of additional arguments for the shape
#'   indicators. Currently used only for monotonicity, where
#'   \code{shape_args = list(direction = "decreasing")} (or
#'   \code{"increasing"}, default \code{"either"}) restricts the
#'   sense.
#' @param phi_init Optional starting value for the chain on \eqn{\phi}.
#'   Defaults to the warm-start frequentist MLE.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical: print progress every 10\% of iterations.
#' @param x A \code{"bpclm"} object (in \code{print}, \code{plot} and
#'   \code{quantile} methods).
#' @param object A \code{"bpclm"} object (in \code{summary}, \code{coef},
#'   \code{fitted} and \code{predict} methods).
#' @param digits Number of significant digits used in \code{print}.
#' @param newdata Optional numeric vector of points at which the fitted
#'   density is evaluated by \code{predict}. If missing, the posterior
#'   mean (or sample) density on the fine grid is returned.
#' @param ... Further arguments. Currently ignored by all methods.
#'
#' @return An object of class \code{"bpclm"}, a list with components:
#'   \describe{
#'     \item{phi}{Posterior mean of \eqn{\phi} (length \eqn{K}).}
#'     \item{phi_chain}{\code{nsim x K} matrix of post-burn-in,
#'       post-thin draws of \eqn{\phi}.}
#'     \item{tau_chain}{Numeric vector of post-burn-in draws of
#'       \eqn{\tau}.}
#'     \item{pi_chain}{\code{nsim x I} matrix of fine-grid probabilities
#'       (each row a draw).}
#'     \item{pi}{Posterior mean of \eqn{\pi} (latent grid density,
#'       summing to 1).}
#'     \item{pi_lower, pi_upper}{Pointwise 90\% (or other, see
#'       \code{cred_level}) credible-interval limits for \eqn{\pi}.}
#'     \item{cred_level}{Credible level used for the bands (default
#'       0.90).}
#'     \item{accept}{Overall acceptance rate of the \eqn{\phi} step
#'       across post-adaptation iterations.}
#'     \item{delta}{Final tuned step size.}
#'     \item{adapt_path}{Numeric vector tracking \eqn{\delta} over the
#'       adaptation phase (for diagnostics).}
#'     \item{warmstart}{The frequentist \code{pclm} fit used for warm
#'       starting (or \code{NULL} if all relevant arguments were
#'       supplied directly).}
#'     \item{grid, grid_mid, basis, C, m, wide_breaks, penalty_order,
#'       a_tau, b_tau, shape, call}{Bookkeeping.}
#'   }
#'
#' @references
#' Lambert, P. and Eilers, P. H. C. (2009). Bayesian density estimation
#' from grouped continuous data. \emph{Computational Statistics and
#' Data Analysis}, 53(4), 1388--1399.
#'
#' Roberts, G. O. and Rosenthal, J. S. (1998). Optimal scaling of
#' discrete approximations to Langevin diffusions. \emph{Journal of the
#' Royal Statistical Society, Series B}, 60(1), 255--268.
#'
#' Atchade, Y. F. and Rosenthal, J. S. (2005). On adaptive Markov chain
#' Monte Carlo algorithms. \emph{Bernoulli}, 11(5), 815--828.
#'
#' @examples
#' \donttest{
#' data(bloodlead)
#' fit <- bpclm(m = bloodlead$count,
#'              wide_breaks = with(bloodlead, cbind(lower, upper)),
#'              a = 0, b = 80,
#'              ngrid = 80, ndx = 17,
#'              niter = 2000, burnin = 500,
#'              shape = "unimodal", seed = 1)
#' summary(fit)
#' }
#'
#' @name bpclm
#' @rdname bpclm
#' @export
bpclm <- function(m, wide_breaks,
                  a = NULL, b = NULL,
                  ngrid = 100L,
                  ndx = 17L, degree = 3L,
                  penalty_order = 3L,
                  niter = 5000L,
                  burnin = NULL,
                  thin = 1L,
                  adapt = 500L,
                  tau_init = NULL,
                  a_tau = 1e-4, b_tau = 1e-4,
                  delta_init = NULL,
                  target_accept = 0.57,
                  Sigma = NULL,
                  shape = NULL,
                  shape_args = list(),
                  phi_init = NULL,
                  seed = NULL,
                  verbose = FALSE) {

  cl <- match.call()

  if (!is.null(seed)) set.seed(seed)

  # ---- Frequentist warm start -------------------------------------------
  warm <- pclm(m = m, wide_breaks = wide_breaks,
               a = a, b = b,
               ngrid = ngrid, ndx = ndx, degree = degree,
               penalty_order = penalty_order)

  basis <- warm$basis
  Bmat  <- basis$B
  Cmat  <- warm$C
  K     <- basis$K
  P     <- diff_penalty(K, r = penalty_order)
  rankP <- K - penalty_order
  m_obs <- as.numeric(m)
  mids  <- warm$grid_mid
  fine_breaks <- warm$grid

  # ---- Sampler initial state --------------------------------------------
  if (is.null(phi_init)) {
    phi <- warm$phi
  } else {
    if (length(phi_init) != K) stop("`phi_init` must have length K = ", K)
    phi <- .centre(as.numeric(phi_init))
  }
  tau <- if (is.null(tau_init)) warm$tau else as.numeric(tau_init)

  # Proposal covariance Sigma
  if (is.null(Sigma)) {
    Sigma <- warm$vcov
  }
  # Symmetrise + small jitter for stability
  Sigma <- (Sigma + t(Sigma)) / 2
  # Cholesky (with retry-with-jitter)
  L <- tryCatch(t(chol(Sigma)),
                error = function(e) {
                  jitter <- 1e-8 * mean(diag(Sigma))
                  t(chol(Sigma + jitter * diag(K)))
                })

  # Step size delta
  if (is.null(delta_init)) {
    delta <- 1.65 ^ 2 / K ^ (1 / 3)  # Roberts & Rosenthal (1998)
  } else {
    delta <- as.numeric(delta_init)
  }
  # Adaptive-tuning bounds
  delta_lo <- 1e-4
  delta_hi <- 1e4
  burnin   <- if (is.null(burnin)) floor(niter / 5) else as.integer(burnin)

  # ---- Helper: log posterior of phi (up to terms in tau only) -----------
  # Returns ell(phi) + log prior(phi | tau)  =  sum_j m_j log gamma_j
  #                                            - 0.5 tau phi' P phi
  # plus the gradient.
  posterior_eval <- function(phi, tau) {
    ev   <- .pclm_eval(phi, B = Bmat, C = Cmat, m = m_obs, compute_FI = FALSE)
    pen  <- 0.5 * tau * as.numeric(crossprod(phi, P %*% phi))
    list(pi   = ev$pi,
         logp = ev$logL - pen,
         grad = ev$grad - tau * as.numeric(P %*% phi))
  }

  # ---- Helper: a single Langevin proposal step --------------------------
  # current_state must contain pi, logp, grad (post-shape-check)
  step_phi <- function(phi, current_state, delta, tau) {
    G   <- current_state$grad
    mu  <- as.numeric(phi + 0.5 * delta * (Sigma %*% G))
    # Sample z ~ N(0, delta * Sigma) by L %*% sqrt(delta) %*% N(0, I)
    z   <- as.numeric(L %*% rnorm(K, sd = sqrt(delta)))
    phi_star <- mu + z
    # Re-centre to the identifiability hyperplane Σphi = 0
    phi_star <- .centre(phi_star)
    new_state <- posterior_eval(phi_star, tau)

    # Shape constraint
    if (!.shape_ok(new_state$pi, shape, shape_args)) {
      return(list(phi = phi, state = current_state, accepted = FALSE,
                  logA = -Inf))
    }
    if (!is.finite(new_state$logp)) {
      return(list(phi = phi, state = current_state, accepted = FALSE,
                  logA = -Inf))
    }

    # Acceptance log-ratio: log p(phi*) - log p(phi)
    #                     - 0.5 (G* + G)' (phi* - phi + (delta/4) Sigma (G* - G))
    Gstar <- new_state$grad
    diff_phi <- phi_star - phi
    diff_G   <- Gstar - G
    correction <- diff_phi + (delta / 4) * as.numeric(Sigma %*% diff_G)
    qratio_log <- -0.5 * sum((Gstar + G) * correction)
    logA <- (new_state$logp - current_state$logp) + qratio_log
    if (is.finite(logA) && log(runif(1)) < logA) {
      list(phi = phi_star, state = new_state, accepted = TRUE,
           logA = logA)
    } else {
      list(phi = phi, state = current_state, accepted = FALSE, logA = logA)
    }
  }

  # ---- Initial state ----------------------------------------------------
  state <- posterior_eval(phi, tau)
  if (!.shape_ok(state$pi, shape, shape_args)) {
    warning("Initial state violates shape constraint(s); attempting to find a ",
            "feasible starting point.")
    feasible <- FALSE
    # Strategy 1: re-fit the frequentist PCLM at progressively larger tau
    # (more smoothing) until the resulting density satisfies the shape
    # constraint. Iterate from the SMALLEST oversmoothing upwards so that
    # we land at the MILDEST oversmoothing that succeeds -- this keeps phi
    # close to the data, rather than collapsing to a near-flat density.
    # When found, also refresh Sigma and L so the proposal geometry
    # matches the chain's new location.
    big_tau_grid <- 10 ^ seq(0, 8, length.out = 17L)  # 1, ~3.16, 10, ..., 1e8
    big_tau_grid <- big_tau_grid[big_tau_grid >= warm$tau]   # only oversmooth
    for (tau_try in big_tau_grid) {
      warm_smoothed <- tryCatch(
        pclm(m = m_obs, wide_breaks = warm$wide_breaks,
             a = warm$grid[1L], b = warm$grid[length(warm$grid)],
             ngrid = length(warm$grid_mid),
             ndx = basis$ndx, degree = basis$degree,
             penalty_order = penalty_order,
             tau = tau_try),
        error = function(e) NULL)
      if (is.null(warm_smoothed)) next
      phi_try <- warm_smoothed$phi
      state_try <- posterior_eval(phi_try, tau)
      if (.shape_ok(state_try$pi, shape, shape_args)) {
        phi <- phi_try; state <- state_try; feasible <- TRUE
        # Refresh proposal covariance to match the new local geometry.
        # Plain `<-` here, because step_phi is a closure over bpclm's
        # frame, so updating Sigma and L locally is what step_phi will
        # see on the next call.
        Sigma <- (warm_smoothed$vcov + t(warm_smoothed$vcov)) / 2
        L <- tryCatch(t(chol(Sigma)),
                      error = function(e) {
                        jit <- 1e-8 * mean(diag(Sigma))
                        t(chol(Sigma + jit * diag(K)))
                      })
        if (is.null(delta_init)) delta <- 1.65 ^ 2 / K ^ (1 / 3)
        break
      }
    }
    # Strategy 2: random perturbations of the smoothed start.
    if (!feasible) {
      for (try_i in seq_len(200L)) {
        phi_try <- .centre(phi + 0.1 * try_i * rnorm(K) / sqrt(K))
        state_try <- posterior_eval(phi_try, tau)
        if (.shape_ok(state_try$pi, shape, shape_args)) {
          phi <- phi_try; state <- state_try; feasible <- TRUE; break
        }
      }
    }
    # Strategy 3: a manifest single-peak prior centred at the highest-
    # density wide bin; effectively a quadratic envelope on phi.
    if (!feasible) {
      bin_dens <- (m_obs / sum(m_obs)) /
                    (warm$wide_breaks[, 2L] - warm$wide_breaks[, 1L])
      mode_x   <- mean(warm$wide_breaks[which.max(bin_dens), ])
      sd_guess <- diff(range(warm$grid)) / 6
      env_phi  <- -0.5 * ((basis$knots[seq_len(K) + basis$degree] - mode_x) /
                            sd_guess) ^ 2
      env_phi  <- .centre(env_phi)
      state_try <- posterior_eval(env_phi, tau)
      if (.shape_ok(state_try$pi, shape, shape_args)) {
        phi <- env_phi; state <- state_try; feasible <- TRUE
      }
    }
    if (!feasible) {
      stop("Could not locate a feasible starting point under the requested ",
           "shape constraint(s). Consider relaxing them or supplying ",
           "`phi_init` manually.")
    }
  }

  # ---- Adaptation phase: update delta after each iteration --------------
  adapt_path <- numeric(adapt)
  if (adapt > 0L) {
    for (t in seq_len(adapt)) {
      step <- step_phi(phi, state, delta, tau)
      phi   <- step$phi
      state <- step$state
      acc_t <- if (step$accepted) 1 else 0
      # Tau Gibbs step (also during adaptation, so its scale settles in)
      tau   <- rgamma(1,
                      shape = a_tau + 0.5 * rankP,
                      rate  = b_tau + 0.5 * as.numeric(crossprod(phi,
                                                                P %*% phi)))
      state <- posterior_eval(phi, tau)
      # Update delta on the sqrt scale (Atchade & Rosenthal 2005)
      sd_old <- sqrt(delta)
      gamma_t <- 1 / t
      sd_new  <- sd_old + gamma_t * (acc_t - target_accept)
      sd_new  <- min(max(sd_new, sqrt(delta_lo)), sqrt(delta_hi))
      delta   <- sd_new ^ 2
      adapt_path[t] <- delta
    }
  }

  # ---- Main run with FIXED delta (so the chain is genuinely Markov) -----
  n_keep <- floor((niter - burnin) / thin)
  if (n_keep <= 0L) stop("`niter - burnin` must give at least one kept draw.")
  phi_chain <- matrix(NA_real_, nrow = n_keep, ncol = K)
  tau_chain <- numeric(n_keep)
  pi_chain  <- matrix(NA_real_, nrow = n_keep, ncol = length(state$pi))

  acc_count <- 0L
  acc_total <- 0L
  store_idx <- 0L
  for (it in seq_len(niter)) {
    step <- step_phi(phi, state, delta, tau)
    phi   <- step$phi
    state <- step$state
    if (it > burnin) {
      acc_count <- acc_count + as.integer(step$accepted)
      acc_total <- acc_total + 1L
    }

    # Tau Gibbs
    tau   <- rgamma(1,
                    shape = a_tau + 0.5 * rankP,
                    rate  = b_tau + 0.5 * as.numeric(crossprod(phi,
                                                              P %*% phi)))
    state <- posterior_eval(phi, tau)

    # Store post-burn-in, post-thin
    if (it > burnin && ((it - burnin) %% thin) == 0L) {
      store_idx <- store_idx + 1L
      if (store_idx <= n_keep) {
        phi_chain[store_idx, ] <- phi
        tau_chain[store_idx]   <- tau
        pi_chain[store_idx, ]  <- state$pi
      }
    }

    if (verbose && (it %% max(1L, niter %/% 10L) == 0L)) {
      message(sprintf("[bpclm] iter %d / %d   tau = %.3g   acc = %.2f",
                      it, niter, tau,
                      if (acc_total > 0) acc_count / acc_total else NA))
    }
  }
  phi_chain <- phi_chain[seq_len(store_idx), , drop = FALSE]
  tau_chain <- tau_chain[seq_len(store_idx)]
  pi_chain  <- pi_chain[seq_len(store_idx),  , drop = FALSE]

  # ---- Posterior summaries ----------------------------------------------
  cred_level <- 0.90
  pi_mean  <- colMeans(pi_chain)
  pi_lo    <- apply(pi_chain, 2L,
                    function(x) quantile(x, (1 - cred_level) / 2))
  pi_hi    <- apply(pi_chain, 2L,
                    function(x) quantile(x, 1 - (1 - cred_level) / 2))

  out <- list(
    phi          = colMeans(phi_chain),
    phi_chain    = phi_chain,
    tau_chain    = tau_chain,
    pi_chain     = pi_chain,
    pi           = pi_mean,
    pi_lower     = pi_lo,
    pi_upper     = pi_hi,
    cred_level   = cred_level,
    accept       = if (acc_total > 0) acc_count / acc_total else NA_real_,
    delta        = delta,
    adapt_path   = adapt_path,
    warmstart    = warm,
    basis        = basis,
    C            = Cmat,
    grid         = fine_breaks,
    grid_mid     = mids,
    m            = m_obs,
    wide_breaks  = warm$wide_breaks,
    penalty_order = penalty_order,
    a_tau        = a_tau,
    b_tau        = b_tau,
    shape        = shape,
    shape_args   = shape_args,
    niter        = niter,
    burnin       = burnin,
    thin         = thin,
    call         = cl
  )
  class(out) <- "bpclm"
  out
}
