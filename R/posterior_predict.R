# =============================================================================
# Posterior predictive distribution of single-year (fine-grid) counts
# -----------------------------------------------------------------------------
# Most users of an ungrouping algorithm want to know not just the smooth
# fitted density but the plausible single-year (or single-fine-cell) counts.
# Two distinct sources of uncertainty are involved:
#
#   * "rate"        : the smooth latent density m_+ * pi_y.  Captures
#                     epistemic uncertainty about the underlying rate from
#                     the smoothness assumption (Bayesian) or none at all
#                     (frequentist point estimate).
#   * "predictive"  : the *realised* per-cell counts.  Given the wide-bin
#                     totals m_j, the realised per-cell counts within
#                     band j follow a multinomial with size m_j and
#                     within-band probabilities pi_y / gamma_j.  This is
#                     the right object when the user wants plausible
#                     single-year breakdowns of the observed band counts.
#
# `posterior_predict()` returns a sample from either distribution, plus
# convenient mean and credible (or prediction) intervals.
#
# For Bayesian fits (`bpclm`), the rate distribution is the posterior
# chain on pi; the predictive distribution combines that chain with the
# multinomial sampling step.  For frequentist point estimates (`pclm`,
# `pclm_exact`), the rate is a single value (no uncertainty); the
# predictive distribution is a parametric (multinomial) sample around
# the point estimate.
# =============================================================================


#' Posterior or parametric predictive distribution of fine-cell counts
#'
#' Generates draws from the posterior or parametric predictive
#' distribution of the latent fine-cell (typically single-year) counts.
#' Two flavours are supported:
#'
#' \describe{
#'   \item{\code{type = "rate"}}{The smooth latent rate \eqn{m_+ \pi}.
#'     For \code{bpclm} fits this returns a chain of rate draws (one per
#'     posterior draw of \eqn{\pi}). For \code{pclm} and
#'     \code{pclm_exact} fits this returns the point estimate
#'     \eqn{m_+ \pi} (no uncertainty). Mean = posterior mean rate; the
#'     credible band is narrow at large \eqn{N}.}
#'   \item{\code{type = "predictive"} (default)}{The actual realised
#'     fine-cell counts. Within each wide bin \eqn{j}, the conditional
#'     distribution of the fine-cell counts given the band total
#'     \eqn{m_j} and the within-band probabilities \eqn{\pi_y/\gamma_j}
#'     is a multinomial with size \eqn{m_j}. For \code{bpclm} this is
#'     the full posterior predictive (sample within each posterior draw
#'     of \eqn{\pi}); for \code{pclm} and \code{pclm_exact} it is a
#'     parametric multinomial bootstrap around the point estimate
#'     (\code{n_draws} samples).}
#' }
#'
#' Both flavours \emph{exactly} preserve the wide-bin totals on every
#' draw under \code{type = "predictive"} (the multinomial draws sum to
#' \eqn{m_j} within each band by construction). Under \code{type =
#' "rate"} the totals are preserved exactly only when the input fit has
#' been calibrated (\code{\link{calibrate}}) or is a
#' \code{\link{pclm_exact}} fit.
#'
#' For an ungrouping use case (band counts \eqn{m_j} treated as the
#' data, fine-cell counts treated as the unknown), the recommended
#' workflow is
#' \preformatted{
#'   fit <- bpclm(m, wide_breaks, ...)
#'   fit <- calibrate(fit)                              # exact band totals
#'   pp  <- posterior_predict(fit, type = "predictive") # plausible counts
#' }
#' which gives credible-interval coverage close to nominal for both the
#' fine-cell counts and any cumulative sum of them.
#'
#' @param fit A fitted \code{"pclm"}, \code{"pclm_exact"} or
#'   \code{"bpclm"} object.
#' @param type Either \code{"predictive"} (the default) or \code{"rate"}.
#' @param level Credible/prediction level for the returned interval
#'   (default 0.9).
#' @param n_draws For frequentist input objects (\code{"pclm"},
#'   \code{"pclm_exact"}), the number of parametric multinomial draws
#'   to generate. Ignored for \code{bpclm} input. Default 2000.
#' @param seed Optional integer for reproducibility.
#' @param ... Currently unused; for future extension.
#'
#' @return An object of class \code{"pclm_posterior_predict"}: a list
#'   with components
#'   \describe{
#'     \item{\code{draws}}{An \code{nsim x ngrid} matrix of posterior
#'       (or parametric) draws of the per-cell counts.}
#'     \item{\code{mean, median, lower, upper}}{Per-cell summaries.}
#'     \item{\code{level, type}}{The arguments used.}
#'     \item{\code{grid, grid_mid}}{Fine-grid breakpoints and
#'       midpoints, copied from the input fit.}
#'     \item{\code{wide_breaks, m}}{The wide-bin definitions and
#'       observed counts.}
#'   }
#'
#' @examples
#' \donttest{
#' # Bayesian workflow with exact preservation and uncertainty:
#' data(bloodlead)
#' fit <- bpclm(m = bloodlead$count,
#'              wide_breaks = with(bloodlead, cbind(lower, upper)),
#'              a = 0, b = 80, ngrid = 80, ndx = 17,
#'              niter = 2000, burnin = 500, adapt = 300, seed = 1)
#' fit <- calibrate(fit)                                # exact band totals
#' pp  <- posterior_predict(fit, type = "predictive")
#' plot(pp)
#' }
#'
#' @seealso \code{\link{calibrate}}, \code{\link{pclm_exact}}.
#' @name posterior_predict
#' @rdname posterior_predict
#' @export
posterior_predict <- function(fit, ...) UseMethod("posterior_predict")


#' @rdname posterior_predict
#' @export
posterior_predict.bpclm <- function(fit,
                                    type = c("predictive", "rate"),
                                    level = 0.9,
                                    seed = NULL,
                                    ...) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  M  <- nrow(fit$pi_chain)
  I  <- ncol(fit$pi_chain)
  N  <- sum(fit$m)

  if (type == "rate") {
    draws <- fit$pi_chain * N
  } else {
    assign <- .bin_assignment(fit$C)
    if (any(is.na(assign))) {
      stop("posterior_predict(type = 'predictive') requires the wide bins ",
           "to partition the support without overlap.  For overlapping or ",
           "misaligned bins, use type = 'rate'.")
    }
    draws <- matrix(NA_real_, nrow = M, ncol = I)
    for (s in seq_len(M)) {
      pi_s <- fit$pi_chain[s, ]
      for (j in seq_along(fit$m)) {
        idx <- which(assign == j)
        g_j <- sum(pi_s[idx])
        if (g_j > 0 && fit$m[j] > 0) {
          draws[s, idx] <- rmultinom(1, size = fit$m[j],
                                       prob = pi_s[idx] / g_j)[, 1L]
        } else {
          draws[s, idx] <- 0
        }
      }
    }
  }

  .package_pp(draws, fit, type = type, level = level)
}


#' @rdname posterior_predict
#' @export
posterior_predict.pclm <- function(fit,
                                    type = c("predictive", "rate"),
                                    level = 0.9,
                                    n_draws = 2000L,
                                    seed = NULL,
                                    ...) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  I <- length(fit$pi)
  N <- sum(fit$m)

  if (type == "rate") {
    # Frequentist point estimate -- single-row "chain" with no uncertainty
    rate  <- fit$pi * N
    draws <- matrix(rate, nrow = 1L, ncol = I)
    return(.package_pp(draws, fit, type = type, level = level))
  }

  # Parametric predictive: multinomial bootstrap around the point estimate
  assign <- .bin_assignment(fit$C)
  if (any(is.na(assign))) {
    stop("posterior_predict(type = 'predictive') requires the wide bins ",
         "to partition the support without overlap.  For overlapping or ",
         "misaligned bins, use type = 'rate'.")
  }
  draws <- matrix(NA_real_, nrow = n_draws, ncol = I)
  for (s in seq_len(n_draws)) {
    for (j in seq_along(fit$m)) {
      idx <- which(assign == j)
      g_j <- sum(fit$pi[idx])
      if (g_j > 0 && fit$m[j] > 0) {
        draws[s, idx] <- rmultinom(1, size = fit$m[j],
                                     prob = fit$pi[idx] / g_j)[, 1L]
      } else {
        draws[s, idx] <- 0
      }
    }
  }
  .package_pp(draws, fit, type = type, level = level)
}


# Internal: package the predictive matrix into a "pclm_posterior_predict"
# object.
.package_pp <- function(draws, fit, type, level) {
  alpha <- (1 - level) / 2
  out <- list(
    draws    = draws,
    mean     = colMeans(draws),
    median   = apply(draws, 2L, median),
    lower    = apply(draws, 2L, quantile, probs = alpha,
                       names = FALSE),
    upper    = apply(draws, 2L, quantile, probs = 1 - alpha,
                       names = FALSE),
    level    = level,
    type     = type,
    grid     = fit$grid,
    grid_mid = fit$grid_mid,
    wide_breaks = fit$wide_breaks,
    m        = fit$m
  )
  class(out) <- "pclm_posterior_predict"
  out
}


#' @rdname posterior_predict
#' @export
print.pclm_posterior_predict <- function(x, n = 6L, ...) {
  cat("Posterior",
      if (x$type == "predictive") "predictive (multinomial)"
      else                         "rate (smooth latent)",
      "draws over", length(x$grid_mid), "fine cells\n")
  cat("Number of draws:", nrow(x$draws), "  |  level:", x$level, "\n")
  cat("Wide bins:", length(x$m), "  |  total counts:", sum(x$m), "\n\n")
  head_idx <- seq_len(min(n, length(x$mean)))
  tab <- data.frame(cell  = head_idx,
                    grid_mid = round(x$grid_mid[head_idx], 3L),
                    mean   = round(x$mean[head_idx],  2L),
                    lower  = round(x$lower[head_idx], 2L),
                    upper  = round(x$upper[head_idx], 2L))
  print(tab, row.names = FALSE)
  if (length(x$mean) > n) cat("... (", length(x$mean) - n, " more cells)\n",
                              sep = "")
  invisible(x)
}


#' @rdname posterior_predict
#' @export
summary.pclm_posterior_predict <- function(object, ...) {
  cat("Posterior",
      if (object$type == "predictive") "predictive (multinomial)"
      else                              "rate (smooth latent)",
      "summary\n")
  cat("Total expected count (sum of means): ",
      round(sum(object$mean), 2L), "\n")
  cat("Wide-bin totals reproduced exactly: ",
      isTRUE(all.equal(.collapse_to_bins(object$mean,
                                          object$wide_breaks,
                                          object$grid),
                        object$m, tolerance = 1e-6)), "\n")
  cat("Per-cell coefficient of variation (median): ",
      round(median(ifelse(object$mean > 0,
                          (object$upper - object$lower) /
                            (2 * 1.645 * object$mean),
                          NA), na.rm = TRUE), 3L), "\n")
  invisible(object)
}


# Internal: collapse a fine-cell vector to wide-bin totals (assumes
# bins partition the support exactly).
.collapse_to_bins <- function(x, wide_breaks, grid) {
  J <- nrow(wide_breaks)
  out <- numeric(J)
  fine_lo <- head(grid, -1L)
  fine_hi <- tail(grid, -1L)
  for (j in seq_len(J)) {
    idx <- fine_lo >= wide_breaks[j, 1L] - 1e-12 &
           fine_hi <= wide_breaks[j, 2L] + 1e-12
    out[j] <- sum(x[idx])
  }
  out
}


#' @param show_bins Logical: overlay the observed wide-bin histogram on
#'   the \code{plot()} method.
#' @rdname posterior_predict
#' @export
plot.pclm_posterior_predict <- function(x,
                                        show_bins = TRUE,
                                        xlab = "Fine-grid value",
                                        ylab = "Counts per cell",
                                        main = NULL,
                                        xlim = NULL, ylim = NULL,
                                        lwd = 2, ...) {
  if (is.null(main))
    main <- if (x$type == "predictive")
              sprintf("Posterior predictive: mean and %d%% PI",
                      round(100 * x$level))
            else
              sprintf("Posterior rate: mean and %d%% CI",
                      round(100 * x$level))
  if (is.null(xlim)) xlim <- range(x$grid)
  if (is.null(ylim)) ylim <- c(0, max(x$upper) * 1.05)

  bin_dens <- x$m / (x$wide_breaks[, 2L] - x$wide_breaks[, 1L])

  plot(NA, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, main = main, ...)
  if (show_bins) {
    for (j in seq_len(nrow(x$wide_breaks)))
      rect(x$wide_breaks[j, 1L], 0,
           x$wide_breaks[j, 2L], bin_dens[j],
           col = grDevices::adjustcolor("grey80", 0.5),
           border = "grey50")
  }
  polygon(c(x$grid_mid, rev(x$grid_mid)),
          c(x$lower, rev(x$upper)),
          col = grDevices::adjustcolor("steelblue", 0.3), border = NA)
  lines(x$grid_mid, x$mean, col = "black", lwd = lwd)
  invisible(list(x = x$grid_mid, mean = x$mean,
                 lower = x$lower, upper = x$upper))
}
