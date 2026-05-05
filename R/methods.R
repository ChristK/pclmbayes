# =============================================================================
# S3 methods for "pclm" and "bpclm" objects
# -----------------------------------------------------------------------------
# Convention: pi_i = f_Y(u_i) * Delta on a fine grid of width Delta. The PDF
# evaluated at the midpoints u_i is therefore pi_i / Delta.
# =============================================================================

# Internal: convert a vector of grid probabilities pi to a density f(u) at the
# midpoints, given the fine-grid breakpoints `grid` (length I + 1).
.pi_to_density <- function(pi, grid) {
  delta <- diff(grid)
  pi / delta
}

# Internal: linear interpolation of the discrete density pi at user-supplied
# evaluation points x. Returns 0 outside the support [a, b].
.density_at <- function(x, mids, pi, grid) {
  a <- grid[1L]
  b <- grid[length(grid)]
  dens <- .pi_to_density(pi, grid)  # length(mids), evaluated at midpoints
  # Inside [a, b] but possibly slightly outside the range of midpoints:
  # extrapolate from the nearest midpoint (rule = 2). Outside [a, b]: 0.
  out <- approx(x = mids, y = dens, xout = x,
                method = "linear", rule = 2)$y
  out[x < a | x > b] <- 0
  out
}

# Internal: cumulative distribution function on the fine grid (vector of
# length I + 1) implied by pi.
.cdf_grid <- function(pi) {
  c(0, cumsum(pi))
}

# Internal: invert the CDF to compute a quantile of probability `p`.
.quantile_from_pi <- function(p, pi, grid) {
  Fg <- .cdf_grid(pi)        # length I + 1, monotonic in [0, 1]
  # Use approx with x = Fg, y = grid breakpoints to interpolate
  approx(x = Fg, y = grid, xout = p, ties = mean, rule = 2)$y
}

# Internal: posterior mean / sd of Y given pi
.moments_from_pi <- function(pi, mids) {
  mu  <- sum(mids * pi)
  s2  <- sum((mids - mu) ^ 2 * pi)
  c(mean = mu, sd = sqrt(s2))
}

# -----------------------------------------------------------------------------
# Methods for "pclm"
# -----------------------------------------------------------------------------

#' @export
print.pclm <- function(x, digits = 4L, ...) {
  cat("Penalised composite link model (frequentist)\n")
  cat("Call: "); print(x$call)
  cat("\nNumber of wide bins:", length(x$m),
      " | total counts:", sum(x$m), "\n")
  cat("Fine grid:", length(x$grid_mid), "intervals on (",
      signif(x$grid[1L], digits), ", ",
      signif(x$grid[length(x$grid)], digits), ")\n", sep = "")
  cat("B-spline basis: K =", x$basis$K, "(degree =", x$basis$degree, ")\n")
  cat("Penalty order r =", x$penalty_order, "\n")
  cat(sprintf("Selected tau = %.4g (%s = %.2f, edf = %.2f)\n",
              x$tau, x$select, x$ic, x$edf))
  cat(sprintf("Converged = %s in %d iterations.  log-likelihood = %.3f\n",
              isTRUE(x$converged), x$iter, x$logL))
  invisible(x)
}

#' @export
summary.pclm <- function(object, probs = c(0.05, 0.10, 0.25, 0.50,
                                           0.75, 0.90, 0.95), ...) {
  print(object)
  cat("\nFitted summary statistics of the latent density:\n")
  mom <- .moments_from_pi(object$pi, object$grid_mid)
  cat(sprintf("  mean = %.4f, sd = %.4f\n", mom["mean"], mom["sd"]))
  qs <- .quantile_from_pi(probs, object$pi, object$grid)
  qtab <- setNames(qs, paste0(format(probs * 100), "%"))
  print(round(qtab, 4L))
  cat("\nGoodness of fit (observed vs fitted counts):\n")
  cmp <- data.frame(lower = object$wide_breaks[, 1L],
                    upper = object$wide_breaks[, 2L],
                    obs   = object$m,
                    exp   = round(object$fitted_counts, 2L))
  print(cmp, row.names = FALSE)
  invisible(object)
}

#' @export
coef.pclm <- function(object, ...) object$phi

#' @export
fitted.pclm <- function(object, ...) object$fitted_counts

#' @export
logLik.pclm <- function(object, ...) {
  structure(object$logL,
            df    = object$edf,
            nobs  = sum(object$m),
            class = "logLik")
}

#' Plot a frequentist PCLM fit
#'
#' Draws the histogram of wide-bin densities together with the fitted
#' latent density on the fine grid.
#'
#' @param x A \code{"pclm"} object.
#' @param add Logical: if \code{TRUE}, add to an existing plot.
#' @param density_col Colour for the fitted density line.
#' @param hist_col Fill colour for the histogram rectangles.
#' @param hist_border Border colour for histogram rectangles.
#' @param xlab,ylab,main,ylim,xlim,lwd,... Standard graphical
#'   parameters.
#'
#' @return Invisibly, the data underlying the plot (a list with
#'   components \code{x} = grid midpoints and \code{y} = density).
#' @export
plot.pclm <- function(x,
                      add = FALSE,
                      density_col = "black",
                      hist_col = grDevices::adjustcolor("grey70", alpha.f = 0.6),
                      hist_border = "grey40",
                      xlab = "y", ylab = "Density",
                      main = "Frequentist PCLM fit",
                      xlim = NULL, ylim = NULL, lwd = 2, ...) {
  mids <- x$grid_mid
  dens <- .pi_to_density(x$pi, x$grid)
  # Histogram heights from the wide bins
  wb       <- x$wide_breaks
  bin_w    <- wb[, 2L] - wb[, 1L]
  bin_dens <- (x$m / sum(x$m)) / bin_w
  if (is.null(xlim)) xlim <- range(c(wb[, 1L], wb[, 2L], mids))
  if (is.null(ylim)) ylim <- c(0, max(c(dens, bin_dens)) * 1.05)

  if (!add) {
    plot.default(NA, type = "n", xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main, ...)
    # Histogram rectangles
    for (j in seq_along(bin_dens)) {
      rect(xleft = wb[j, 1L], ybottom = 0,
           xright = wb[j, 2L], ytop = bin_dens[j],
           col = hist_col, border = hist_border)
    }
  }
  lines(mids, dens, col = density_col, lwd = lwd)
  invisible(list(x = mids, y = dens))
}

#' @export
predict.pclm <- function(object, newdata, ...) {
  if (missing(newdata)) {
    .pi_to_density(object$pi, object$grid)
  } else {
    .density_at(as.numeric(newdata),
                mids = object$grid_mid,
                pi   = object$pi,
                grid = object$grid)
  }
}

#' Quantile method for PCLM fits
#'
#' Inverts the latent CDF on the fine grid by linear interpolation and
#' returns the requested quantiles.
#'
#' @param x A \code{"pclm"} object.
#' @param probs Probabilities at which to compute quantiles.
#' @param ... Ignored.
#' @return Numeric vector of quantiles (named by \code{probs}).
#' @export
quantile.pclm <- function(x, probs = c(0.25, 0.5, 0.75), ...) {
  q <- .quantile_from_pi(probs, x$pi, x$grid)
  setNames(q, paste0(format(probs * 100), "%"))
}

# -----------------------------------------------------------------------------
# Methods for "bpclm"
# -----------------------------------------------------------------------------

#' @export
print.bpclm <- function(x, digits = 4L, ...) {
  cat("Bayesian penalised composite link model\n")
  cat("Call: "); print(x$call)
  cat("\nNumber of wide bins:", length(x$m),
      " | total counts:", sum(x$m), "\n")
  cat("Fine grid:", length(x$grid_mid), "intervals on (",
      signif(x$grid[1L], digits), ", ",
      signif(x$grid[length(x$grid)], digits), ")\n", sep = "")
  cat("B-spline basis: K =", x$basis$K, "(degree =", x$basis$degree, ")\n")
  cat("Penalty order r =", x$penalty_order, "\n")
  cat(sprintf("MCMC: niter = %d, burnin = %d, thin = %d, kept = %d\n",
              x$niter, x$burnin, x$thin, length(x$tau_chain)))
  cat(sprintf("Final delta = %.3g  |  acceptance rate = %.2f\n",
              x$delta, x$accept))
  if (!is.null(x$shape) && length(x$shape))
    cat("Shape constraint(s):", paste(x$shape, collapse = ", "), "\n")
  cat(sprintf("Posterior mean tau = %.3g (sd %.3g)\n",
              mean(x$tau_chain), sd(x$tau_chain)))
  invisible(x)
}

#' @export
summary.bpclm <- function(object,
                          probs = c(0.05, 0.10, 0.25, 0.50,
                                    0.75, 0.90, 0.95),
                          cred = 0.90, ...) {
  print(object)

  alpha <- (1 - cred) / 2
  # Posterior summaries of moments
  mom_chain <- t(apply(object$pi_chain, 1L, .moments_from_pi,
                       mids = object$grid_mid))
  mu_summ <- c(mean = mean(mom_chain[, 1L]),
               lo   = quantile(mom_chain[, 1L], alpha,     names = FALSE),
               hi   = quantile(mom_chain[, 1L], 1 - alpha, names = FALSE))
  sd_summ <- c(mean = mean(mom_chain[, 2L]),
               lo   = quantile(mom_chain[, 2L], alpha,     names = FALSE),
               hi   = quantile(mom_chain[, 2L], 1 - alpha, names = FALSE))

  cat(sprintf("\nPosterior of mean(Y):  %.4f  (%.0f%% CI: %.4f, %.4f)\n",
              mu_summ["mean"], 100 * cred, mu_summ["lo"], mu_summ["hi"]))
  cat(sprintf("Posterior of sd(Y):    %.4f  (%.0f%% CI: %.4f, %.4f)\n",
              sd_summ["mean"], 100 * cred, sd_summ["lo"], sd_summ["hi"]))

  # Posterior of selected quantiles
  q_chain <- t(apply(object$pi_chain, 1L,
                     function(pp) .quantile_from_pi(probs, pp, object$grid)))
  q_mean  <- colMeans(q_chain)
  q_lo    <- apply(q_chain, 2L, quantile, probs = alpha,     names = FALSE)
  q_hi    <- apply(q_chain, 2L, quantile, probs = 1 - alpha, names = FALSE)
  cat(sprintf("\nPosterior summaries of quantiles (mean and %.0f%% CI):\n",
              100 * cred))
  qtab <- data.frame(p    = probs,
                     mean = round(q_mean, 4L),
                     lo   = round(q_lo,   4L),
                     hi   = round(q_hi,   4L))
  print(qtab, row.names = FALSE)
  invisible(list(mu = mu_summ, sd = sd_summ, quantiles = qtab))
}

#' @export
coef.bpclm <- function(object, ...) object$phi

#' @export
fitted.bpclm <- function(object, ...) {
  # Fitted bin counts under posterior mean of pi
  ga <- as.numeric(object$C %*% object$pi)
  sum(object$m) * ga
}

#' Plot a Bayesian PCLM fit
#'
#' Histogram of wide-bin densities, the posterior mean of the latent
#' density, and a pointwise credible band.
#'
#' @inheritParams plot.pclm
#' @param x A \code{"bpclm"} object.
#' @param band_col Fill colour for the credible band (alpha-blended).
#' @param cred Credible level for the band (default
#'   \code{x$cred_level}).
#' @return Invisibly, a list with components \code{x}, \code{y},
#'   \code{lower}, \code{upper}.
#' @export
plot.bpclm <- function(x,
                       add = FALSE,
                       cred = NULL,
                       density_col = "black",
                       band_col = grDevices::adjustcolor("steelblue",
                                                          alpha.f = 0.25),
                       hist_col = grDevices::adjustcolor("grey70",
                                                          alpha.f = 0.6),
                       hist_border = "grey40",
                       xlab = "y", ylab = "Density",
                       main = "Bayesian PCLM fit (posterior mean + credible band)",
                       xlim = NULL, ylim = NULL, lwd = 2, ...) {
  if (is.null(cred)) cred <- x$cred_level
  mids   <- x$grid_mid
  delta  <- diff(x$grid)
  dens   <- x$pi / delta
  if (!isTRUE(all.equal(cred, x$cred_level))) {
    alpha <- (1 - cred) / 2
    lo <- apply(x$pi_chain, 2L, quantile, probs = alpha,
                names = FALSE) / delta
    hi <- apply(x$pi_chain, 2L, quantile, probs = 1 - alpha,
                names = FALSE) / delta
  } else {
    lo <- x$pi_lower / delta
    hi <- x$pi_upper / delta
  }
  # Histogram
  wb       <- x$wide_breaks
  bin_w    <- wb[, 2L] - wb[, 1L]
  bin_dens <- (x$m / sum(x$m)) / bin_w
  if (is.null(xlim)) xlim <- range(c(wb[, 1L], wb[, 2L], mids))
  if (is.null(ylim)) ylim <- c(0, max(c(dens, hi, bin_dens)) * 1.05)

  if (!add) {
    plot.default(NA, type = "n", xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main, ...)
    for (j in seq_along(bin_dens)) {
      rect(xleft = wb[j, 1L], ybottom = 0,
           xright = wb[j, 2L], ytop = bin_dens[j],
           col = hist_col, border = hist_border)
    }
  }
  # Credible band
  polygon(c(mids, rev(mids)), c(lo, rev(hi)),
          col = band_col, border = NA)
  lines(mids, dens, col = density_col, lwd = lwd)
  invisible(list(x = mids, y = dens, lower = lo, upper = hi))
}

#' @export
predict.bpclm <- function(object, newdata,
                          summary = c("mean", "median", "sample"),
                          cred = NULL, ...) {
  summary <- match.arg(summary)
  delta  <- diff(object$grid)
  dens   <- object$pi / delta
  if (missing(newdata)) {
    if (summary == "sample") {
      return(t(object$pi_chain) / delta)  # I x nsim
    }
    if (summary == "median") {
      med <- apply(object$pi_chain, 2L, median) / delta
      return(med)
    }
    return(dens)
  }
  x <- as.numeric(newdata)
  if (summary == "sample") {
    M    <- nrow(object$pi_chain)
    out  <- matrix(NA_real_, nrow = length(x), ncol = M)
    for (mi in seq_len(M)) {
      out[, mi] <- .density_at(x, object$grid_mid,
                               object$pi_chain[mi, ], object$grid)
    }
    return(out)
  }
  if (summary == "median") {
    med <- apply(object$pi_chain, 2L, median)
    return(.density_at(x, object$grid_mid, med, object$grid))
  }
  .density_at(x, object$grid_mid, object$pi, object$grid)
}

#' Quantile method for Bayesian PCLM fits
#'
#' Returns posterior summaries of any quantile of the latent
#' distribution. The quantile of the posterior-mean density is returned
#' if \code{summary = "mean"}; with \code{summary = "sample"} the full
#' posterior sample of quantile values is returned.
#'
#' @param x A \code{"bpclm"} object.
#' @param probs Probabilities at which to compute quantiles.
#' @param summary One of \code{"mean"} (posterior mean of the quantile),
#'   \code{"median"}, or \code{"sample"} (returns the full posterior
#'   sample as a matrix).
#' @param cred Credible level (default \code{x$cred_level}). When
#'   \code{summary = "mean"}, also reports lower/upper credible limits.
#' @param ... Ignored.
#' @return If \code{summary = "sample"}: a \code{nsim x length(probs)}
#'   matrix. Otherwise a \code{length(probs) x 4} data frame with
#'   columns \code{p}, \code{estimate}, \code{lower}, \code{upper}.
#' @export
quantile.bpclm <- function(x,
                           probs = c(0.25, 0.5, 0.75),
                           summary = c("mean", "median", "sample"),
                           cred = NULL, ...) {
  summary <- match.arg(summary)
  if (is.null(cred)) cred <- x$cred_level
  alpha <- (1 - cred) / 2
  q_chain <- t(apply(x$pi_chain, 1L,
                     function(pp) .quantile_from_pi(probs, pp, x$grid)))
  if (summary == "sample") {
    colnames(q_chain) <- paste0(format(probs * 100), "%")
    return(q_chain)
  }
  est <- if (summary == "median") apply(q_chain, 2L, median)
         else                     colMeans(q_chain)
  lo  <- apply(q_chain, 2L, quantile, probs = alpha,     names = FALSE)
  hi  <- apply(q_chain, 2L, quantile, probs = 1 - alpha, names = FALSE)
  data.frame(p = probs, estimate = est, lower = lo, upper = hi)
}
