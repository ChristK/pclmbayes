# =============================================================================
# Tests for S3 methods on "pclm" and "bpclm" objects
# =============================================================================

# ---- Set up small fits to reuse -----------------------------------------
set.seed(2026)
data(bloodlead, envir = environment(), package = "pclmbayes")

fit_p <- pclm(m = bloodlead$count,
              wide_breaks = with(bloodlead, cbind(lower, upper)),
              a = 0, b = 80, ngrid = 60,
              ndx = 13L, degree = 3L, penalty_order = 3L)

fit_b <- bpclm(m = bloodlead$count,
               wide_breaks = with(bloodlead, cbind(lower, upper)),
               a = 0, b = 80, ngrid = 60,
               ndx = 13L, degree = 3L, penalty_order = 3L,
               niter = 500L, burnin = 100L, adapt = 150L,
               seed = 11)

# ---- coef(), fitted(), logLik() on pclm --------------------------------
co <- coef(fit_p)
expect_true(is.numeric(co))
expect_equal(length(co), fit_p$basis$K)
expect_equal(sum(co), 0, tolerance = 1e-8)

fv <- fitted(fit_p)
expect_true(is.numeric(fv))
expect_equal(length(fv), length(fit_p$m))
expect_true(all(fv >= 0))

ll <- logLik(fit_p)
expect_inherits(ll, "logLik")
expect_equal(attr(ll, "nobs"), sum(fit_p$m))
expect_true(is.finite(as.numeric(ll)))

# ---- predict.pclm with newdata -----------------------------------------
nd <- c(-5, 0, 10, 40, 80, 100)
pr <- predict(fit_p, newdata = nd)
expect_equal(length(pr), length(nd))
expect_true(all(pr >= 0))
# Outside [a, b] => 0
expect_equal(pr[1], 0)
expect_equal(pr[6], 0)

# Default (no newdata) returns density on the fine grid
pr0 <- predict(fit_p)
expect_equal(length(pr0), length(fit_p$grid_mid))

# ---- predict.bpclm: all summary modes, with and without newdata --------
# Default = "mean"
d_mean <- predict(fit_b)
expect_equal(length(d_mean), length(fit_b$grid_mid))
expect_true(all(d_mean >= 0))

d_med <- predict(fit_b, summary = "median")
expect_equal(length(d_med), length(fit_b$grid_mid))

d_samp <- predict(fit_b, summary = "sample")
expect_true(is.matrix(d_samp))
expect_equal(dim(d_samp),
             c(length(fit_b$grid_mid), nrow(fit_b$pi_chain)))

# With newdata
nd <- c(-1, 5, 25, 50, 80, 90)
pn_mean <- predict(fit_b, newdata = nd)
expect_equal(length(pn_mean), length(nd))
expect_equal(pn_mean[1], 0)
expect_equal(pn_mean[6], 0)

pn_med <- predict(fit_b, newdata = nd, summary = "median")
expect_equal(length(pn_med), length(nd))

pn_samp <- predict(fit_b, newdata = nd, summary = "sample")
expect_true(is.matrix(pn_samp))
expect_equal(dim(pn_samp), c(length(nd), nrow(fit_b$pi_chain)))

# ---- quantile.bpclm summary = "sample" ---------------------------------
qs <- quantile(fit_b, probs = c(0.25, 0.5, 0.75), summary = "sample")
expect_true(is.matrix(qs))
expect_equal(dim(qs), c(nrow(fit_b$pi_chain), 3L))
expect_equal(colnames(qs), c("25%", "50%", "75%"))

qm <- quantile(fit_b, probs = c(0.25, 0.5, 0.75), summary = "median")
expect_true(is.data.frame(qm))
expect_equal(nrow(qm), 3L)

# ---- summary methods exercise ------------------------------------------
out_p <- capture.output(summary(fit_p))
expect_true(length(out_p) > 0L)

out_b <- capture.output(s_b <- summary(fit_b))
expect_true(length(out_b) > 0L)
expect_true(is.list(s_b))

# fitted.bpclm returns expected counts
fb <- fitted(fit_b)
expect_equal(length(fb), length(fit_b$m))
expect_true(all(fb >= 0))
