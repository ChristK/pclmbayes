# bspline_basis dimensions and partition-of-unity
x <- seq(0, 10, length.out = 50)
bs <- bspline_basis(x, a = 0, b = 10, ndx = 7L, degree = 3L)
expect_inherits(bs, "pclm_basis")
expect_equal(bs$K, 7L + 3L)
expect_equal(dim(bs$B), c(50L, 10L))
# B-splines of given degree on equally-spaced knots form a partition of
# unity at any interior point of the support.
row_sums <- rowSums(bs$B)
expect_true(all(abs(row_sums - 1) < 1e-8))

# diff_penalty has the expected null space and rank
P3 <- diff_penalty(K = 8L, r = 3L)
expect_equal(dim(P3), c(8L, 8L))
expect_equal(qr(P3)$rank, 8L - 3L)
# r = 3 => null space contains polynomials up to degree 2
k <- 1:8
expect_true(all(abs(P3 %*% rep(1, 8)) < 1e-10))
expect_true(all(abs(P3 %*% k)        < 1e-10))
expect_true(all(abs(P3 %*% k^2)      < 1e-10))

# bin_matrix gives 0/1 weights when bins align with the fine grid
fine_breaks <- seq(0, 10, by = 1)
wide_breaks <- c(0, 2, 5, 10)
C <- bin_matrix(wide_breaks, fine_breaks)
expect_equal(dim(C), c(3L, 10L))
expect_true(all(C %in% c(0, 1)))
expect_true(all(colSums(C) == 1))
expect_equal(rowSums(C), c(2, 3, 5))
pi_uni <- rep(1 / 10, 10)
expect_equal(as.numeric(C %*% pi_uni), c(0.2, 0.3, 0.5),
             tolerance = 1e-12)

# bin_matrix uses partial overlap when bins misalign
fine_breaks <- seq(0, 1, length.out = 11)
C <- bin_matrix(matrix(c(0.05, 0.55), nrow = 1), fine_breaks)
expect_equal(dim(C), c(1L, 10L))
expect_equal(unname(C[1, 1]), 0.5,  tolerance = 1e-12)
expect_equal(unname(sum(C[1, 2:5])), 4, tolerance = 1e-12)
expect_equal(unname(C[1, 6]), 0.5,  tolerance = 1e-12)
expect_equal(unname(sum(C[1, 7:10])), 0)

# latent_density returns a valid discrete pdf summing to 1
x <- seq(0, 10, length.out = 51)
mids <- (head(x, -1) + tail(x, -1)) / 2
bs <- bspline_basis(mids, a = 0, b = 10, ndx = 7L, degree = 3L)
set.seed(42)
phi <- rnorm(bs$K, sd = 0.5)
pi  <- latent_density(phi, bs)
expect_equal(length(pi), 50L)
expect_true(all(pi >= 0))
expect_equal(sum(pi), 1, tolerance = 1e-12)
