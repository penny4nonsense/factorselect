# tests/testthat/test-onatski_2010.R

test_that(".onatski_2010 returns a list with correct elements", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_true(is.list(result))
  expect_named(result, c("k", "delta", "beta", "differences", "n_iter"))
})

test_that("k is a non-negative integer within range", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_length(result$k, 1)
  expect_true(result$k >= 0 && result$k <= 8)
})

test_that("delta is positive", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_true(result$delta > 0)
})

test_that("delta equals 2 * abs(beta)", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_equal(result$delta, 2 * abs(result$beta))
})

test_that("differences vector has length kmax", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_length(result$differences, 8)
})

test_that("differences are positive", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_true(all(result$differences >= 0))
})

test_that("n_iter is returned correctly", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8, n_iter = 4)
  expect_equal(result$n_iter, 4)
})

test_that("insufficient eigenvalues throws informative error", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  expect_error(.onatski_2010(eig$values, kmax = 8), "kmax \\+ 5")
})

test_that("different n_iter values run without error", {
  set.seed(42)
  X   <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X   <- .prepare_matrix(X)
  eig <- .extract_eigenvalues(X, kmax = 13)
  expect_no_error(.onatski_2010(eig$values, kmax = 8, n_iter = 1))
  expect_no_error(.onatski_2010(eig$values, kmax = 8, n_iter = 4))
  expect_no_error(.onatski_2010(eig$values, kmax = 8, n_iter = 10))
})

test_that("more iterations does not break the result", {
  set.seed(42)
  X       <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X       <- .prepare_matrix(X)
  eig     <- .extract_eigenvalues(X, kmax = 13)
  result2 <- .onatski_2010(eig$values, kmax = 8, n_iter = 2)
  result4 <- .onatski_2010(eig$values, kmax = 8, n_iter = 4)
  # Both should be valid estimates
  expect_true(result2$k >= 0 && result2$k <= 8)
  expect_true(result4$k >= 0 && result4$k <= 8)
})

test_that("recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 13)
  result <- .onatski_2010(eig$values, kmax = 8)
  expect_equal(result$k, k_true)
})

test_that("recovers true number of factors — small N large TT", {
  set.seed(42)
  N <- 25; TT <- 200; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 10)
  result <- .onatski_2010(eig$values, kmax = 5)
  expect_equal(result$k, k_true)
})

test_that("recovers true number of factors — large N small TT", {
  set.seed(42)
  N <- 200; TT <- 50; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 10)
  result <- .onatski_2010(eig$values, kmax = 5)
  expect_equal(result$k, k_true)
})

test_that("works with N equal to TT", {
  set.seed(42)
  X   <- simulate_factor_model(N = 100, TT = 100, k = 3, sd = 0.5)
  X   <- .prepare_matrix(X)
  eig <- .extract_eigenvalues(X, kmax = 13)
  expect_no_error(.onatski_2010(eig$values, kmax = 8))
})

test_that("works with N > TT", {
  set.seed(42)
  X   <- simulate_factor_model(N = 200, TT = 50, k = 2, sd = 0.5)
  X   <- .prepare_matrix(X)
  eig <- .extract_eigenvalues(X, kmax = 10)
  expect_no_error(.onatski_2010(eig$values, kmax = 5))
})
