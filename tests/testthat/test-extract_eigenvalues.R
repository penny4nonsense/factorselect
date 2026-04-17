# tests/testthat/test-extract_eigenvalues.R

test_that(".extract_eigenvalues returns a list with correct elements", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_true(is.list(result))
  expect_named(result, c("values", "vectors", "orientation"))
})

test_that("eigenvalues are returned in descending order", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_true(all(diff(result$values) <= 0))
})

test_that("correct number of eigenvalues returned (kmax + 1)", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_equal(length(result$values), 6)
})

test_that("correct number of eigenvectors returned (kmax + 1)", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_equal(ncol(result$vectors), 6)
})

test_that("eigenvalues are positive", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_true(all(result$values > 0))
})

test_that("orientation is N when N <= T", {
  # T = 20, N = 10, so N <= T
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_equal(result$orientation, "N")
})

test_that("orientation is T when N > T", {
  # T = 10, N = 20, so N > T
  X <- matrix(rnorm(200), 10, 20)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_equal(result$orientation, "T")
})

test_that("eigenvector matrix has correct number of rows when orientation is N", {
  # N x N decomposition so vectors should have N rows = ncol(X)
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_equal(nrow(result$vectors), ncol(X))  # N = ncol(X)
})

test_that("eigenvector matrix has correct number of rows when orientation is T", {
  # T x T decomposition so vectors should have T rows = nrow(X)
  X <- matrix(rnorm(200), 10, 20)
  X <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 5)
  expect_equal(nrow(result$vectors), nrow(X))  # T = nrow(X)
})

test_that("kmax too large generates warning and reduces kmax", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  expect_warning(
    .extract_eigenvalues(X, kmax = 50),
    "kmax"
  )
})

test_that("reduced kmax still returns correct structure", {
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)
  result <- suppressWarnings(.extract_eigenvalues(X, kmax = 50))
  expect_true(is.list(result))
  expect_named(result, c("values", "vectors", "orientation"))
})

test_that("eigenvalues match base R eigen on small matrix", {
  set.seed(42)
  X <- matrix(rnorm(200), 20, 10)
  X <- .prepare_matrix(X)

  # Our function
  result <- .extract_eigenvalues(X, kmax = 5)

  # Base R directly
  S    <- tcrossprod(X) / nrow(X)
  expected <- eigen(S, symmetric = TRUE)$values[1:6]

  expect_equal(result$values, expected, tolerance = 1e-10)
})

test_that("eigenvalues recover true factor structure", {
  set.seed(42)
  N <- 100; T <- 200; k_true <- 3

  # Simulate factor model with 3 strong factors
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E

  X      <- .prepare_matrix(X)
  result <- .extract_eigenvalues(X, kmax = 8)

  # First 3 eigenvalues should be much larger than the rest
  gap <- result$values[k_true] / result$values[k_true + 1]
  expect_gt(gap, 5)
})
