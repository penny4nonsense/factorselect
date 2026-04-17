# tests/testthat/test-ahn_horenstein.R

test_that(".ahn_horenstein returns a list with correct elements", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_true(is.list(result))
  expect_named(result, c("k_er", "k_gr", "er", "gr"))
})

test_that("k_er is a single positive integer within range", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_length(result$k_er, 1)
  expect_true(result$k_er >= 1 && result$k_er <= 5)
})

test_that("k_gr is a single positive integer within range", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_length(result$k_gr, 1)
  expect_true(result$k_gr >= 1 && result$k_gr <= 5)
})

test_that("er statistic vector has length kmax", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_length(result$er, 5)
})

test_that("gr statistic vector has length kmax", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_length(result$gr, 5)
})

test_that("er statistic values are positive", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_true(all(result$er > 0))
})

test_that("gr statistic values are positive", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_true(all(result$gr > 0))
})

test_that("k_er is the argmax of er", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_equal(result$k_er, which.max(result$er))
})

test_that("k_gr is the argmax of gr", {
  set.seed(42)
  X      <- matrix(rnorm(200), 20, 10)
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(nrow(X), ncol(X)))
  expect_equal(result$k_gr, which.max(result$gr))
})

test_that("ER recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; T <- 200; k_true <- 3
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .ahn_horenstein(eig$values, kmax = 8, n = min(N, T))
  expect_equal(result$k_er, k_true)
})

test_that("GR recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; T <- 200; k_true <- 3
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .ahn_horenstein(eig$values, kmax = 8, n = min(N, T))
  expect_equal(result$k_gr, k_true)
})

test_that("ER recovers true number of factors — small N large T", {
  set.seed(42)
  N <- 25; T <- 200; k_true <- 2
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .ahn_horenstein(eig$values, kmax = 8, n = min(N, T))
  expect_equal(result$k_er, k_true)
})

test_that("GR recovers true number of factors — large N small T", {
  set.seed(42)
  N <- 200; T <- 10; k_true <- 2
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  X      <- .prepare_matrix(X)
  eig    <- .extract_eigenvalues(X, kmax = 5)
  result <- .ahn_horenstein(eig$values, kmax = 5, n = min(N, T))
  expect_equal(result$k_gr, k_true)
})
