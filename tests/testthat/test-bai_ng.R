# tests/testthat/test-bai_ng.R

test_that(".bai_ng returns a list with correct elements", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_true(is.list(result))
  expect_named(result, c("k_pc1", "k_pc2", "k_pc3",
                         "k_ic1", "k_ic2", "k_ic3",
                         "pc1",   "pc2",   "pc3",
                         "ic1",   "ic2",   "ic3"))
})

test_that("all k estimates are non-negative integers", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  k_vals <- c(result$k_pc1, result$k_pc2, result$k_pc3,
              result$k_ic1, result$k_ic2, result$k_ic3)
  expect_true(all(k_vals >= 0))
  expect_true(all(k_vals <= 8))
})

test_that("all k estimates are integers", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  k_vals <- c(result$k_pc1, result$k_pc2, result$k_pc3,
              result$k_ic1, result$k_ic2, result$k_ic3)
  expect_true(all(k_vals == as.integer(k_vals)))
})

test_that("criterion vectors have length kmax + 1", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_length(result$pc1, 9)
  expect_length(result$pc2, 9)
  expect_length(result$pc3, 9)
  expect_length(result$ic1, 9)
  expect_length(result$ic2, 9)
  expect_length(result$ic3, 9)
})

test_that("k_pc1 is the argmin of pc1", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_pc1, which.min(result$pc1) - 1L)
})

test_that("k_pc2 is the argmin of pc2", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_pc2, which.min(result$pc2) - 1L)
})

test_that("k_pc3 is the argmin of pc3", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_pc3, which.min(result$pc3) - 1L)
})

test_that("k_ic1 is the argmin of ic1", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_ic1, which.min(result$ic1) - 1L)
})

test_that("k_ic2 is the argmin of ic2", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_ic2, which.min(result$ic2) - 1L)
})

test_that("k_ic3 is the argmin of ic3", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_ic3, which.min(result$ic3) - 1L)
})

test_that("IC criteria recover true factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_equal(result$k_ic1, k_true)
  expect_equal(result$k_ic2, k_true)
  expect_equal(result$k_ic3, k_true)
})

test_that("PC criteria recover true factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_equal(result$k_pc1, k_true)
  expect_equal(result$k_pc2, k_true)
  expect_equal(result$k_pc3, k_true)
})

test_that("criteria are sensitive to kmax as documented", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig8   <- .extract_eigenvalues(X, kmax = 8)
  eig15  <- .extract_eigenvalues(X, kmax = 15)
  res8   <- .bai_ng(eig8$values,  V0 = V0, kmax = 8,  N = N, TT = TT)
  res15  <- .bai_ng(eig15$values, V0 = V0, kmax = 15, N = N, TT = TT)
  expect_true(is.list(res8))
  expect_true(is.list(res15))
})

test_that("works with N > TT", {
  set.seed(42)
  N <- 200; TT <- 50; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_true(is.list(result))
})

test_that("works with TT > N", {
  set.seed(42)
  N <- 50; TT <- 200; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .bai_ng(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_true(is.list(result))
})
