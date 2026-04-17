# tests/testthat/test-lam_yao.R

test_that(".lam_yao returns a list with correct elements", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_true(is.list(result))
  expect_named(result, c("k", "ratios", "eigenvalues"))
})

test_that("k is a single positive integer within range", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_length(result$k, 1)
  expect_true(result$k >= 1 && result$k <= 8)
})

test_that("ratios vector has length kmax", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_length(result$ratios, 8)
})

test_that("eigenvalues vector has length kmax + 1", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_length(result$eigenvalues, 9)
})

test_that("eigenvalues of M are in descending order", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_true(all(diff(result$eigenvalues) <= 0))
})

test_that("eigenvalues of M are non-negative", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_true(all(result$eigenvalues >= 0))
})

test_that("ratios are positive", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_true(all(result$ratios > 0))
})

test_that("k is the argmax of ratios", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_equal(result$k, which.max(result$ratios))
})

test_that("h = 1 and h = 3 both run without error", {
  set.seed(42)
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X <- .prepare_matrix(X, standardize = FALSE)
  expect_no_error(.lam_yao(X, kmax = 8, h = 1))
  expect_no_error(.lam_yao(X, kmax = 8, h = 3))
})

test_that("larger h gives different result than h = 1", {
  set.seed(42)
  X       <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1)
  X       <- .prepare_matrix(X, standardize = FALSE)
  result1 <- .lam_yao(X, kmax = 8, h = 1)
  result5 <- .lam_yao(X, kmax = 8, h = 5)
  expect_false(identical(result1$eigenvalues, result5$eigenvalues))
})

test_that("recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 1)
  expect_equal(result$k, k_true)
})

test_that("recovers true number of factors — small N large TT", {
  set.seed(42)
  N <- 25; TT <- 200; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 5, h = 1)
  expect_equal(result$k, k_true)
})

test_that("recovers true number of factors — large N small TT", {
  set.seed(42)
  N <- 200; TT <- 25; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 5, h = 1)
  expect_equal(result$k, k_true)
})

test_that("works with N equal to TT", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 100, k = 3, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  expect_no_error(.lam_yao(X, kmax = 8, h = 1))
})

test_that("works with multiple lags h = 5", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  result <- .lam_yao(X, kmax = 8, h = 5)
  expect_true(result$k >= 1 && result$k <= 8)
})
