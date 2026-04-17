# tests/testthat/test-onatski_2009.R

test_that(".onatski_2009 returns a list with correct elements", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_true(is.list(result))
  expect_named(result, c("k", "ratios", "eigenvalues", "critical_value", "alpha"))
})

test_that("k is a non-negative integer within range", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_length(result$k, 1)
  expect_true(result$k >= 0 && result$k <= 8)
})

test_that("ratios vector has length kmax", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_length(result$ratios, 8)
})

test_that("eigenvalues vector has length kmax + 2", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_length(result$eigenvalues, 10)
})

test_that("eigenvalues are in descending order", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_true(all(diff(result$eigenvalues) <= 0))
})

test_that("eigenvalues are non-negative", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_true(all(result$eigenvalues >= 0))
})

test_that("alpha is returned correctly", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_equal(result$alpha, 0.05)
})

test_that("critical_value is positive", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_true(result$critical_value > 0)
})

test_that("invalid alpha throws informative error", {
  set.seed(42)
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X <- .prepare_matrix(X)
  expect_error(.onatski_2009(X, kmax = 8, alpha = 0.25), "alpha must be one of")
})

test_that("larger alpha gives smaller critical value", {
  set.seed(42)
  X       <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X       <- .prepare_matrix(X)
  result1 <- .onatski_2009(X, kmax = 8, alpha = 0.01)
  result2 <- .onatski_2009(X, kmax = 8, alpha = 0.10)
  expect_gt(result1$critical_value, result2$critical_value)
})

test_that("odd T is handled by dropping last observation", {
  set.seed(42)
  X <- simulate_factor_model(N = 100, TT = 201, k = 3, sd = 0.5)
  X <- .prepare_matrix(X)
  expect_no_error(.onatski_2009(X, kmax = 8, alpha = 0.05))
})

test_that("all valid alpha levels run without error", {
  set.seed(42)
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
  X <- .prepare_matrix(X)
  alpha_levels <- c(0.15, 0.10, 0.09, 0.08, 0.07, 0.06,
                    0.05, 0.04, 0.03, 0.02, 0.01)
  for (a in alpha_levels) {
    expect_no_error(.onatski_2009(X, kmax = 8, alpha = a))
  }
})

test_that("recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
  expect_equal(result$k, k_true)
})

test_that("recovers true number of factors — small N large TT", {
  set.seed(42)
  N <- 25; TT <- 200; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 5, alpha = 0.05)
  expect_equal(result$k, k_true)
})

test_that("recovers true number of factors — large N small TT", {
  set.seed(42)
  N <- 200; TT <- 50; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X      <- .prepare_matrix(X)
  result <- .onatski_2009(X, kmax = 5, alpha = 0.05)
  expect_equal(result$k, k_true)
})

test_that("stricter alpha produces equal or fewer factors", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X       <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 0.5)
  X       <- .prepare_matrix(X)
  result1 <- .onatski_2009(X, kmax = 8, alpha = 0.01)
  result2 <- .onatski_2009(X, kmax = 8, alpha = 0.15)
  expect_true(result1$k <= result2$k)
})

test_that("works with N equal to TT", {
  set.seed(42)
  X <- simulate_factor_model(N = 100, TT = 100, k = 3, sd = 0.5)
  X <- .prepare_matrix(X)
  expect_no_error(.onatski_2009(X, kmax = 8, alpha = 0.05))
})

test_that("works with N > TT", {
  set.seed(42)
  X <- simulate_factor_model(N = 200, TT = 50, k = 2, sd = 0.5)
  X <- .prepare_matrix(X)
  expect_no_error(.onatski_2009(X, kmax = 5, alpha = 0.05))
})
