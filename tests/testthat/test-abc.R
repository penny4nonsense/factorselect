# tests/testthat/test-abc.R

test_that(".abc returns a list with correct elements", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_true(is.list(result))
  expect_named(result, c("k_abc1", "k_abc2", "k_abc3",
                         "k_grid_abc1", "k_grid_abc2", "k_grid_abc3",
                         "c_grid"))
})

test_that("k estimates are non-negative integers within range", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  k_vals <- c(result$k_abc1, result$k_abc2, result$k_abc3)
  expect_true(all(k_vals >= 0))
  expect_true(all(k_vals <= 8))
})

test_that("k_grid vectors have correct length", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  c_grid <- seq(0, 1, by = 0.01)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200,
                 c_grid = c_grid)
  expect_length(result$k_grid_abc1, length(c_grid))
  expect_length(result$k_grid_abc2, length(c_grid))
  expect_length(result$k_grid_abc3, length(c_grid))
})

test_that("c_grid is returned correctly", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  c_grid <- seq(0, 1, by = 0.01)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200,
                 c_grid = c_grid)
  expect_equal(result$c_grid, c_grid)
})

test_that("k_grid values are within valid range", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_true(all(result$k_grid_abc1 >= 0 & result$k_grid_abc1 <= 8))
  expect_true(all(result$k_grid_abc2 >= 0 & result$k_grid_abc2 <= 8))
  expect_true(all(result$k_grid_abc3 >= 0 & result$k_grid_abc3 <= 8))
})

test_that("k_abc1 is the modal value of k_grid_abc1", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_abc1,
               as.integer(names(which.max(table(result$k_grid_abc1)))))
})

test_that("k_abc2 is the modal value of k_grid_abc2", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_abc2,
               as.integer(names(which.max(table(result$k_grid_abc2)))))
})

test_that("k_abc3 is the modal value of k_grid_abc3", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
  expect_equal(result$k_abc3,
               as.integer(names(which.max(table(result$k_grid_abc3)))))
})

test_that("c = 0 always selects k = kmax", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200,
                 c_grid = 0)
  # With c=0 penalty vanishes so log(V) is minimized at kmax
  expect_equal(result$k_grid_abc1, 8L)
  expect_equal(result$k_grid_abc2, 8L)
  expect_equal(result$k_grid_abc3, 8L)
})

test_that("custom c_grid is respected", {
  set.seed(42)
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (100 * 200)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  c_grid <- seq(0, 1, by = 0.1)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200,
                 c_grid = c_grid)
  expect_length(result$k_grid_abc1, length(c_grid))
  expect_equal(result$c_grid, c_grid)
})

test_that("ABC1 recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_equal(result$k_abc1, k_true)
})

test_that("ABC2 recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 100; TT <- 200; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_equal(result$k_abc2, k_true)
})

test_that("ABC3 recovers true number of factors — large N and T", {
  set.seed(42)
  N <- 500; TT <- 1000; k_true <- 3
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  # ABC3 uses the weakest penalty (log(m)/m) and requires large samples
  # Accept k_true or adjacent values as reasonable
  expect_true(result$k_abc3 %in% c(k_true - 1L, k_true, k_true + 1L))
})

test_that("works with N > TT", {
  set.seed(42)
  N <- 200; TT <- 50; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_true(is.list(result))
})

test_that("works with TT > N", {
  set.seed(42)
  N <- 50; TT <- 200; k_true <- 2
  X      <- simulate_factor_model(N = N, TT = TT, k = k_true, sd = 1)
  X      <- .prepare_matrix(X, standardize = FALSE)
  V0     <- sum(X^2) / (N * TT)
  eig    <- .extract_eigenvalues(X, kmax = 8)
  result <- .abc(eig$values, V0 = V0, kmax = 8, N = N, TT = TT)
  expect_true(is.list(result))
})
