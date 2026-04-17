# tests/testthat/test-simulate_factor_model.R

test_that("simulate_factor_model returns a matrix", {
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  expect_true(is.matrix(X))
})

test_that("simulate_factor_model returns correct dimensions", {
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  expect_equal(dim(X), c(200, 100))
})

test_that("simulate_factor_model returns numeric matrix", {
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  expect_true(is.numeric(X))
})

test_that("simulate_factor_model returns no missing values", {
  X <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  expect_false(anyNA(X))
})

test_that("seed produces reproducible results", {
  X1 <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  X2 <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  expect_equal(X1, X2)
})

test_that("different seeds produce different results", {
  X1 <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 42)
  X2 <- simulate_factor_model(N = 100, TT = 200, k = 3, seed = 99)
  expect_false(identical(X1, X2))
})

test_that("NULL seed runs without error", {
  expect_no_error(simulate_factor_model(N = 100, TT = 200, k = 3, seed = NULL))
})

test_that("N < 2 throws error", {
  expect_error(simulate_factor_model(N = 1, TT = 200, k = 1), "N must")
})

test_that("TT < 2 throws error", {
  expect_error(simulate_factor_model(N = 100, TT = 1, k = 1), "TT must")
})

test_that("k < 1 throws error", {
  expect_error(simulate_factor_model(N = 100, TT = 200, k = 0), "k must")
})

test_that("k >= min(N, TT) throws error", {
  expect_error(simulate_factor_model(N = 10, TT = 200, k = 10), "k must be less")
})

test_that("sd <= 0 throws error", {
  expect_error(simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0), "sd must")
  expect_error(simulate_factor_model(N = 100, TT = 200, k = 3, sd = -1), "sd must")
})

test_that("higher sd produces higher residual variance", {
  X_low  <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.1, seed = 42)
  X_high <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 2.0, seed = 42)
  expect_gt(var(as.vector(X_high)), var(as.vector(X_low)))
})

test_that("works with small N large TT", {
  expect_no_error(simulate_factor_model(N = 10, TT = 500, k = 2, seed = 42))
})

test_that("works with large N small TT", {
  expect_no_error(simulate_factor_model(N = 500, TT = 10, k = 2, seed = 42))
})

test_that("works with N equal to TT", {
  X <- simulate_factor_model(N = 100, TT = 100, k = 3, seed = 42)
  expect_equal(dim(X), c(100, 100))
})

test_that("factor structure is recoverable by select_factors", {
  X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5, seed = 42)
  result <- select_factors(X, method = "ahn_horenstein", kmax = 8)
  expect_equal(result$k[["ahn_horenstein"]], 3L)
})
