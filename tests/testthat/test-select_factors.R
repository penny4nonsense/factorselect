# tests/testthat/test-select_factors.R

test_that("select_factors returns a factor_select object", {
  set.seed(42)
  N <- 100; T <- 200; k_true <- 3
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  result <- select_factors(X)
  expect_s3_class(result, "factor_select")
})

test_that("select_factors returns correct list elements", {
  set.seed(42)
  N <- 100; T <- 200; k_true <- 3
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  result <- select_factors(X)
  expect_named(result, c("k", "method", "kmax", "eigenvalues", "details", "call"))
})

test_that("default method is ahn_horenstein", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X)
  expect_equal(result$method, "ahn_horenstein")
})

test_that("default kmax is set automatically", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X)
  expect_true(!is.null(result$kmax))
  expect_true(result$kmax > 0)
})

test_that("user supplied kmax is respected", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X, kmax = 5)
  expect_equal(result$kmax, 5)
})

test_that("eigenvalues are returned", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X, kmax = 5)
  expect_true(is.numeric(result$eigenvalues))
  expect_true(length(result$eigenvalues) > 0)
})

test_that("k is a named integer vector", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X)
  expect_true(is.integer(result$k))
  expect_named(result$k, "ahn_horenstein")
})

test_that("k is within valid range", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X, kmax = 8)
  expect_true(all(result$k >= 1))
  expect_true(all(result$k <= 8))
})

test_that("unknown method throws informative error", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  expect_error(select_factors(X, method = "made_up"), "Unknown method")
})

test_that("unimplemented method throws informative error", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  expect_error(select_factors(X, method = "onatski_2009"), "not yet implemented")
})

test_that("details list contains entry for each method", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X, method = "ahn_horenstein")
  expect_true("ahn_horenstein" %in% names(result$details))
})

test_that("print method runs without error", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X)
  expect_output(print(result))
})

test_that("summary method runs without error", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X)
  expect_output(summary(result))
})

test_that("plot method runs without error", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  result <- select_factors(X)
  expect_silent(plot(result))
})

test_that("select_factors recovers true number of factors", {
  set.seed(42)
  N <- 100; T <- 200; k_true <- 3
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  result <- select_factors(X, method = "ahn_horenstein", kmax = 8)
  expect_equal(result$k[["ahn_horenstein"]], k_true)
})

test_that("select_factors works with small N large T", {
  set.seed(42)
  N <- 25; T <- 200; k_true <- 2
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  result <- select_factors(X, method = "ahn_horenstein", kmax = 8)
  expect_equal(result$k[["ahn_horenstein"]], k_true)
})

test_that("select_factors works with large N small T", {
  set.seed(42)
  N <- 200; T <- 25; k_true <- 2
  Lambda <- matrix(rnorm(N * k_true), N, k_true)
  F_mat  <- matrix(rnorm(T * k_true), T, k_true)
  E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
  X      <- F_mat %*% t(Lambda) + E
  result <- select_factors(X, method = "ahn_horenstein", kmax = 8)
  expect_equal(result$k[["ahn_horenstein"]], k_true)
})

test_that("demean argument is passed through correctly", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  expect_no_error(select_factors(X, demean = "individual"))
  expect_no_error(select_factors(X, demean = "time"))
  expect_no_error(select_factors(X, demean = "none"))
  expect_no_error(select_factors(X, demean = "both"))
})

test_that("standardize = FALSE runs without error", {
  set.seed(42)
  X <- matrix(rnorm(2000), 200, 100)
  expect_no_error(select_factors(X, standardize = FALSE))
})
