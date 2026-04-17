# tests/testthat/test-prepare_matrix.R

test_that(".prepare_matrix returns a matrix", {
  X <- matrix(rnorm(200), 20, 10)
  result <- .prepare_matrix(X)
  expect_true(is.matrix(result))
})

test_that(".prepare_matrix returns correct dimensions", {
  X <- matrix(rnorm(200), 20, 10)
  result <- .prepare_matrix(X)
  expect_equal(dim(result), c(20, 10))
})

test_that("individual demeaning removes column means", {
  X <- matrix(rnorm(200), 20, 10)
  result <- .prepare_matrix(X, demean = "individual", standardize = FALSE)
  col_means <- colMeans(result)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("time demeaning removes row means", {
  X <- matrix(rnorm(200), 20, 10)
  result <- .prepare_matrix(X, demean = "time", standardize = FALSE)
  row_means <- rowMeans(result)
  expect_true(all(abs(row_means) < 1e-10))
})

test_that("both demeaning removes column and row means", {
  X <- matrix(rnorm(200), 20, 10)
  result <- .prepare_matrix(X, demean = "both", standardize = FALSE)
  expect_true(all(abs(colMeans(result)) < 1e-10))
  expect_true(all(abs(rowMeans(result)) < 1e-10))
})

test_that("standardize scales columns to unit variance", {
  X <- matrix(rnorm(200), 20, 10)
  result <- .prepare_matrix(X, demean = "none", standardize = TRUE)
  col_sds <- apply(result, 2, sd)
  expect_true(all(abs(col_sds - 1) < 1e-10))
})

test_that("demean = none leaves means unchanged", {
  X <- matrix(rnorm(200), 20, 10)
  col_means_before <- colMeans(X)
  result <- .prepare_matrix(X, demean = "none", standardize = FALSE)
  expect_equal(colMeans(result), col_means_before)
})

test_that("non-matrix input is coerced to matrix", {
  X <- as.data.frame(matrix(rnorm(200), 20, 10))
  result <- .prepare_matrix(X)
  expect_true(is.matrix(result))
})

test_that("non-numeric input throws error", {
  X <- matrix(letters[1:20], 4, 5)
  expect_error(.prepare_matrix(X), "numeric")
})

test_that("missing values throw error", {
  X <- matrix(rnorm(200), 20, 10)
  X[1, 1] <- NA
  expect_error(.prepare_matrix(X), "missing")
})

test_that("single row matrix throws error", {
  X <- matrix(rnorm(10), 1, 10)
  expect_error(.prepare_matrix(X))
})

test_that("single column matrix throws error", {
  X <- matrix(rnorm(10), 10, 1)
  expect_error(.prepare_matrix(X))
})

test_that("zero variance column generates warning", {
  X <- matrix(rnorm(200), 20, 10)
  X[, 1] <- 5  # constant column
  expect_warning(
    .prepare_matrix(X, demean = "none", standardize = TRUE),
    "zero variance"
  )
})
