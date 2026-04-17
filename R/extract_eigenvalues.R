#' Extract Leading Eigenvalues from a Panel Data Matrix
#'
#' Computes the leading eigenvalues of the sample covariance matrix using
#' a truncated eigendecomposition. Automatically selects the smaller of the
#' N x N or T x T covariance matrix for efficiency. Uses RSpectra when
#' available for large matrices, falling back to base R otherwise.
#'
#' @param X Numeric matrix of dimensions T x N, typically preprocessed by
#'   \code{\link{.prepare_matrix}}.
#' @param kmax Integer. Number of leading eigenvalues to compute. Should be
#'   set generously (e.g., 8-15) to allow estimators to evaluate the full
#'   candidate range.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{values}{Numeric vector of length \code{kmax + 1} containing the
#'       leading eigenvalues in descending order. The extra eigenvalue is
#'       required by ratio-based estimators.}
#'     \item{vectors}{Numeric matrix of corresponding eigenvectors.}
#'     \item{orientation}{Character string, either \code{"N"} or \code{"T"},
#'       indicating which covariance matrix was decomposed.}
#'   }
#'
#' @details
#'   When N <= T, decomposes the N x N matrix \eqn{XX'/T}.
#'   When N > T, decomposes the T x T matrix \eqn{X'X/N}.
#'   This ensures the cheaper decomposition is always used.
#'
#'   RSpectra's \code{eigs_sym()} is used when available and when
#'   \code{min(N, T) > 100}, as the truncated decomposition only provides
#'   meaningful speedup at larger scales.
#'
#' @references
#'   Ahn, S.C. and Horenstein, A.R. (2013). Eigenvalue Ratio Test for the
#'   Number of Factors. \emph{Econometrica}, 81(3), 1203-1227.
#'
#' @seealso \code{\link{.prepare_matrix}}, \code{\link{.ahn_horenstein}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' N <- 100; T <- 200
#' X <- matrix(rnorm(T * N), T, N)
#' X <- .prepare_matrix(X)
#' eig <- .extract_eigenvalues(X, kmax = 8)
#' eig$values
#' }
.extract_eigenvalues <- function(X, kmax) {

  T_dim <- nrow(X)
  N_dim <- ncol(X)

  # Validate kmax
  max_possible <- min(N_dim, T_dim) - 1
  if (kmax >= min(N_dim, T_dim)) {
    warning(
      "kmax (", kmax, ") is too large. ",
      "Reducing to min(N, T) - 1 = ", max_possible, "."
    )
    kmax <- max_possible
  }

  # Always decompose the smaller covariance matrix
  if (N_dim <= T_dim) {
    S           <- crossprod(X) / T_dim    # N x N  — X'X
    orientation <- "N"
  } else {
    S           <- tcrossprod(X) / N_dim   # T x T  — XX'
    orientation <- "T"
  }

  # Need kmax + 1 eigenvalues for ratio-based estimators
  k_compute <- kmax + 1

  # Use RSpectra if available and matrix is large enough to benefit
  use_rspectra <- requireNamespace("RSpectra", quietly = TRUE) &&
    nrow(S) > 100

  if (use_rspectra) {
    eig <- RSpectra::eigs_sym(S, k = k_compute)
  } else {
    eig_full    <- eigen(S, symmetric = TRUE)
    eig         <- list(
      values  = eig_full$values[1:k_compute],
      vectors = eig_full$vectors[, 1:k_compute, drop = FALSE]
    )
  }

  list(
    values      = eig$values,
    vectors     = eig$vectors,
    orientation = orientation
  )
}
