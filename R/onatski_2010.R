#' Onatski (2010) Edge Distribution Estimator
#'
#' Estimates the number of factors using the Edge Distribution (ED)
#' estimator of Onatski (2010). The estimator exploits the fact that
#' idiosyncratic eigenvalues of the sample covariance matrix cluster
#' around a single point, while systematic eigenvalues diverge to
#' infinity. The threshold separating the two groups is estimated
#' iteratively using the square root shape of the edge of the
#' eigenvalue distribution.
#'
#' @param eigenvalues Numeric vector of eigenvalues in descending order,
#'   typically obtained from \code{\link{.extract_eigenvalues}}. Must
#'   contain at least \code{kmax + 5} elements to allow the OLS
#'   regression in the calibration step.
#' @param kmax Integer. Maximum number of factors to consider.
#' @param n_iter Integer. Maximum number of iterations for the
#'   calibration procedure. Defaults to \code{4} as recommended
#'   by Onatski (2010).
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{k}{Integer. Estimated number of factors.}
#'     \item{delta}{Numeric. The estimated threshold \eqn{\delta = 2|\hat{\beta}|}.}
#'     \item{beta}{Numeric. The estimated slope coefficient \eqn{\hat{\beta}}
#'       from the OLS regression in the final iteration.}
#'     \item{differences}{Numeric vector of length kmax. Successive
#'       eigenvalue differences \eqn{\lambda_i - \lambda_{i+1}}.}
#'     \item{n_iter}{Integer. Number of iterations performed.}
#'   }
#'
#' @details
#'   The ED estimator of Onatski (2010) is based on the theoretical result
#'   that idiosyncratic eigenvalues cluster around the upper edge
#'   \eqn{u(\mathcal{F}^{c,A,B})} of the limiting spectral distribution,
#'   while systematic eigenvalues diverge. Near the edge, the density of
#'   the limiting spectral distribution behaves like a square root function,
#'   implying that eigenvalue differences \eqn{\lambda_i - \lambda_{i+1}}
#'   for idiosyncratic eigenvalues behave approximately as
#'   \eqn{(an)^{-2/3}}.
#'
#'   The calibration procedure estimates \eqn{\hat{\beta} = (an)^{-2/3}}
#'   by regressing five consecutive eigenvalues \eqn{\lambda_j, \ldots,
#'   \lambda_{j+4}} on a constant and \eqn{(j-1)^{2/3}, \ldots,
#'   (j+3)^{2/3}}, where \eqn{j} is initialized at \eqn{r_{max} + 1}
#'   and updated iteratively.
#'
#'   The estimator requires \code{eigenvalues} to contain at least
#'   \code{kmax + 5} elements so that the OLS window \eqn{j, \ldots, j+4}
#'   is always available.
#'
#' @references
#'   Onatski, A. (2010). Determining the Number of Factors From Empirical
#'   Distribution of Eigenvalues. \emph{The Review of Economics and
#'   Statistics}, 92(4), 1004-1016.
#'
#' @seealso \code{\link{.extract_eigenvalues}}, \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
#' X      <- .prepare_matrix(X)
#' eig    <- .extract_eigenvalues(X, kmax = 13)  # need kmax + 5 eigenvalues
#' result <- .onatski_2010(eig$values, kmax = 8)
#' result$k
#' }
.onatski_2010 <- function(eigenvalues, kmax, n_iter = 4L) {

  # Need at least kmax + 5 eigenvalues for the OLS window
  if (length(eigenvalues) < kmax + 5) {
    stop("eigenvalues must contain at least kmax + 5 = ", kmax + 5,
         " elements. Re-run .extract_eigenvalues() with kmax = ",
         kmax + 4, " or larger.")
  }

  # Successive eigenvalue differences lambda_i - lambda_{i+1}
  diffs <- diff(eigenvalues)               # negative since descending
  diffs <- -diffs[1:kmax]                  # positive differences

  # Iterative calibration procedure
  j     <- kmax + 1L                       # initialize j = rmax + 1
  k_hat <- 0L
  beta  <- NA_real_
  delta <- NA_real_

  for (iter in seq_len(n_iter)) {

    # OLS regression of lambda_j, ..., lambda_{j+4}
    # on constant and (j-1)^(2/3), ..., (j+3)^(2/3)
    idx  <- j:(j + 4L)
    y    <- eigenvalues[idx]
    x    <- ((idx - 1L)^(2/3))
    x_dm <- x - mean(x)
    y_dm <- y - mean(y)

    # OLS slope
    beta  <- sum(x_dm * y_dm) / sum(x_dm^2)
    delta <- 2 * abs(beta)

    # Estimate r(delta)
    exceeds <- which(diffs >= delta)
    k_hat   <- if (length(exceeds) > 0) max(exceeds) else 0L

    # Update j
    j <- k_hat + 1L

    # Ensure j + 4 doesn't exceed available eigenvalues
    if (j + 4L > length(eigenvalues)) {
      j <- length(eigenvalues) - 4L
    }
  }

  list(
    k           = k_hat,
    delta       = delta,
    beta        = beta,
    differences = diffs,
    n_iter      = n_iter
  )
}
