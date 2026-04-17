#' Lam and Yao (2012) Eigenvalue Ratio Estimator
#'
#' Estimates the number of factors using the eigenvalue ratio estimator of
#' Lam and Yao (2012). Unlike estimators based on the contemporaneous
#' covariance matrix, this estimator uses lagged auto-covariance matrices,
#' exploiting the fact that the factor loading space is spanned by the
#' eigenvectors of the summed lagged auto-covariance matrix M corresponding
#' to its nonzero eigenvalues.
#'
#' @param X Numeric matrix of dimensions T x N, typically preprocessed by
#'   \code{\link{.prepare_matrix}}.
#' @param kmax Integer. Maximum number of factors to consider.
#' @param h Integer. Number of lags to use in constructing the auto-covariance
#'   matrix M. Defaults to \code{1}. The paper suggests small values are
#'   sufficient; increasing h may improve performance when factors have
#'   strong serial correlation.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{k}{Integer. Selected number of factors.}
#'     \item{ratios}{Numeric vector of length kmax. Full eigenvalue ratio
#'       sequence of M.}
#'     \item{eigenvalues}{Numeric vector of length kmax + 1. Leading
#'       eigenvalues of M in descending order.}
#'   }
#'
#' @details
#'   The estimator constructs the N x N matrix:
#'   \deqn{M = \sum_{k=1}^{h} \hat{\Sigma}_k \hat{\Sigma}_k'}
#'   where \eqn{\hat{\Sigma}_k = T^{-1} \sum_{t=k+1}^{T} x_t x_{t-k}'} is
#'   the lag-k sample auto-covariance matrix.
#'
#'   The factor loading space is spanned by the eigenvectors of M
#'   corresponding to its nonzero eigenvalues, and the number of nonzero
#'   eigenvalues equals the number of factors r (Lam and Yao, 2012,
#'   Proposition 1). In finite samples, the ratio of adjacent eigenvalues
#'   of M spikes at r because eigenvalues r+1 onward are theoretically zero.
#'
#'   The number of factors is estimated as:
#'   \deqn{\hat{r} = \arg\max_{1 \leq k \leq k_{max}} \frac{\lambda_k(M)}{\lambda_{k+1}(M)}}
#'
#' @references
#'   Lam, C. and Yao, Q. (2012). Factor Modelling for High-Dimensional
#'   Time Series: Inference for the Number of Factors.
#'   \emph{The Annals of Statistics}, 40(2), 694-726.
#'
#' @seealso \code{\link{.ahn_horenstein}}, \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
#' X      <- .prepare_matrix(X)
#' result <- .lam_yao(X, kmax = 8, h = 1)
#' result$k
#' }
.lam_yao <- function(X, kmax, h = 1) {

  TT <- nrow(X)
  N  <- ncol(X)

  # Compute M = sum_{k=1}^{h} Sigma_k %*% t(Sigma_k)
  M <- matrix(0, N, N)

  for (lag in 1:h) {
    # Lag-k auto-covariance matrix: N x N
    # Sigma_k = (1/T) * t(X[(lag+1):T, ]) %*% X[1:(T-lag), ]
    X_lead <- X[(lag + 1):TT, ]   # T-lag x N
    X_lag  <- X[1:(TT - lag), ]   # T-lag x N
    Sigma_k <- crossprod(X_lead, X_lag) / TT   # N x N
    M <- M + Sigma_k %*% t(Sigma_k)
  }

  # Eigendecomposition of M — symmetric positive semidefinite
  eig <- eigen(M, symmetric = TRUE)

  # Leading kmax + 1 eigenvalues
  eigs <- eig$values[1:(kmax + 1)]

  # Direct ratio of adjacent eigenvalues
  ratios <- eigs[1:kmax] / eigs[2:(kmax + 1)]

  list(
    k           = which.max(ratios),
    ratios      = ratios,
    eigenvalues = eigs
  )
}
