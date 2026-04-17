#' Ahn-Horenstein Eigenvalue Ratio Estimator
#'
#' Estimates the number of factors using the eigenvalue ratio (ER) and
#' growth ratio (GR) statistics of Ahn and Horenstein (2013). The ratio
#' approach provides robustness to perturbations in the eigenvalue spectrum
#' and performs well when only one dimension (N or T) is large.
#'
#' @param eigenvalues Numeric vector of eigenvalues in descending order of
#'   length kmax + 1, typically obtained from \code{\link{.extract_eigenvalues}}.
#'   Must be positive.
#' @param kmax Integer. Maximum number of factors to consider. The function
#'   evaluates the ratio statistics for k = 1, ..., kmax.
#' @param n Integer. The value of min(N, T), used to compute the mock
#'   eigenvalue boundary term following Ahn and Horenstein (2013) Corollary 1.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{k_er}{Integer. Selected number of factors based on the ER statistic.}
#'     \item{k_gr}{Integer. Selected number of factors based on the GR statistic.}
#'     \item{er}{Numeric vector of length kmax. Full ER statistic sequence.}
#'     \item{gr}{Numeric vector of length kmax. Full GR statistic sequence.}
#'   }
#'
#' @details
#'   The ER statistic is defined as the ratio of successive eigenvalue
#'   differences:
#'   \deqn{ER(k) = \delta_k / \delta_{k+1}}
#'   where \eqn{\delta_k} is the k-th successive difference in the eigenvalue
#'   sequence. The GR statistic replaces raw differences with log growth rates:
#'   \deqn{GR(k) = \log(1 + \delta_k / \lambda_k) / \log(1 + \delta_{k+1} / \lambda_{k+1})}
#'   The boundary case k = 0 is handled by assigning \eqn{\lambda_1 / \log(n)}
#'   as the initial difference term, following Ahn and Horenstein (2013).
#'
#'   The number of factors is selected as the argmax of each statistic over
#'   k = 1, ..., kmax.
#'
#' @references
#'   Ahn, S.C. and Horenstein, A.R. (2013). Eigenvalue Ratio Test for the
#'   Number of Factors. \emph{Econometrica}, 81(3), 1203-1227.
#'
#' @seealso \code{\link{.extract_eigenvalues}}, \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' N <- 100; T <- 200; k_true <- 3
#' Lambda <- matrix(rnorm(N * k_true), N, k_true)
#' F_mat  <- matrix(rnorm(T * k_true), T, k_true)
#' E      <- matrix(rnorm(N * T, sd = 0.5), T, N)
#' X      <- F_mat %*% t(Lambda) + E
#' X      <- .prepare_matrix(X)
#' eig    <- .extract_eigenvalues(X, kmax = 8)
#' result <- .ahn_horenstein(eig$values, kmax = 8, n = min(N, T))
#' result$k_er
#' result$k_gr
#' }
.ahn_horenstein <- function(eigenvalues, kmax, n) {

  eigs <- eigenvalues[1:(kmax + 1)]

  # Mock eigenvalue for k=0: V(0)/log(m) where m = min(N,T)
  V0       <- sum(eigs)
  mu_tilde <- c(V0 / log(n), eigs)   # prepend mock eigenvalue

  # ER statistic: ratio of adjacent eigenvalues
  # ER(k) = mu(k) / mu(k+1), k = 1,...,kmax
  er <- mu_tilde[2:(kmax + 1)] / mu_tilde[3:(kmax + 2)]

  # GR statistic: ratio of log growth rates of V(k)
  # V(k) = sum of eigenvalues from k+1 to m
  # mu*_NTk = mu_NTk / V(k)
  V   <- rev(cumsum(rev(eigs)))         # V(k) for k = 0,...,kmax
  V   <- c(V0, V)                       # prepend V(0)
  mu_star <- mu_tilde[2:(kmax + 2)] / V[1:(kmax + 1)]

  gr <- log(1 + mu_star[1:kmax]) / log(1 + mu_star[2:(kmax + 1)])

  list(
    k_er = which.max(er),
    k_gr = which.max(gr),
    er   = er,
    gr   = gr
  )
}
