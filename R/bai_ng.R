#' Bai and Ng (2002) Information Criteria for Number of Factors
#'
#' Estimates the number of factors using the six penalty-based criteria
#' of Bai and Ng (2002). Includes three PC criteria (minimize penalized
#' residual variance) and three IC criteria (minimize penalized log residual
#' variance).
#'
#' @param eigenvalues Numeric vector of eigenvalues in descending order of
#'   length kmax + 1, typically obtained from \code{\link{.extract_eigenvalues}}.
#' @param kmax Integer. Maximum number of factors to consider.
#' @param N Integer. Number of cross-sectional units.
#' @param TT Integer. Number of time periods.
#' @param V0 Numeric scalar. Total mean squared value of the panel,
#'   \code{sum(X^2) / (N * TT)}, computed from unstandardized demeaned data.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{k_pc1}{Integer. Selected number of factors by PC_p1.}
#'     \item{k_pc2}{Integer. Selected number of factors by PC_p2.}
#'     \item{k_pc3}{Integer. Selected number of factors by PC_p3.}
#'     \item{k_ic1}{Integer. Selected number of factors by IC_p1.}
#'     \item{k_ic2}{Integer. Selected number of factors by IC_p2.}
#'     \item{k_ic3}{Integer. Selected number of factors by IC_p3.}
#'     \item{pc1}{Numeric vector of length kmax. Full PC_p1 criterion sequence.}
#'     \item{pc2}{Numeric vector of length kmax. Full PC_p2 criterion sequence.}
#'     \item{pc3}{Numeric vector of length kmax. Full PC_p3 criterion sequence.}
#'     \item{ic1}{Numeric vector of length kmax. Full IC_p1 criterion sequence.}
#'     \item{ic2}{Numeric vector of length kmax. Full IC_p2 criterion sequence.}
#'     \item{ic3}{Numeric vector of length kmax. Full IC_p3 criterion sequence.}
#'   }
#'
#' @details
#'   The six criteria are defined as follows. Let \eqn{V(k)} denote the
#'   residual variance from a k-factor model, \eqn{m = \min(N, T)}, and
#'   \eqn{\hat{\sigma}^2 = V(k_{max})}.
#'
#'   PC criteria (minimize penalized residual variance):
#'   \deqn{PC_{p1}(k) = V(k) + k\hat{\sigma}^2 \frac{N+T}{NT} \ln\left(\frac{NT}{N+T}\right)}
#'   \deqn{PC_{p2}(k) = V(k) + k\hat{\sigma}^2 \frac{N+T}{NT} \ln(m)}
#'   \deqn{PC_{p3}(k) = V(k) + k\hat{\sigma}^2 \frac{\ln(m)}{m}}
#'
#'   IC criteria (minimize penalized log residual variance):
#'   \deqn{IC_{p1}(k) = \ln(V(k)) + k \frac{N+T}{NT} \ln\left(\frac{NT}{N+T}\right)}
#'   \deqn{IC_{p2}(k) = \ln(V(k)) + k \frac{N+T}{NT} \ln(m)}
#'   \deqn{IC_{p3}(k) = \ln(V(k)) + k \frac{\ln(m)}{m}}
#'
#'   \eqn{V(k)} is computed from the eigenvalues of \eqn{XX'/(NT)} as:
#'   \deqn{V(k) = \frac{1}{NT} \sum_{j=k+1}^{m} \lambda_j}
#'   which is the mean residual variance after removing the first k factors.
#'
#'   All six criteria are minimized over \eqn{k = 0, 1, \ldots, k_{max}}.
#'   Note that \eqn{k = 0} is included to allow for the possibility of no
#'   factors.
#'
#'   These estimators require both N and T to be large for consistent
#'   estimation. They may perform poorly when either dimension is small.
#'   For more robust estimation, consider \code{\link{.ahn_horenstein}}.
#'
#' @references
#'   Bai, J. and Ng, S. (2002). Determining the Number of Factors in
#'   Approximate Factor Models. \emph{Econometrica}, 70(1), 191-221.
#'
#' @seealso \code{\link{.ahn_horenstein}}, \code{\link{.extract_eigenvalues}},
#'   \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
#' X      <- .prepare_matrix(X)
#' eig    <- .extract_eigenvalues(X, kmax = 8)
#' result <- .bai_ng(eig$values, kmax = 8, N = 100, TT = 200)
#' result$k_ic1
#' result$k_ic2
#' result$k_ic3
#' }
.bai_ng <- function(eigenvalues, V0, kmax, N, TT) {

  m    <- min(N, TT)
  NT   <- N * TT

  # V(k) = V(0) - (1/N) * sum of top k eigenvalues of XX'/T
  cumulative <- cumsum(eigenvalues[1:kmax])
  V          <- V0 - c(0, cumulative) / N

  # sigma^2 hat = V(kmax)
  sigma2 <- V[kmax + 1]

  # Penalty terms
  g1 <- ((N + TT) / NT) * log(NT / (N + TT))
  g2 <- ((N + TT) / NT) * log(m)
  g3 <- log(m) / m

  ks <- 0:kmax

  # PC criteria
  pc1 <- V + ks * sigma2 * g1
  pc2 <- V + ks * sigma2 * g2
  pc3 <- V + ks * sigma2 * g3

  # IC criteria
  ic1 <- log(V) + ks * g1
  ic2 <- log(V) + ks * g2
  ic3 <- log(V) + ks * g3

  list(
    k_pc1 = which.min(pc1) - 1L,
    k_pc2 = which.min(pc2) - 1L,
    k_pc3 = which.min(pc3) - 1L,
    k_ic1 = which.min(ic1) - 1L,
    k_ic2 = which.min(ic2) - 1L,
    k_ic3 = which.min(ic3) - 1L,
    pc1   = pc1,
    pc2   = pc2,
    pc3   = pc3,
    ic1   = ic1,
    ic2   = ic2,
    ic3   = ic3
  )
}
