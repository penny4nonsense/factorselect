#' Alessi, Barigozzi and Capasso (2010) Tuned Information Criteria
#'
#' Estimates the number of factors using the tuning-stability procedure of
#' Alessi, Barigozzi and Capasso (2010) applied to the three IC penalty
#' functions of Bai and Ng (2002). For each penalty function, a grid of
#' tuning constants is used and the most stable estimate across the grid
#' is selected as the final estimate.
#'
#' @param eigenvalues Numeric vector of eigenvalues in descending order of
#'   length kmax + 1, typically obtained from \code{\link{.extract_eigenvalues}}.
#' @param V0 Numeric scalar. Total mean squared value of the panel,
#'   \code{sum(X^2) / (N * TT)}, computed from unstandardized demeaned data.
#' @param kmax Integer. Maximum number of factors to consider.
#' @param N Integer. Number of cross-sectional units.
#' @param TT Integer. Number of time periods.
#' @param c_grid Numeric vector. Grid of tuning constants over which to
#'   evaluate stability. Defaults to \code{seq(0, 1, by = 0.01)}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{k_abc1}{Integer. Selected number of factors using ABC with IC1 penalty.}
#'     \item{k_abc2}{Integer. Selected number of factors using ABC with IC2 penalty.}
#'     \item{k_abc3}{Integer. Selected number of factors using ABC with IC3 penalty.}
#'     \item{k_grid_abc1}{Integer vector of length \code{length(c_grid)}.
#'       Selected k for each value of c using IC1 penalty.}
#'     \item{k_grid_abc2}{Integer vector of length \code{length(c_grid)}.
#'       Selected k for each value of c using IC2 penalty.}
#'     \item{k_grid_abc3}{Integer vector of length \code{length(c_grid)}.
#'       Selected k for each value of c using IC3 penalty.}
#'     \item{c_grid}{Numeric vector. The tuning constant grid used.}
#'   }
#'
#' @details
#'   The ABC estimator applies the tuning-stability procedure of Hallin and
#'   Liska (2007) to the IC criteria of Bai and Ng (2002). For each tuning
#'   constant \eqn{c} in the grid, a modified criterion is minimized:
#'   \deqn{IC_j(k, c) = \ln(V(k)) + k \cdot c \cdot g_j(N, T)}
#'   where \eqn{g_j} is the penalty function from \eqn{IC_{pj}} of Bai and
#'   Ng (2002), for j = 1, 2, 3. The final estimate is the modal value of
#'   \eqn{\hat{k}(c)} across the grid — the value of k that is selected
#'   most frequently as c varies.
#'
#'   As with \code{\link{.bai_ng}}, this estimator requires unstandardized
#'   data. The argument \code{V0} should be computed from demeaned but
#'   unstandardized data.
#'
#'   The ABC estimator generally outperforms the raw Bai & Ng IC criteria
#'   in finite samples, particularly when errors are cross-sectionally
#'   correlated.
#'
#' @references
#'   Alessi, L., Barigozzi, M. and Capasso, M. (2010). Improved Penalization
#'   for Determining the Number of Factors in Approximate Factor Models.
#'   \emph{Statistics and Probability Letters}, 80, 1806-1813.
#'
#'   Bai, J. and Ng, S. (2002). Determining the Number of Factors in
#'   Approximate Factor Models. \emph{Econometrica}, 70(1), 191-221.
#'
#'   Hallin, M. and Liska, R. (2007). Determining the Number of Factors in
#'   the Generalized Dynamic Factor Model. \emph{Journal of the American
#'   Statistical Association}, 102, 603-617.
#'
#' @seealso \code{\link{.bai_ng}}, \code{\link{.extract_eigenvalues}},
#'   \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 1)
#' X      <- .prepare_matrix(X, standardize = FALSE)
#' V0     <- sum(X^2) / (100 * 200)
#' eig    <- .extract_eigenvalues(X, kmax = 8)
#' result <- .abc(eig$values, V0 = V0, kmax = 8, N = 100, TT = 200)
#' result$k_abc1
#' result$k_abc2
#' result$k_abc3
#' }
.abc <- function(eigenvalues, V0, kmax, N, TT,
                 c_grid = seq(0, 1, by = 0.01)) {

  m  <- min(N, TT)
  NT <- N * TT

  # V(k) for k = 0, 1, ..., kmax — same as Bai & Ng
  cumulative <- cumsum(eigenvalues[1:kmax])
  V          <- V0 - c(0, cumulative) / N

  # Base penalty functions from Bai & Ng IC criteria
  g1 <- ((N + TT) / NT) * log(NT / (N + TT))
  g2 <- ((N + TT) / NT) * log(m)
  g3 <- log(m) / m

  ks <- 0:kmax

  # For each c in grid, compute k_hat for each penalty function
  k_grid_abc1 <- integer(length(c_grid))
  k_grid_abc2 <- integer(length(c_grid))
  k_grid_abc3 <- integer(length(c_grid))

  for (i in seq_along(c_grid)) {
    c <- c_grid[i]
    k_grid_abc1[i] <- which.min(log(V) + ks * c * g1) - 1L
    k_grid_abc2[i] <- which.min(log(V) + ks * c * g2) - 1L
    k_grid_abc3[i] <- which.min(log(V) + ks * c * g3) - 1L
  }

  # Modal value across grid
  .mode <- function(x) as.integer(names(which.max(table(x))))

  list(
    k_abc1      = .mode(k_grid_abc1),
    k_abc2      = .mode(k_grid_abc2),
    k_abc3      = .mode(k_grid_abc3),
    k_grid_abc1 = k_grid_abc1,
    k_grid_abc2 = k_grid_abc2,
    k_grid_abc3 = k_grid_abc3,
    c_grid      = c_grid
  )
}
