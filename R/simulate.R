#' Simulate Data from an Approximate Factor Model
#'
#' Generates a simulated panel data matrix from a static approximate factor
#' model. Useful for testing and benchmarking factor number estimators.
#'
#' @param N Integer. Number of cross-sectional units.
#' @param TT Integer. Number of time periods. Named \code{TT} to avoid
#'   conflict with the base R function \code{T} (which evaluates to
#'   \code{TRUE}).
#' @param k Integer. True number of factors.
#' @param sd Numeric. Standard deviation of the idiosyncratic error term.
#'   Defaults to \code{1}. Lower values produce stronger signal relative
#'   to noise.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'   Defaults to \code{NULL} (no seed set).
#'
#' @return A numeric matrix of dimensions \code{TT x N} generated from:
#'   \deqn{X = F \Lambda' + E}
#'   where \eqn{F} is a \code{TT x k} matrix of factors drawn from
#'   \eqn{N(0,1)}, \eqn{\Lambda} is an \code{N x k} matrix of loadings
#'   drawn from \eqn{N(0,1)}, and \eqn{E} is a \code{TT x N} matrix of
#'   idiosyncratic errors drawn from \eqn{N(0, sd^2)}.
#'
#' @details
#'   The data generating process follows the standard approximate factor
#'   model of Chamberlain and Rothschild (1983) as used in the simulation
#'   exercises of Ahn and Horenstein (2013). Factors and loadings are
#'   independent standard normal draws. Errors are i.i.d. normal with
#'   mean zero and standard deviation \code{sd}.
#'
#'   The signal-to-noise ratio is controlled by \code{sd} — smaller values
#'   produce a cleaner factor structure that is easier for estimators to
#'   recover. The default \code{sd = 1} matches the baseline simulation
#'   design of Ahn and Horenstein (2013) with \code{theta = 1}.
#'
#' @references
#'   Ahn, S.C. and Horenstein, A.R. (2013). Eigenvalue Ratio Test for the
#'   Number of Factors. \emph{Econometrica}, 81(3), 1203-1227.
#'
#'   Chamberlain, G. and Rothschild, M. (1983). Arbitrage, Factor Structure,
#'   and Mean-Variance Analysis on Large Asset Markets.
#'   \emph{Econometrica}, 51(5), 1281-1304.
#'
#' @seealso \code{\link{select_factors}}
#'
#' @export
#'
#' @examples
#' # Simulate a factor model with 3 factors
#' X <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5, seed = 42)
#' dim(X)
#'
#' # Pass directly to select_factors
#' result <- select_factors(X)
#' result$k
simulate_factor_model <- function(N, TT, k, sd = 1, seed = NULL) {

  # Input validation
  if (!is.numeric(N) || length(N) != 1 || N < 2)
    stop("N must be a single integer >= 2.")
  if (!is.numeric(TT) || length(TT) != 1 || TT < 2)
    stop("TT must be a single integer >= 2.")
  if (!is.numeric(k) || length(k) != 1 || k < 1)
    stop("k must be a single integer >= 1.")
  if (k >= min(N, TT))
    stop("k must be less than min(N, TT).")
  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0)
    stop("sd must be a single positive number.")

  if (!is.null(seed)) set.seed(seed)

  Lambda <- matrix(rnorm(N  * k),  N,  k)
  F_mat  <- matrix(rnorm(TT * k),  TT, k)
  E      <- matrix(rnorm(TT * N, sd = sd), TT, N)

  F_mat %*% t(Lambda) + E
}
