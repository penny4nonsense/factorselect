#' Demean and Scale a Matrix for Factor Analysis
#'
#' Removes individual means, time means, or both from a numeric matrix,
#' and optionally scales to unit variance. This is the standard preprocessing
#' step required before eigendecomposition in factor number estimation.
#'
#' @param X Numeric matrix of dimensions T x N (time periods x units).
#' @param demean Character string specifying the demeaning method. One of:
#'   \describe{
#'     \item{"both"}{Remove both individual (column) and time (row) means.
#'       This is the recommended default for macro panels.}
#'     \item{"individual"}{Remove individual (column) means only.}
#'     \item{"time"}{Remove time (row) means only.}
#'     \item{"none"}{No demeaning applied.}
#'   }
#' @param standardize Logical. If \code{TRUE} (default), scale each column
#'   to unit variance after demeaning.
#'
#' @return A demeaned (and optionally scaled) numeric matrix of the same
#'   dimensions as \code{X}.
#'
#' @details
#'   When \code{demean = "both"}, the function iterates individual and time
#'   demeaning to convergence (two passes is sufficient for practical purposes).
#'   This follows the within-transformation used in panel data models.
#'
#' @references
#'   Bai, J. and Ng, S. (2002). Determining the Number of Factors in
#'   Approximate Factor Models. \emph{Econometrica}, 70(1), 191-221.
#'
#' @seealso \code{\link{.extract_eigenvalues}}, \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rnorm(200 * 100, mean = 5), 200, 100)
#' X_clean <- .prepare_matrix(X, demean = "both", standardize = TRUE)
#' }
#'
#' @importFrom stats sd
.prepare_matrix <- function(X, demean = c("both", "individual", "time", "none"),
                            standardize = TRUE) {

  demean <- match.arg(demean)

  # Input validation
  if (!is.matrix(X))        X <- as.matrix(X)
  if (!is.numeric(X))       stop("X must be numeric.")
  if (anyNA(X))             stop("X contains missing values. A balanced panel is required.")
  if (nrow(X) < 2)          stop("X must have at least 2 rows.")
  if (ncol(X) < 2)          stop("X must have at least 2 columns.")

  # Demeaning
  X <- switch(demean,
              individual = {
                X - matrix(colMeans(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
              },
              time = {
                X - rowMeans(X)
              },
              both = {
                # Two-pass within transformation
                X <- X - matrix(colMeans(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
                X <- X - rowMeans(X)
                X
              },
              none = X
  )

  # Standardize columns to unit variance
  if (standardize) {
    col_sd <- apply(X, 2, sd)
    zero_var <- col_sd < .Machine$double.eps
    if (any(zero_var)) {
      warning(sum(zero_var), " column(s) have zero variance after demeaning and will not be scaled.")
      col_sd[zero_var] <- 1
    }
    X <- X / matrix(col_sd, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  }

  X
}
