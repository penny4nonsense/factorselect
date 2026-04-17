#' Onatski (2009) Test for the Number of Factors
#'
#' Estimates the number of factors using the sequential hypothesis testing
#' procedure of Onatski (2009), applied to the static approximate factor
#' model version described in Section 4 of that paper. The test statistic
#' is based on ratios of differences of adjacent eigenvalues of a
#' complex-valued transformation of the data.
#'
#' @param X Numeric matrix of dimensions T x N, typically preprocessed by
#'   \code{\link{.prepare_matrix}}. Must have an even number of rows.
#' @param kmax Integer. Maximum number of factors to consider. Defines the
#'   upper bound k1 in the sequential testing procedure.
#' @param alpha Numeric. Significance level for the sequential test.
#'   Defaults to \code{0.05}. Must be one of \code{0.01, 0.02, 0.03, 0.04,
#'   0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{k}{Integer. Estimated number of factors from the sequential
#'       testing procedure.}
#'     \item{ratios}{Numeric vector of length kmax. The ratio statistic
#'       \eqn{(\tilde{\gamma}_i - \tilde{\gamma}_{i+1}) /
#'       (\tilde{\gamma}_{i+1} - \tilde{\gamma}_{i+2})} for each i.}
#'     \item{eigenvalues}{Numeric vector of length kmax + 2. Leading
#'       eigenvalues of the complex covariance matrix in descending order.}
#'     \item{critical_value}{Numeric. Critical value used for the test
#'       at the specified significance level.}
#'     \item{alpha}{Numeric. The significance level used.}
#'   }
#'
#' @details
#'   The static approximate factor model version of the Onatski (2009) test
#'   (Section 4) proceeds as follows:
#'
#'   \enumerate{
#'     \item Split the T x N data matrix into two halves of length T/2.
#'     \item Form complex-valued vectors
#'       \eqn{\tilde{X}_j = X_j + i X_{j + T/2}} for \eqn{j = 1, \ldots, T/2}.
#'     \item Compute eigenvalues \eqn{\tilde{\gamma}_i} of
#'       \eqn{\frac{2}{T} \sum_{j=1}^{T/2} \tilde{X}_j \tilde{X}_j^*}.
#'     \item Sequentially test \eqn{H_0: r = k_0} versus
#'       \eqn{H_1: k_0 < r \leq k_{max}} for \eqn{k_0 = 0, 1, \ldots}
#'       using the statistic
#'       \eqn{\tilde{R} = \max_{k_0 < i \leq k_{max}}
#'       (\tilde{\gamma}_i - \tilde{\gamma}_{i+1}) /
#'       (\tilde{\gamma}_{i+1} - \tilde{\gamma}_{i+2})}.
#'     \item Stop when \eqn{H_0} is not rejected. The estimate is the
#'       current \eqn{k_0}.
#'   }
#'
#'   Critical values are taken from Table I of Onatski (2009) and depend
#'   on the significance level \code{alpha} and the number of factors
#'   tested under the alternative \eqn{k_1 - k_0 = k_{max} - k_0}.
#'
#'   If T is odd, the last observation is dropped to ensure equal-length
#'   halves.
#'
#' @references
#'   Onatski, A. (2009). Testing Hypotheses About the Number of Factors
#'   in Large Factor Models. \emph{Econometrica}, 77(5), 1447-1479.
#'
#' @seealso \code{\link{.ahn_horenstein}}, \code{\link{select_factors}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X      <- simulate_factor_model(N = 100, TT = 200, k = 3, sd = 0.5)
#' X      <- .prepare_matrix(X)
#' result <- .onatski_2009(X, kmax = 8, alpha = 0.05)
#' result$k
#' }
.onatski_2009 <- function(X, kmax, alpha = 0.05) {

  # Critical values from Table I of Onatski (2009)
  # Rows: significance levels (15%, 10%, 9%, ..., 1%)
  # Columns: k1 - k0 = 1, 2, 3, 4, 5, 6, 7, 8
  .cv_table <- matrix(c(
    2.75, 3.62, 4.15, 4.54, 4.89, 5.20, 5.45, 5.70,
    3.33, 4.31, 4.91, 5.40, 5.77, 6.13, 6.42, 6.66,
    3.50, 4.49, 5.13, 5.62, 6.03, 6.39, 6.67, 6.92,
    3.69, 4.72, 5.37, 5.91, 6.31, 6.68, 6.95, 7.25,
    3.92, 4.99, 5.66, 6.24, 6.62, 7.00, 7.32, 7.59,
    4.20, 5.31, 6.03, 6.57, 7.00, 7.41, 7.74, 8.04,
    4.52, 5.73, 6.46, 7.01, 7.50, 7.95, 8.29, 8.59,
    5.02, 6.26, 6.97, 7.63, 8.16, 8.61, 9.06, 9.36,
    5.62, 6.91, 7.79, 8.48, 9.06, 9.64, 10.11, 10.44,
    6.55, 8.15, 9.06, 9.93, 10.47, 11.27, 11.75, 12.13,
    8.74, 10.52, 11.67, 12.56, 13.42, 14.26, 14.88, 15.25
  ), nrow = 11, byrow = TRUE)

  .alpha_levels <- c(0.15, 0.10, 0.09, 0.08, 0.07, 0.06,
                     0.05, 0.04, 0.03, 0.02, 0.01)

  # Validate alpha
  if (!alpha %in% .alpha_levels) {
    stop("alpha must be one of: ",
         paste(.alpha_levels, collapse = ", "))
  }

  alpha_row <- which(.alpha_levels == alpha)

  TT <- nrow(X)
  N  <- ncol(X)

  # Drop last observation if T is odd
  if (TT %% 2 != 0) {
    X  <- X[1:(TT - 1), ]
    TT <- TT - 1
  }

  T2 <- TT / 2

  # Form complex matrix X_tilde: T/2 x N
  X1 <- X[1:T2, ]
  X2 <- X[(T2 + 1):TT, ]
  X_tilde <- X1 + 1i * X2   # T/2 x N complex

  # Compute eigenvalues of (2/T) * t(Conj(X_tilde)) %*% X_tilde
  # This is an N x N complex Hermitian matrix
  S <- (2 / TT) * Conj(t(X_tilde)) %*% X_tilde

  # Need kmax + 2 eigenvalues for the ratio statistic
  eig     <- eigen(S, symmetric = FALSE)
  eig_val <- Re(eig$values)   # eigenvalues are real for Hermitian matrix
  eig_val <- sort(eig_val, decreasing = TRUE)
  eigs    <- eig_val[1:(kmax + 2)]

  # Compute ratio statistics for all positions
  # ratio[i] = (gamma[i] - gamma[i+1]) / (gamma[i+1] - gamma[i+2])
  ratios <- (eigs[1:kmax] - eigs[2:(kmax + 1)]) /
    (eigs[2:(kmax + 1)] - eigs[3:(kmax + 2)])

  # Sequential testing procedure (Section 5.3)
  # Test H0: r = k0 vs H1: k0 < r <= kmax for k0 = 0, 1, ...
  k_hat <- 0L

  for (k0 in 0:(kmax - 1)) {

    # Number of factors under alternative: k1 - k0 = kmax - k0
    k1_minus_k0 <- kmax - k0

    # Cap at 8 (maximum column in Table I)
    k1_minus_k0_capped <- min(k1_minus_k0, 8L)

    # Critical value from table
    cv <- .cv_table[alpha_row, k1_minus_k0_capped]

    # Test statistic: max ratio over k0 < i <= kmax
    R_stat <- max(ratios[(k0 + 1):kmax])

    # If we fail to reject H0, stop
    if (R_stat <= cv) {
      k_hat <- k0
      break
    }

    # If we reject all the way to kmax - 1, estimate is kmax
    if (k0 == kmax - 1) {
      k_hat <- kmax
    }
  }

  list(
    k              = k_hat,
    ratios         = ratios,
    eigenvalues    = eigs,
    critical_value = .cv_table[alpha_row, min(kmax, 8L)],
    alpha          = alpha
  )
}
