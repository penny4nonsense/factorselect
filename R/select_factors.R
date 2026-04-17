#' Select the Number of Factors in an Approximate Factor Model
#'
#' A unified interface for estimating the number of factors in a large
#' dimensional approximate factor model. Preprocesses the data and dispatches
#' to one or more factor number estimators.
#'
#' @param X A numeric matrix of dimensions T x N (time periods x units), or
#'   an object coercible to a numeric matrix. Must be a balanced panel with no
#'   missing values.
#' @param method Character vector specifying which estimator(s) to use. One or
#'   more of \code{"ahn_horenstein"}, \code{"bai_ng"}, \code{"onatski_2009"},
#'   \code{"onatski_2010"}, \code{"abc"}, \code{"lam_yao"}. Defaults to
#'   \code{"ahn_horenstein"}.
#' @param kmax Integer. Maximum number of factors to consider. Defaults to
#'   \code{NULL}, in which case it is set to \code{min(floor(sqrt(min(N, T))),
#'   8)}.
#' @param demean Character string passed to \code{.prepare_matrix()}. One of
#'   \code{"both"}, \code{"individual"}, \code{"time"}, \code{"none"}.
#'   Defaults to \code{"both"} as recommended by Ahn and Horenstein (2013).
#' @param standardize Logical. Whether to standardize columns to unit variance
#'   before estimation. Defaults to \code{TRUE}. Note that \code{bai_ng},
#'   \code{abc}, and \code{lam_yao} always use unstandardized data regardless
#'   of this setting.
#' @param h Integer. Number of lags to use for the \code{lam_yao} estimator.
#'   Defaults to \code{1}. Ignored for all other methods.
#' @param alpha Numeric. Significance level for the \code{onatski_2009}
#'   sequential test. Defaults to \code{0.05}. Must be one of
#'   \code{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15}.
#'   Ignored for all other methods.
#'
#' @return An object of class \code{"factor_select"}, which is a named list
#'   with the following elements:
#'   \describe{
#'     \item{k}{Named integer vector of selected factor numbers, one per
#'       method requested.}
#'     \item{method}{Character vector of methods used.}
#'     \item{kmax}{Integer. Maximum number of factors considered.}
#'     \item{eigenvalues}{Numeric vector of leading eigenvalues.}
#'     \item{details}{Named list of full output from each estimator.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @details
#'   The data are first preprocessed via \code{.prepare_matrix()} and then
#'   a single eigendecomposition is performed via \code{.extract_eigenvalues()},
#'   which is shared across all requested estimators for efficiency.
#'
#'   The default method is \code{"ahn_horenstein"}, which is recommended for
#'   most applications. It is robust to perturbations in the eigenvalue
#'   spectrum and performs well when only one of N or T is large.
#'
#'   The \code{"bai_ng"}, \code{"abc"}, and \code{"lam_yao"} methods always
#'   use unstandardized data because their penalty terms and auto-covariance
#'   structure depend on the actual scale of the data.
#'
#' @references
#'   Ahn, S.C. and Horenstein, A.R. (2013). Eigenvalue Ratio Test for the
#'   Number of Factors. \emph{Econometrica}, 81(3), 1203-1227.
#'
#'   Bai, J. and Ng, S. (2002). Determining the Number of Factors in
#'   Approximate Factor Models. \emph{Econometrica}, 70(1), 191-221.
#'
#'   Alessi, L., Barigozzi, M. and Capasso, M. (2010). Improved Penalization
#'   for Determining the Number of Factors in Approximate Factor Models.
#'   \emph{Statistics and Probability Letters}, 80, 1806-1813.
#'
#'   Lam, C. and Yao, Q. (2012). Factor Modelling for High-Dimensional
#'   Time Series: Inference for the Number of Factors.
#'   \emph{The Annals of Statistics}, 40(2), 694-726.
#'
#' @seealso \code{\link{.ahn_horenstein}}, \code{\link{.bai_ng}},
#'   \code{\link{.abc}}, \code{\link{.lam_yao}},
#'   \code{\link{.prepare_matrix}}, \code{\link{.extract_eigenvalues}}
#'
#' @importFrom graphics abline legend
#' @export
#'
#' @examples
#' set.seed(42)
#' N <- 100; T <- 200; k_true <- 3
#' Lambda <- matrix(rnorm(N * k_true), N, k_true)
#' F_mat  <- matrix(rnorm(T * k_true), T, k_true)
#' E      <- matrix(rnorm(T * N, sd = 0.5), T, N)
#' X      <- F_mat %*% t(Lambda) + E
#' select_factors(X)
select_factors <- function(X,
                           method      = "ahn_horenstein",
                           kmax        = NULL,
                           demean      = c("both", "individual", "time", "none"),
                           standardize = TRUE,
                           h           = 1L,
                           alpha       = 0.05) {

  call   <- match.call()
  demean <- match.arg(demean)

  # Validate method
  valid_methods <- c("ahn_horenstein", "bai_ng", "onatski_2009",
                     "onatski_2010", "abc", "lam_yao")
  bad <- setdiff(method, valid_methods)
  if (length(bad) > 0) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         "\nValid methods: ", paste(valid_methods, collapse = ", "))
  }

  # Validate h
  if (!is.numeric(h) || length(h) != 1 || h < 1) {
    stop("h must be a single positive integer.")
  }

  # Validate alpha
  valid_alphas <- c(0.15, 0.10, 0.09, 0.08, 0.07, 0.06,
                    0.05, 0.04, 0.03, 0.02, 0.01)
  if (!alpha %in% valid_alphas) {
    stop("alpha must be one of: ", paste(valid_alphas, collapse = ", "))
  }

  # Preprocess — standardized version for most estimators
  X_clean <- .prepare_matrix(X, demean = demean, standardize = standardize)

  T_dim <- nrow(X_clean)
  N_dim <- ncol(X_clean)
  n     <- min(N_dim, T_dim)

  # Set default kmax
  if (is.null(kmax)) {
    kmax <- min(floor(sqrt(n)), 8)
  }

  # Single eigendecomposition shared across standardized estimators
  eig <- .extract_eigenvalues(X_clean, kmax = kmax)

  # Bai & Ng, ABC and Lam-Yao require unstandardized data
  if (any(c("bai_ng", "abc", "lam_yao") %in% method)) {
    X_bn   <- .prepare_matrix(X, demean = demean, standardize = FALSE)
    V0     <- sum(X_bn^2) / (N_dim * T_dim)
    eig_bn <- .extract_eigenvalues(X_bn, kmax = kmax)
  }

  # Dispatch to each requested estimator
  details  <- list()
  k        <- integer(length(method))
  names(k) <- method

  for (m in method) {
    result <- switch(m,
                     ahn_horenstein = {
                       res  <- .ahn_horenstein(eig$values, kmax = kmax, n = n)
                       k[m] <- res$k_gr
                       res
                     },
                     bai_ng = {
                       res  <- .bai_ng(eig_bn$values, V0 = V0, kmax = kmax,
                                       N = N_dim, TT = T_dim)
                       k[m] <- res$k_ic1
                       res
                     },
                     abc = {
                       res  <- .abc(eig_bn$values, V0 = V0, kmax = kmax,
                                    N = N_dim, TT = T_dim)
                       k[m] <- res$k_abc1
                       res
                     },
                     lam_yao = {
                       res  <- .lam_yao(X_bn, kmax = kmax, h = h)
                       k[m] <- res$k
                       res
                     },
                     onatski_2009 = {
                       res  <- .onatski_2009(X_clean, kmax = kmax, alpha = alpha)
                       k[m] <- res$k
                       res
                     },
                     onatski_2010 = {
                       eig_ed <- .extract_eigenvalues(X_clean, kmax = kmax + 4)
                       res    <- .onatski_2010(eig_ed$values, kmax = kmax)
                       k[m]   <- res$k
                       res
                     }
    )
    details[[m]] <- result
  }

  structure(
    list(
      k           = k,
      method      = method,
      kmax        = kmax,
      eigenvalues = eig$values,
      details     = details,
      call        = call
    ),
    class = "factor_select"
  )
}


#' Print Method for factor_select Objects
#'
#' @param x A \code{factor_select} object.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.factor_select <- function(x, ...) {
  cat("Factor Number Selection\n")
  cat("=======================\n")
  cat("Call: "); print(x$call)
  cat("\nkmax:", x$kmax, "\n")
  cat("\nEstimated number of factors:\n")
  for (m in x$method) {
    cat(sprintf("  %-20s %d\n", m, x$k[m]))
  }
  invisible(x)
}


#' Summary Method for factor_select Objects
#'
#' @param object A \code{factor_select} object.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.factor_select <- function(object, ...) {
  print(object)
  cat("\nLeading eigenvalues:\n")
  cat(round(object$eigenvalues, 4), "\n")
  invisible(object)
}


#' Plot Method for factor_select Objects
#'
#' Produces a scree plot of the leading eigenvalues with the selected
#' number of factors marked.
#'
#' @param x A \code{factor_select} object.
#' @param main Character string. Plot title. Defaults to \code{"Scree Plot"}.
#' @param ... Further arguments passed to \code{plot()}.
#' @importFrom graphics abline legend
#' @export
plot.factor_select <- function(x, main = "Scree Plot", ...) {
  eigs <- x$eigenvalues
  k    <- x$k

  plot(seq_along(eigs), eigs,
       type = "b",
       xlab = "Factor Number",
       ylab = "Eigenvalue",
       main = main,
       pch  = 19, ...)

  cols <- c("red", "blue", "darkgreen", "purple", "orange", "brown")
  for (i in seq_along(x$method)) {
    abline(v   = k[x$method[i]],
           col = cols[i],
           lty = 2,
           lwd = 1.5)
  }

  legend("topright",
         legend = paste(x$method, "=", k),
         col    = cols[seq_along(x$method)],
         lty    = 2,
         lwd    = 1.5,
         bty    = "n")

  invisible(x)
}
