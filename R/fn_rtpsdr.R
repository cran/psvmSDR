#'Real time sufficient dimension reduction through principal least squares SVM
#'@description
#' This function implements a real-time version of principal SDR based on
#' least squares SVM loss. It is intended for streaming or sequential data
#' settings where new observations arrive continuously and re-fitting the full
#' SDR model would be computationally expensive.
#'
#' After an initial psdr or rtpsdr fit is obtained, this function updates the
#' working matrix M, slice statistics, and eigen-decomposition efficiently
#' using only the new batch of data. The method supports both regression and
#' binary classification, automatically choosing the appropriate LS-SVM variant.
#'
#' The returned object includes cumulative sample size, updated mean vector,
#' slice coefficients, intermediate matrices required for updates, and the
#' resulting central subspace basis.
#'@param x x in new data
#'@param y y in new data, y is continuous
#'@param obj the latest output object from the \code{rtpsdr}
#'@param h unified control for slicing or weighting; accepts either an integer or a numeric vector.
#'@param lambda hyperparameter for the loss function. default is set to 1.
#' @author Jungmin Shin, \email{c16267@gmail.com}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#' @references Artemiou, A. and Dong, Y. (2016)
#' \emph{Sufficient dimension reduction via principal lq support vector machine,
#'  Electronic Journal of Statistics 10: 783–805}.\cr
#'  Artemiou, A., Dong, Y. and Shin, S. J. (2021)
#' \emph{Real-time sufficient dimension reduction through principal least
#'  squares support vector machines, Pattern Recognition 112: 107768}.\cr
#'  Kim, B. and Shin, S. J. (2019)
#' \emph{Principal weighted logistic regression for sufficient dimension
#' reduction in binary classification, Journal of the Korean Statistical Society 48(2): 194–206}.\cr
#'  Li, B., Artemiou, A. and Li, L. (2011)
#' \emph{Principal support vector machines for linear and
#' nonlinear sufficient dimension reduction, Annals of Statistics 39(6): 3182–3210}.\cr
#' Soale, A.-N. and Dong, Y. (2022)
#' \emph{On sufficient dimension reduction via principal asymmetric
#'  least squares, Journal of Nonparametric Statistics 34(1): 77–94}.\cr
#'  Wang, C., Shin, S. J. and Wu, Y. (2018)
#' \emph{Principal quantile regression for sufficient dimension
#'  reduction with heteroscedasticity, Electronic Journal of Statistics 12(2): 2114–2140}.\cr
#'  Shin, S. J., Wu, Y., Zhang, H. H. and Liu, Y. (2017)
#' \emph{Principal weighted support vector machines for sufficient dimension reduction in
#'  binary classification, Biometrika 104(1): 67–81}. \cr
#'  Li, L. (2007)
#' \emph{Sparse sufficient dimension reduction, Biometrika 94(3): 603–613}.
#' @return An object of class \code{c("rtpsdr","psdr")} containing:
#' \itemize{
#'   \item \code{x}, \code{y}: latest batch data
#'   \item \code{M}: working matrix
#'   \item \code{evalues}, \code{evectors}: eigen-decomposition of \code{M} (central subspace basis)
#'   \item \code{N}: cumulative sample size
#'   \item \code{Xbar}: cumulative mean vector
#'   \item \code{r}: slice-specific coefficient matrix
#'   \item \code{A}: new A part for update. See Artemiou et. al., (2021)
#'   \item \code{loss}: "lssvm" (continuous) or "wlssvm" (binary)
#'   \item \code{fit}: metadata (mode="realtime", H, cutpoints, weight_cutpoints, lambda, etc.)
#' }
#'
#'@seealso \code{\link{psdr}}, \code{\link{npsdr}}
#'@examples
#'\donttest{
#' set.seed(1)
#' p <- 5; m <- 300; B <- 3
#' obj <- NULL
#' for (b in 1:B) {
#'   x <- matrix(rnorm(m*p), m, p)
#'   y <- x[,1]/(0.5+(x[,2]+1)^2) + 0.2*rnorm(m)
#'   obj <- rtpsdr(x, y, obj=obj, h=8, lambda=1)
#' }
#' print(obj)
#' summary(obj)
#'}
#'@import stats
#' @importFrom utils head
#'@export rtpsdr

# -----------------------------------------------------------------------------------------------------------
# Real time Principal (weighted) least squares SDR
# -----------------------------------------------------------------------------------------------------------

rtpsdr <- function(x, y, obj = NULL, h = 10, lambda = 1) {

  .check_input(x, y)
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0)
    stop("lambda must be a positive scalar.")

  ytype <- .response_type(y)
  if (!ytype %in% c("binary", "continuous"))
    stop("Unsupported response type for rtpsdr(): ", ytype)

  n_new <- nrow(x); p <- ncol(x)

  if (ytype == "binary") {
    if (is.list(h)) {
      weight_cutpoints <- if (!is.null(h$weight)) h$weight else NULL
    } else if (is.numeric(h) && length(h) > 1) {
      if (any(h < 0 | h > 1))
        stop("When y is binary, h must contain numeric values in [0, 1].")
      weight_cutpoints <- sort(unique(h))
    } else if (length(h) == 1L && is.numeric(h)) {
      weight_cutpoints <- seq(0, 1, length = as.integer(h) + 2)[2:(as.integer(h) + 1)]
    } else {
      stop("Invalid input for h.")
    }

    weight_list <- weight_cutpoints
    pi.grid <- weight_list
    H <- length(weight_list)
    qy <- NULL

  } else {
    if (is.list(h)) {
      cutpoints <- if (!is.null(h$slice)) h$slice else NULL
    } else if (is.numeric(h) && length(h) > 1) {
      if (any(h <= 0 | h >= 1)) {
        h[h <= 0] <- 0.01
        h[h >= 1] <- 0.99
      }
      cutpoints <- sort(unique(h))
    } else if (length(h) == 1L && is.numeric(h)) {
      step <- 1 / as.integer(h)
      cutpoints <- seq(step, 1 - step, by = step)
    } else {
      stop("Invalid input for h.")
    }

    pi.grid <- cutpoints
    H <- length(pi.grid) + 1L
    weight_list <- seq(0, 1, length = H + 2)[2:(H + 1)]
    qprob <- (1:(H - 1)) / H
    qy <- stats::quantile(y, qprob)
  }

  if (is.null(obj)) {
    init_loss <- if (ytype == "continuous") "lssvm" else "wlssvm"

    init_fit <- psdr(
      x, y,
      loss = init_loss,
      h = h, lambda = lambda,
      plot = FALSE
    )

    init_fit$fit$mode <- "realtime"
    init_fit$fit$H <- H
    init_fit$fit$cutpoints <- if (ytype == "continuous") pi.grid else NULL
    init_fit$fit$weight_cutpoints <- if (ytype == "binary") weight_list else NULL
    init_fit$fit$lambda <- lambda
    init_fit$N <- n_new
    init_fit$Xbar <- colMeans(x)

    class(init_fit) <- c("rtpsdr", "psdr", class(init_fit))
    return(invisible(init_fit))
  }

  if (!inherits(obj, "rtpsdr") && !inherits(obj, "psdr"))
    stop("obj must be a previous rtpsdr/psdr object.")

  loss <- obj$loss
  if (ytype == "continuous" && loss != "lssvm") loss <- "lssvm"
  if (ytype == "binary" && loss != "wlssvm") loss <- "wlssvm"

  N_old <- if (!is.null(obj$N)) obj$N else nrow(obj$x)
  Xbar_old <- if (!is.null(obj$Xbar)) obj$Xbar else colMeans(obj$x)
  Xbar_new <- (Xbar_old * N_old + colMeans(x) * n_new) / (N_old + n_new)

  prep <- .rt_prepare(x, Xbar = Xbar_new)
  x.star <- prep$x.star
  cov.new.Xstar <- prep$cov.x.star * (n_new - 1) / n_new  # bias correction

  A_old <- obj$A
  r_old <- obj$r
  if (is.null(A_old) || is.null(r_old))
    stop("obj must contain fields A and r from a previous rtpsdr/psdr fit.")

  if (ytype == "continuous") {
    qprob <- (1:(H - 1)) / H
    qy <- stats::quantile(y, qprob)
  }

  A_new <- vector("list", length = H)
  r_new <- matrix(0, nrow = H, ncol = p + 1)
  Ipp1 <- diag(p + 1)

  # --- Case 1: Continuous response (lssvm) ---
  if (loss == "lssvm") {
    for (s in 1:H) {
      thr_idx <- min(s, length(qy))
      y.tilde <- ifelse(y < qy[thr_idx], -1, 1)

      B_new <- n_new * cov.new.Xstar / lambda + t(x.star) %*% x.star
      C_new <- t(x.star) %*% y.tilde
      AB <- A_old[[s]] %*% B_new
      s_part <- Ipp1 + AB
      K <- Ipp1 - A_old[[s]] %*% B_new %*% .safe_solve(s_part, Ipp1)

      r_new[s, ] <- t(K %*% (r_old[s, ] + A_old[[s]] %*% C_new))
      A_new[[s]] <- K %*% A_old[[s]]
    }

    # --- Case 2: Binary response (wlssvm) ---
  } else if (loss == "wlssvm") {
    .to_pm_binary <- function(y) {
      y0 <- stats::na.omit(y)
      if (is.logical(y0)) return(ifelse(y, 1, -1))
      if (is.factor(y0) || is.character(y0)) {
        lev <- if (is.factor(y0)) levels(y0) else sort(unique(as.character(y0)))
        return(ifelse(as.character(y) == lev[2], 1, -1))
      }
      u <- sort(unique(y0))
      if (length(u) != 2L) stop("Binary coercion requires exactly 2 unique values.")
      ifelse(y == u[2], 1, -1)
    }
    y.bin <- .to_pm_binary(y)

    for (i in 1:H) {
      wgt <- weight_list[i]
      W <- diag(ifelse(y.bin == 1, 1 - wgt, wgt), nrow = n_new, ncol = n_new)
      B_new <- n_new * cov.new.Xstar / lambda + t(x.star) %*% W %*% x.star
      C_new <- t(x.star) %*% W %*% y.bin
      AB <- A_old[[i]] %*% B_new
      s_part <- Ipp1 + AB
      S <- Ipp1 - A_old[[i]] %*% B_new %*% .safe_solve(s_part, Ipp1)

      r_new[i, ] <- t(S %*% (r_old[i, ] + A_old[[i]] %*% C_new))
      A_new[[i]] <- S %*% A_old[[i]]
    }

  } else {
    stop("rtpsdr currently supports only lssvm (continuous) and wlssvm (binary).")
  }

  Mn <- t(r_new[, 1:p, drop = FALSE]) %*% r_new[, 1:p, drop = FALSE]
  eMn <- eigen(Mn)

  # Metadata
    meta <- list(
      mode = "realtime",
      H = H,
      cutpoints = if (ytype == "continuous") pi.grid else NULL,
      weight_cutpoints = if (ytype == "binary") weight_list else NULL,
      lambda = lambda,
      ytype = ytype
    )

  # Output object
  out <- list(
    loss = loss,
    x = x, y = y,
    M = Mn,
    evalues = eMn$values,
    evectors = eMn$vectors,
    r = r_new,
    A = A_new,
    N = N_old + n_new,
    Xbar = Xbar_new,
    fit = meta
  )
  class(out) <- c("rtpsdr", "psdr", class(out))
  invisible(out)
}

# ------------------------------------------------------------------------------
# print / summary methods for rtpsdr (consistent with psdr)
# ------------------------------------------------------------------------------

#' @export
print.rtpsdr <- function(x, ...) {
  cat("--- Real-time Principal SDR (linear) ---\n")
  cat("Loss:", x$loss,
      " | N:", x$N,
      " | p:", if (!is.null(x$x)) ncol(x$x) else NA,
      " | Mode:", if (!is.null(x$fit$mode)) x$fit$mode else "realtime", "\n")
  if (!is.null(x$fit)) {
    cat("H:", x$fit$H,
        " | lambda:", x$fit$lambda, "\n")
  }
  if (!is.null(x$evalues)) {
    topk <- head(round(x$evalues, 4), 5)
    cat("Eigenvalues (first 5): ", paste(topk, collapse = ", "), "\n\n", sep = "")
  }
  if (!is.null(x$evectors)) {
    cat("Eigenvectors (columns are SDR directions):\n")
    print(round(x$evectors, 4))
  }
  invisible(x)
}

#' @export
summary.rtpsdr <- function(object, digits = 4, ...) {
  cat("=== Summary of rtpsdr Object ===\n")
  cat("Loss:", object$loss,
      " | N:", object$N,
      " | p:", if (!is.null(object$x)) ncol(object$x) else NA,
      " | Mode:", if (!is.null(object$fit$mode)) object$fit$mode else "realtime", "\n")
  if (!is.null(object$fit)) {
    cat("H:", object$fit$H,
        " | lambda:", object$fit$lambda, "\n")
    if (!is.null(object$fit$cutpoints)) {
      cat("cutpoints:", paste(round(object$fit$cutpoints, 3), collapse = ", "), "\n")
    }
    if (!is.null(object$fit$weight_cutpoints) && object$loss == "wlssvm") {
      cat("weight_cutpoints:", paste(round(object$fit$weight_cutpoints, 3), collapse = ", "), "\n")
    }
  }
  cat("\n--- Eigen Decomposition of Working Matrix (M) ---\n")
  cat("Top eigenvalues (up to 10):\n")
  print(round(head(object$evalues, 10), digits))
  cat("\nEstimated Eigenvectors (columns = central subspace basis):\n")
  print(round(object$evectors, digits))

  invisible(object)
}


