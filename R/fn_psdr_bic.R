#' Structural dimension selection for principal SDR
#' @description
#'
#' This function selects the structural dimension d of a fitted psdr model
#' using the BIC-type criterion proposed by Li, Artemiou and Li (2011).
#' The criterion evaluates cumulative eigenvalues of the working matrix,
#' applying a penalty that depends on the tuning parameter rho and the sample size.
#'
#' Selects the structural dimension \eqn{d} of a principal SDR model using
#' the BIC-type criterion of Li et al. (2011):
#'
#' \deqn{
#'   G(d) = \sum_{j=1}^{d} v_j \;-\;
#'   \rho \frac{d \log n}{\sqrt{n}} \, v_1 ,
#' }
#'
#' where \eqn{v_j} are the eigenvalues of the working matrix M.
#'
#' To improve robustness, cross-validation is used to choose \eqn{\rho}
#' based on the stability of the selected structural dimension across folds.
#' Specifically, for each candidate \eqn{\rho}, the data are split into
#' \eqn{K} folds, and a dimension estimate
#' \eqn{\hat{d}^{(k)}(\rho)} is obtained from fold \eqn{k}.
#'
#' The CV stability metric is defined as
#'
#' \deqn{
#'   \mathrm{Var}_{CV}(\rho)
#'   = \frac{1}{K} \sum_{k=1}^{K}
#'     \left\{ \hat{d}^{(k)}(\rho)
#'     - \overline{d}(\rho) \right\}^{2},
#' }
#'
#' where
#'
#' \deqn{
#'   \overline{d}(\rho) = \frac{1}{K} \sum_{k=1}^{K}
#'   \hat{d}^{(k)}(\rho).
#' }
#'
#' The value of \eqn{\rho} that minimizes
#' \eqn{\mathrm{Var}_{CV}(\rho)} is selected, yielding a dimension estimate that
#' is both theoretically justified (via the BIC-type criterion) and empirically
#' stable (via cross-validation).
#'
#' The function returns the selected \eqn{\rho}, the corresponding estimated
#' dimension \eqn{d}, the matrix of BIC-type criterion values, and the CV-based
#' stability metrics.
#'
#' @param obj A fitted \code{psdr} object.
#' @param rho_grid Numeric vector of candidate \eqn{\rho} values. Default \code{seq(0.001, 0.05, length=10)}.
#' @param cv_folds Number of cross-validation folds for stability evaluation. Default is 5.
#' @param plot Logical; if TRUE, plots the BIC-type criterion curve and CV stability.
#' @param seed Random seed for reproducibility.
#' @param ... Additional graphical arguments for plot.
#' @return A list of class \code{"psdr_bic"} containing:
#' \itemize{
#'   \item \code{rho_star} - selected rho minimizing cross-validated variation
#'   \item \code{d_hat} - estimated structural dimension
#'   \item \code{G_values} - matrix of BIC-type scores for each rho
#'   \item \code{cv_variation} - variation (variance) of d_hat across folds
#'   \item \code{fold_dhat} - per-fold estimated dimensions
#' }
#'@references  Li, B., Artemiou, A. and Li, L. (2011)
#' \emph{Principal support vector machines for linear and
#' nonlinear sufficient dimension reduction, Annals of Statistics 39(6): 3182–3210}.
#' @author Jungmin Shin, \email{c16267@gmail.com}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#'@seealso \code{\link{psdr}}
#'@examples
#'\donttest{
#' set.seed(1)
#' n <- 200; p <- 5;
#' x <- matrix(rnorm(n*p), n, p)
#' y <- x[,1]/(0.5+(x[,2]+1)^2)+0.2*rnorm(n)
#' fit <- psdr(x, y, loss="svm")
#' bic_out <- psdr_bic(fit, rho_grid=seq(0.05, 0.1, length=5), cv_folds=5)
#' bic_out$d_hat
#' }
#'@import stats graphics
#'@export psdr_bic

psdr_bic <- function(obj, rho_grid = seq(0.001, 0.05, length = 10), cv_folds = 5, plot = TRUE, seed = 123, ...) {

  if (!inherits(obj, "psdr"))
    stop("psdr_bic() requires a fitted psdr object.")

  set.seed(seed)

  n <- nrow(obj$x)
  v <- obj$evalues

  compute_G <- function(v, rho, n) {
    d_seq <- seq_len(length(v))
    sapply(d_seq, function(d) sum(v[1:d]) - rho * (d * log(n) / sqrt(n)) * v[1])
  }

  fold_ids <- sample(rep(1:cv_folds, length.out = n))
  fold_dhat <- matrix(NA_integer_, nrow = cv_folds, ncol = length(rho_grid))
  colnames(fold_dhat) <- paste0("rho=", format(rho_grid, digits = 3))

  for (f in seq_len(cv_folds)) {
    idx_train <- which(fold_ids != f)
    x_train <- obj$x[idx_train, , drop = FALSE]
    y_train <- obj$y[idx_train]

    fit_fold <- psdr(x_train, y_train, loss = obj$loss,
                     h = length(obj$fit$per_slice$slice),
                     lambda = obj$fit$lambda, plot = FALSE)

    v_fold <- fit_fold$evalues
    for (i in seq_along(rho_grid)) {
      G_vals <- compute_G(v_fold, rho_grid[i], length(y_train))
      fold_dhat[f, i] <- which.max(G_vals)
    }
  }

  cv_variation <- apply(fold_dhat, 2, var, na.rm = TRUE)
  rho_star <- rho_grid[which.min(cv_variation)]

  G_final <- compute_G(v, rho_star, n)
  d_hat <- which.max(G_final)

  result <- list(
    rho_star = rho_star,
    d_hat = d_hat,
    G_values = sapply(rho_grid, function(r) compute_G(v, r, n)),
    cv_variation = cv_variation,
    fold_dhat = fold_dhat
  )
  class(result) <- "psdr_bic"

  if (plot) plot(result, ...)
  return(result)
}


#' @noRd
#' @export
plot.psdr_bic <- function(x, ...) {
  v_rho <- x$rho_star
  d_len <- length(as.vector(x$G_values[, 1]))
  n_rho <- length(x$cv_variation)

  if (is.matrix(x$G_values)) {
    idx_star <- which.min(abs(as.numeric(sub("rho=", "", colnames(x$fold_dhat))) - v_rho))
    G_curve <- x$G_values[, idx_star]
  } else {
    G_curve <- x$G_values
  }

  par(mfrow = c(1, 1))
  plot(1:length(G_curve), G_curve, type = "l", lwd = 2, col = "black",
       xlab = expression(d), ylab = "BIC-type criterion", ...)
  abline(v = x$d_hat, col = "red", lty = 2, lwd = 2)
  title(main = bquote("BIC-type curve at optimal " ~ rho ~ " from CV"))
  grid()

  invisible(x)
}
