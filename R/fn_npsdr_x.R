#' Reconstruct estimated sufficient predictors for new data
#'
#' @description
#' Computes the nonlinear sufficient predictors \eqn{\hat{\phi}(\mathbf{x})} for a
#' new data matrix using a previously fitted \code{npsdr} object.
#'
#' This function evaluates the learned kernel-based sufficient dimension
#' reduction (SDR) mapping on new observations. Given a fitted nonlinear SDR
#' model \eqn{\hat{\phi}} estimated from \code{npsdr()}, the function computes:
#'
#' \deqn{
#'   \hat{Z} = \hat{\phi}(X_{\text{new}}) = \Psi(X_{\text{new}})^{\top}
#'   \, \hat{V}_{1:d},
#' }
#'
#' where \eqn{\Psi(\cdot)} is the kernel feature map constructed from the training
#' data, and \eqn{\hat{V}_{1:d}} contains the first \eqn{d} eigenvectors of the
#' estimated working matrix \eqn{M}. These eigenvectors span the estimated
#' central subspace in the kernel-transformed space.
#'
#' This enables users to extract sufficient predictors for downstream tasks
#' such as visualization, classification, regression, or clustering on new data,
#' without re-estimating the SDR model.
#'
#'@param object The object from function \code{npsdr}
#'@param newdata new data \eqn{\mathbf{X}}
#'@param d structural dimensionality. d=2 is default.
#'@return the value of the estimated nonlinear mapping \eqn{\phi(\cdot)} is applied to
#' newdata \eqn{X} with dimension d is returned.
#' @author Jungmin Shin, \email{c16267@gmail.com}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#'@seealso \code{\link{npsdr}}
#'@examples
#'\donttest{
#'set.seed(1)
#'n <- 200; n.new <- 300
#'p <- 5;
#'h <- 20;
#'x <- matrix(rnorm(n*p, 0, 2), n, p)
#'y <- 0.5*sqrt((x[,1]^2+x[,2]^2))*(log(x[,1]^2+x[,2]^2))+ 0.2*rnorm(n)
#'new.x <- matrix(rnorm(n.new*p, 0, 2), n.new, p)
#'obj_kernel <- npsdr(x, y)
#'z_new <- npsdr_x(object=obj_kernel, newdata=new.x)
#'dim(z_new)
#'}
#'@import stats graphics
#'@export npsdr_x


npsdr_x <- function(object, newdata, d = 2) {

  if (!inherits(object, "npsdr"))
    stop("Input must be an 'npsdr' object.")
  if (!is.matrix(newdata) && !is.data.frame(newdata))
    stop("newdata must be a matrix or data.frame.")
  if (!is.numeric(newdata))
    newdata <- as.matrix(newdata)

  psi.obj <- object$obj.psi
  if (is.null(psi.obj))
    stop("obj.psi is missing not a valid fitted npsdr model.")
  v <- object$evectors
  if (is.null(v))
    stop("evectors not found eigenvectors missing in model output.")

  x_train <- psi.obj$scaled.x
  bw <- psi.obj$bw
  w <- psi.obj$w
  l <- psi.obj$l
  m <- attr(x_train, "scaled:center")
  s <- attr(x_train, "scaled:scale")
  p <- ncol(x_train)
  n <- nrow(x_train)

  new.x <- scale(as.matrix(newdata), center = m, scale = s)
  n.new <- nrow(new.x)

  Kern <- kernel.function(new.x, y = x_train, param.kernel = bw)
  Kern.c <- t(Kern - matrix(colMeans(Kern), nrow = n.new, ncol = n, byrow = TRUE))

  f <- matrix(0, n.new, length(l))
  for (j in 1:n.new) {
    f[j, ] <- crossprod(rep(1, n), w * Kern.c[, j]) / l
  }

  d <- min(d, ncol(v))
  pred <- f %*% v[, 1:d, drop = FALSE]
  colnames(pred) <- paste0("phi", 1:d)
  rownames(pred) <- rownames(newdata)

  return(pred)
}
