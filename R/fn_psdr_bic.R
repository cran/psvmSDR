#'Order estimation via BIC-type criterion
#' @description
#'Estimation of a structural dimensionality. Choose the k which maximizes a BIC (Bayesian information criterion) value.
#'@param obj The psdr object
#'@param rho Parameter for BIC criterion. Default is 0.01.
#'@param plot Boolean. If TRUE, the plot of BIC values are depicted.
#'@param ... Additional arguments to be passed to generic \code{plot} function.
#'@return Estimated BIC scores for determining the optimal structural dimension will be returned with plot.
#'@references  Li, B., Artemiou, A. and Li, L. (2011)
#' \emph{Principal support vector machines for linear and
#' nonlinear sufficient dimension reduction, Annals of Statistics 39(6): 3182–3210}.
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#'@seealso \code{\link{psdr}}
#'@examples
#'\donttest{
#'set.seed(1234)
#'n <- 200; p <- 10;
#'x <- matrix(rnorm(n*p, 0, 1), n, p)
#'y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + rnorm(n, 0, .2)
#'obj <- psdr(x, y, loss="svm")
#'d.hat <- psdr_bic(obj)
#'print(d.hat)
#'}
#'
#'@import stats graphics
#'@export psdr_bic


psdr_bic <- function(obj, rho=0.01, plot=TRUE, ...){
  p <- nrow(obj$evectors)
  v <- obj$evalues
  n <- nrow(obj$x)
  temp_Mat <- diag(p); temp_Mat[lower.tri(temp_Mat)] <- 1
  temp_vec <- c(1:p)
  value <- as.numeric(temp_Mat %*% v - temp_vec * (rho * log(n) / n^(1/2) * v[1]))
  class(value) <- "Gn"
  #structure(class = "Gn", value)

  if(plot == TRUE){
    plot.Gn(value, ...)
    return(value)
  }else if(plot == FALSE){
    return(value)
  }
}


#' @noRd
#' @export
print.Gn <- function(x, ...) {
  attr(x, "class") <- NULL
  obj <- x
  #writeLines("BIC result:")
  print(obj, ...)
  invisible(obj)
}


#' @noRd
#' @export
plot.Gn <- function(x, ...){
  obj <- x
  x.vec <- seq(1, length(obj), by=1)
  graphics::plot(x.vec, obj, xlab=expression(d), ylab="BIC", type='b', ...)
  abline(v=which.max(obj)-0.01, lty=2, col='red1', lwd=2)
  grid(nx = NULL, ny = NULL, lty = 1, col = "grey75", lwd = 1)
}

