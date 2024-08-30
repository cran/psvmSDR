#' Scatter plot with sufficient predictors from psdr() function
#'
#' @param x object from the function \code{psdr()}
#' @param d number of sufficient predictors. Default is 1.
#' @param lowess draw a locally weighted scatterplot smoothing curve. Default is TRUE.
#' @param ... Additional arguments to be passed to generic \code{plot} function.
#' @return A scatter plot with sufficient predictors.
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#' @seealso \code{\link{psdr_bic}}, \code{\link{psdr}}
#' @examples
#'\donttest{
#' set.seed(1)
#' n <- 200; p <- 5;
#' x <- matrix(rnorm(n*p, 0, 2), n, p)
#' y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
#' obj <- psdr(x, y)
#' plot(obj, d=2, lowess=TRUE)
#' }
#' @import stats
#' @importFrom  graphics lines points plot layout legend par arrows abline axis
#' @export


plot.psdr <- function(x, ..., d=1, lowess=TRUE) {
  object <- obj <- x
  dim <- d
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (!inherits(obj, "psdr"))
    stop("use only with \"psdr\" objects")
  temp <- object$evectors
  if(sum(unique(object$y)) != 0){
    obj_psdr <- object$x %*% temp
    if(d <= 2){
      par(mfrow=c(1,dim))
    }else{
      par(mfrow=c(ceiling(sqrt(dim)), ceiling(sqrt(dim))))
    }
    for(j in 1:dim){
      if(lowess == TRUE){
        graphics::plot(obj_psdr[,j], object$y, pch=16, xlab =  paste("Sufficient predictor ",j), ylab  = expression(Y), ...)
        graphics::lines(lowess(obj_psdr[,j], object$y), col='red', lwd=1)
        grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
      }else{
        graphics::plot(obj_psdr[,j], object$y, pch=16, xlab =  paste("Sufficient predictor ",j), ylab  = expression(Y), ...)
        grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
      }
    }
    #par(mfrow=c(1,1))
  }else{
    x.lsvm <- object$x %*% temp
    par(mar=c(5,5,5,5), oma=c(1,1,1,1))
    graphics::plot(x.lsvm[,1], x.lsvm[,2], type = "n", xlab = paste("Sufficient predictor ",1),
                   ylab = paste("Sufficient predictor ",2), ...)
    graphics::points(x.lsvm[object$y == 1,1], x.lsvm[object$y == 1,2], col = 4, pch = 16, ...)
    graphics::points(x.lsvm[object$y != 1,1], x.lsvm[object$y != 1,2], col = 2, pch = 16, ...)
    grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
  }
}



