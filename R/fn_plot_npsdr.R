#' Scatter plot with sufficient predictors from npsdr() function
#' @param x object from the function \code{npsdr()}
#' @param d number of sufficient predictors. Default is 1.
#' @param lowess draw a lowess curve. Default is TRUE.
#' @param ... Additional arguments to be passed to generic \code{plot} function.
#' @return A scatter plot with sufficient predictors.
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#' @seealso \code{\link{npsdr_x}}, \code{\link{npsdr}}
#' @import stats
#' @examples
#'\donttest{
#'set.seed(1)
#'n <- 200;
#'p <- 5;
#'x <- matrix(rnorm(n*p, 0, 2), n, p)
#'y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
#'obj_kernel <- npsdr(x, y, plot=FALSE)
#'plot(obj_kernel)
#' }
#' @import stats graphics
#' @export

plot.npsdr <- function(x, ..., d=1, lowess=TRUE) {
  object <- obj <- x
  dim <- d
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (!inherits(obj, "npsdr"))
    stop("use only with \"npsdr\" objects")
  temp <- object$evectors
  if(sum(unique(object$y)) != 0){
    #obj_psdr <- object$x %*% temp
    obj_npsdr <- phix(object$x, object)
    if(d <= 2){
      par(mfrow=c(1,dim))
    }else{
      par(mfrow=c(ceiling(sqrt(dim)), ceiling(sqrt(dim))))
    }
    for(j in 1:dim){
      if(lowess == TRUE){
        graphics::plot(obj_npsdr[,j], object$y, xlab = paste("Sufficient predictor ",j), ylab  = expression(Y), pch=16, ...)
        graphics::lines(lowess(obj_npsdr[,j], object$y), col='red', lwd=1)
        grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
      }else{
        graphics::plot(obj_npsdr[,j], object$y, xlab = paste("Sufficient predictor ",j), ylab  = expression(Y), pch=16, ...)
        grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
      }
    }
    #par(mfrow=c(1,1))
  }else{
    x.nlsvm <- phix(object$x, object)
    par(mar=c(5,5,5,5), oma=c(1,1,1,1))
    graphics::plot(x.nlsvm[,1], x.nlsvm[,2], type = "n", xlab = paste("Sufficient predictor ",1),
                   ylab = paste("Sufficient predictor ",2), ...)
    graphics::points(x.nlsvm[object$y == 1,1], x.nlsvm[object$y == 1,2], col = 4, pch = 16, ...)
    graphics::points(x.nlsvm[object$y != 1,1], x.nlsvm[object$y != 1,2], col = 2, pch = 16, ...)
    grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
  }
}


