#'Reconstruct the estimated sufficient predictors for a given data matrix
#'@description
#'Returning the estimated sufficient predictors \eqn{\hat{\phi}(\mathbf{x})} for a given \eqn{\mathbf{x}}
#'@param object The object from function \code{npsdr}
#'@param newdata new data \eqn{\mathbf{X}}
#'@param d structural dimensionality. d=2 is default.
#'@return the value of the estimated nonlinear mapping \eqn{\phi(\cdot)} is applied to
#' newdata \eqn{X} with dimension d is returned.
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
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
#'npsdr_x(object=obj_kernel, newdata=new.x)
#'}
#'@import stats graphics
#'@export npsdr_x

npsdr_x <- function(object, newdata, d = 2){
  obj <- object
  if (!inherits(obj, "npsdr"))
    stop("use only with \"npsdr\" objects")
  x.obj <- obj$obj.psi$scaled.x
  m <- attr(x.obj, "scaled:center")
  s <- attr(x.obj, "scaled:scale")
  bw <- obj$obj.psi$bw

  w <- obj$obj.psi$w
  l <- obj$obj.psi$l
  k <- length(l)

  v <- obj$evectors

  p <- length(m)
  n <- nrow(x.obj)
  x <- matrix(x.obj, n, p)

  n.new <- length(newdata)/p

  new.x <- matrix(newdata, n.new, p)
  new.x <- t((t(new.x) - m)/s)
  Kern <- kernel.function(new.x, y = x, param.kernel = bw)
  #Kern <- svmpath::radial.kernel(new.x, x, bw) # n * n.new
  Kern.c <- t(Kern - apply(Kern, 1, mean))
  f <- matrix(0, n.new,  k)
  for (j in 1:nrow(new.x)) {
    f[j,] <- crossprod(rep(1, n), w * Kern.c[,j])/l
  }
  pred <- f %*% v[,1:d, drop = F]
  #class(pred) <- "npsdr"
  return(pred)
}



