#'Real time sufficient dimension reduction through principal least squares SVM
#'@description
#'In stream data, where we need to constantly update the estimation as new data are collected,
#'the use of all available data can create computational challenges even for computationally efficient algorithms.
#'Therefore it is important to develop real time SDR algorithms that work efficiently in the case that there are data streams.
#'After getting an initial estimator with the currently available data,
#'the basic idea of real-time method is to update the estimator efficiently as new data are collected.
#'This function realizes real time least squares SVM SDR method for a both regression and classification problem
#'It is efficient algorithms for either adding new data or removing old data are provided.
#'@param x x in new data
#'@param y y in new data, y is continuous
#'@param obj the latest output object from the \code{rtpsdr}
#'@param h a number of slices. default is set to 10.
#'@param lambda hyperparameter for the loss function. default is set to 1.
#' @author Jungmin Shin, \email{jungminshin@korea.ac.kr}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
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
#' @return An object with S3 class "rtpsdr". Details are listed below.
#' \item{\code{x}}{input data matrix}
#' \item{\code{y}}{iniput response vector}
#' \item{\code{Mn}}{The estimated working matrix, which is obtained by the cumulative
#' outer product of the estimated parameters over H}
#' \item{\code{evalues}}{Eigenvalues of the Mn}
#' \item{\code{evectors}}{Eigenvectors of the Mn, the first d leading eigenvectors consists
#' the basis of the central subspace}
#' \item{\code{N}}{total number of observation \eqn{n_1 + n_2}}
#' \item{\code{Xbar}}{mean of total \eqn{\mathbf{x}}}
#' \item{\code{r}}{updated estimated coefficients matrix}
#' \item{\code{A}}{new A part for update. See Artemiou et. al., (2021)}
#'@seealso \code{\link{psdr}}, \code{\link{npsdr}}
#'@examples
#'\donttest{
#'p <- 5
#'m <- 500 # batch size
#'N <- 10  # number of batches
#'obj <- NULL
#'for (iter in 1:N){
#'  set.seed(iter)
#'  x <- matrix(rnorm(m*p), m, p)
#'  y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2 * rnorm(m)
#'  obj <- rtpsdr(x = x, y = y, obj=obj)
#'}
#'print(obj)
#'}
#'@import stats
#'@export rtpsdr

# -----------------------------------------------------------------------------------------------------------
# Real time Principal (weighted) least squares SDR
# -----------------------------------------------------------------------------------------------------------

rtpsdr <- function(x, y, obj=NULL, h=10, lambda=1)
{
  if(is.null(obj) == TRUE ){
    if(sum(unique(y)) == 0){
      obj_rtpsdr <- psdr(x, y, loss = "wlssvm", h, lambda)
    }else{
      obj_rtpsdr <- psdr(x, y, loss="lssvm", h, lambda)
    }
    class(obj_rtpsdr) <- "psdr"
  }else{
    obj_rtpsdr <- rtpsdr.update(x, y, obj, h, lambda)
  }
  class(obj_rtpsdr) <- "psdr"
  return(obj_rtpsdr)
}

#' @noRd
rtpsdr.update <- function(x, y, obj, h=10, lambda=1)
{
  tmp_obj <- obj
  X <- tmp_obj$x #old data: X, new data: x
  A <- tmp_obj$A
  r <- tmp_obj$r
  n <- tmp_obj$N
  Xbar <- tmp_obj$Xbar


  p <- ncol(x) #new data
  m <- nrow(x) #new data
  H <- h

  if(sum(unique(y)) != 0){#plssvm
    new.Xbar <- apply(x, 2, mean)
    Xbar <- (new.Xbar*m+Xbar*n)/(n+m)
    x <- cbind(t(t(x)-Xbar), -1)
    cov.new.X <- cov(x)*(m-1)/m
    qprob <- (1:(H-1))/H
    qy <- stats::quantile(y, qprob)

    A.new <- vector(mode = "list", length = length(qprob))
    r.new <- matrix(0, ncol=p+1, nrow=length(qprob))


    for (s in 1:length(qprob)){
      y.tilde <- rep(1, nrow(x))
      y.tilde[y < qy[s]] <- -1  #s #new y
      B_new <- m*cov.new.X/lambda + t(x) %*%  x
      C_new <- t(x) %*% y.tilde
      s_part <- diag(rep(1, p+1))+A[[s]]%*%B_new
      K <- diag(rep(1, p+1))-A[[s]] %*% B_new %*% solve(s_part) #s^w

      r.new[s,] <- t(K %*% (r[s,] + A[[s]] %*% C_new))
      A.new[[s]] <- K %*% A[[s]]
    }
  }else{ #pwlssvm
    Xbar <- (apply(x, 2, mean)*m+Xbar*n)/(n+m)
    x <- cbind(t(t(x)-Xbar), -1)
    weight_list <- seq(0, 1, length=H+2)[2:(H+1)]
    A.new <- vector(mode = "list", length = H)
    r.new <- matrix(0, ncol=p+1, nrow=H)

    cov.new.X <- cov(x)*(m-1)/m
    for (i in 1:H){
      new.W <- diag(ifelse(y==1, 1-weight_list[i], weight_list[i]))
      B_new <- m*cov.new.X/lambda + t(x) %*% new.W %*% x
      C_new <- t(x) %*% new.W %*% y
      s_part <- diag(rep(1, p+1))+A[[i]]%*%B_new
      s <- diag(rep(1, p+1))-A[[i]] %*% B_new %*% solve(s_part)

      r.new[i,] <- t(s %*% (r[i,] + A[[i]] %*% C_new))
      A.new[[i]] <- s %*% A[[i]]
    }
  }
  Mn <- t(r[,1:p])%*%r[,1:p]
  eigen.Mn <- eigen(Mn)
  tmp <- list("M"=Mn, "evalues" = eigen.Mn$values, "evectors" = eigen.Mn$vectors,
                r=r.new, A=A.new, "x" = x[ ,-ncol(x)], "y" = y, N=n+m, Xbar=Xbar, "loss"=obj$loss )
  class(tmp) <- "psdr"
  return(tmp)
}




