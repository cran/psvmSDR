#'A unified Principal sufficient dimension reduction method via kernel trick
#'@description
#'Principal Sufficient Dimension Reduction method
#'@param x data matrix
#'@param y either continuous or (+1,-1) typed binary response vector
#'@param loss pre-specified loss functions belongs to \code{"svm", "logit", "l2svm", "wsvm", "qr", "asls", "wlogit", "wl2svm", "lssvm", "wlssvm"}, and user-defined loss function object also can be used formed by inside double (or single) quotation mark. Default is 'svm'.
#'@param h the number of slices. default value is 10
#'@param lambda hyperparameter for the loss function. default value is 1
#'@param b number of basis functions for a kernel trick, floor(length(y)/3) is default
#'@param eps threshold for stopping iteration with respect to the magnitude of derivative, default value is 1.0e-4
#'@param max.iter maximum iteration number for the optimization process. default value is 30
#'@param eta learning rate for gradient descent method. default value is 0.1
#' @param mtype type of margin, either "m" or "r" refer margin and residual, respectively (See, Table 1 in the pacakge manuscript). When one use user-defined loss function this argument should be specified. Default is "m".
#' @param plot If \code{TRUE} then it produces scatter plots of \eqn{Y} versus the first sufficient predictor. The default is FALSE.
#' @return An object with S3 class "npsdr". Details are listed below.
#' \item{\code{evalues}}{Eigenvalues of the estimated working matrix M.}
#' \item{\code{evectors}}{Eigenvectors of the estimated working matrix M, the first d leading eigenvectors consists
#' the basis of the central subspace.}
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
#'@seealso \code{\link{npsdr_x}}, \code{\link{psdr}}, \code{\link{rtpsdr}}
#'@examples
#'\donttest{
#'set.seed(1)
#'n <- 200;
#'p <- 5;
#'x <- matrix(rnorm(n*p, 0, 2), n, p)
#'y <- 0.5*sqrt((x[,1]^2+x[,2]^2))*(log(x[,1]^2+x[,2]^2))+ 0.2*rnorm(n)
#'obj_kernel <- npsdr(x, y, plot=FALSE)
#'print(obj_kernel)
#'plot(obj_kernel)
#'}
#'@import stats graphics
#'@importFrom graphics abline grid par
#'@importFrom stats cov dist rnorm
#'@export npsdr



npsdr <- function(x, y, loss="svm", h=10, lambda=1, b=floor(length(y)/3),
                  eps=1.0e-5, max.iter=100, eta=0.1, mtype, plot=TRUE)
{
  if(sum(as.character(loss) == c("lssvm", "wlssvm")) == 0){
    if(!is.matrix(x) & !is.data.frame(x))
      stop("x must be a matrix or dataframe.")
    if(ncol(as.matrix(y)) != 1)
      stop("y must be a univariate.")
    if(is.null(eta) == T ){
      #message("eta is a learning rate which should be specified as a positve value. It is set to 0.1")
      eta <- 0.1
    }
    if(is.null(loss)==T){
      #message("Loss function should be specified, hinge loss was applied.")
      loss <- 'svm'
    }

    if(is.null(max.iter) == T){
      #message("max.iter is set to 100 as a default.")
      max.iter <- 100
    }
  }
  if(is.null(lambda) == T){
    #message("lambda is set to 1 as a default.")
    lambda <- 1
  }
  if (length(y) != nrow(x)){     #check lengths
    stop("The response and predictors have different number of observations.")
  }

  psi.gen <- get.psi(x, y, b)
  Psi.new <- psi.gen$w   #n*k
  n <- nrow(Psi.new)
  p <- ncol(Psi.new)
  H <- h
  init.theta <- rnorm(sd=1,n=p)
  step <- 1/H
  pi.grid <- seq(step, 1-step, by = step)

  # generate y.tilde

  qprob <- (1:(H-1))/H
  qy <- stats::quantile(y, qprob)

  theta.new <- rep(0,p)
  w.init <- matrix(init.theta, nrow=p, ncol=length(qprob))
  w.final <- matrix(0, nrow=p, ncol=length(qprob))
  eigen.mat <- diag(1,p,p)

  type.list <- c("svm","logit","l2svm", "wsvm", "qr","asls", "wlogit","wl2svm")
  type.list2 <- c("lssvm","wlssvm")

  if(sum(as.character(loss) == type.list) != 0){
    if(as.character(loss) == "svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- -Psi[,k]*y.tilde*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*(A[k,]%*%w[,s])  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
          # if(max(abs(deriv)) < eps)
          #   break

        }
        w.final[,s] <- w[,s] #s,s
      }
    }
    if(as.character(loss) == "l2svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- 2*(1-margin.v)*(-Psi[,k]*y.tilde)*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "wl2svm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      for (s in 1:length(qprob)) {
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.bi <- y.new
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- 2*(1-margin.v)*(-Psi[,k]*y.bi)*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }

          w[,s] <- theta.new  #c
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "logit"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.tilde #s
            deriv <- -Psi[,k]*y.tilde*(1/(1+exp(margin.v)))#k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }

    if(as.character(loss) == "wlogit"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      for (s in 1:length(qprob)){
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.bi <- y.new
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- weight*(-Psi[,k])*y.bi*(1/(1+exp(margin.v)))   #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }

    if(as.character(loss) == "wsvm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      y.new <- y
      for (s in 1:length(qprob)) {
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.bi <- y.new
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi

          for (k in 1:p){
            margin.v <- Psi %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -weight*Psi[,k]*y.bi*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*(1/nrow(Psi))*A[k,]%*%w[,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }

          w[,s] <- theta.new  #c
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "qr"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(pi.grid)) {
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.new <- y
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            u <- y.new - Psi%*%w[,s]   #s
            #derivative.j <- 2*w[k,s]+lambda*sum((-z[,k]*pi.grid[s])*as.numeric(I(u>0))+(z[,k]*(1-pi.grid[s]))*as.numeric(I(u<0)))*(1/length(y.new)) #k,s,k,s,k,s
            derivative.j <- 2*(1/nrow(Psi))*w[k,s]+lambda*(1/length(y.new))*sum(-Psi[,k]*{pi.grid[s]*as.numeric(I(u>0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          #print(derivative.j)
          w[,s] <- theta.new  #s
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }
    if(as.character(loss) == "asls"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(pi.grid)) {
        for(iter in 1:max.iter){
          Psi <- Psi.new
          y.new <- y
          n <- nrow(Psi)
          w <- w.init
          A <- t(Psi)%*%Psi
          for (k in 1:p){
            u <- y.new - Psi%*%w[,s]   #s
            derivative.j <- 2*(1/nrow(Psi))*w[k,s]+lambda*(1/length(y.new))*sum(-Psi[,k]*2*u*{pi.grid[s]*as.numeric(I(u>0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          #print(derivative.j)
          w[,s] <- theta.new  #s
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }

    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + w.final[,h, drop = F] %*% t(w.final[,h, drop = F])
    result <- eigen(Mn)
    v <- result$vectors
    u <- result$values
    newlist <- list("x" = x, "y" = y, "Mn"=Mn, "evalues" = u, "evectors" = v, "obj.psi" = psi.gen)
    class(newlist) <- c("npsdr", class(newlist))
    structure(class = "npsdr", newlist)


    if(plot == TRUE){
      plot.npsdr(newlist, d=1, lowess=FALSE)
      return(newlist)
    }else{
      return(newlist)
    }
  }

  else if(sum(as.character(loss) == type.list) == 0 & sum(as.character(loss)==type.list2)==0){
    ft <- E(loss)
    # grid.m <- seq(-2,2,length=100)
    # plot(grid.m, ft(grid.m,a,c), type="l", xlab="margin", ylab="loss")
    # message("loss function must be a convex function")
    interval <- 1.0e-5

    ## arbitrary loss continuous##

    if(sum(unique(y)) != 0){
      ## main for loop ##
      for(s in 1:length(qprob)){
        y.tilde.new <- rep(1, nrow(Psi.new))
        y.tilde.new[y < qy[s]] <- -1  #s
        w <- w.init
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.tilde <- y.tilde.new
          n <- nrow(Psi)
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            star <- (fn_arbitrary_nonlinear_loss(Psi, y.tilde,theta=w[,s]+interval*eigen.mat[,k],lambda,loss,prob=pi.grid[s], mtype="m")-
                       fn_arbitrary_nonlinear_loss(Psi, y.tilde,theta=w[,s]-interval*eigen.mat[,k],lambda,loss,prob=pi.grid[s], mtype="m"))
            derivative.vec[k] <- sign(star)*exp( log(abs(star))-log(2*interval))
            theta.new[k] <- w[k,s] - eta * derivative.vec[k] ###k,k,s,k
          }
          w[,s] <- theta.new   #s
          if(max(abs(derivative.vec)) < eps)
          break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s] ##s,s
      }
    }

    ## arbitrary loss binary##
    if(sum(unique(y))==0){
      y.new <- y
      for(s in 1:length(qprob)){
        w <- w.init
        for(iter in 1:max.iter){
            Psi <- Psi.new
            y.bi <- y.new
          n <- nrow(Psi)
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            derivative.vec[k] <- (fn_arbitrary_nonlinear_binary_loss(Psi, y.bi, prob=pi.grid[s],theta=w[,s]+interval*eigen.mat[,k],lambda,loss, mtype="m")-
                                    fn_arbitrary_nonlinear_binary_loss(Psi, y.bi ,prob=pi.grid[s],theta=w[,s]-interval*eigen.mat[,k],lambda,loss, mtype="m")) / (2*interval)
            theta.new[k] <- w[k,s] - eta * derivative.vec[k]
          }
          w[,s] <- theta.new
          if(mean(abs(derivative.vec)) < eps)
          break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s]
      }
    }

    Working_mat <- matrix(0, p, p)
    for (h in 1:length(qprob)) Working_mat <- Working_mat + w.final[,h, drop = F] %*% t(w.final[,h, drop = F])
    result <- eigen(Working_mat)
    v <- result$vectors
    u <- result$values
    newlist <- list("x" = x, "y" = y, "Mn"=Working_mat, "evalues" = u, "evectors" = v, "obj.psi" = psi.gen)
    class(newlist) <- c("npsdr", class(newlist))
    structure(class = "npsdr", newlist)

    if(plot == TRUE){
      plot.npsdr(newlist, d=1, lowess=FALSE)
      return(newlist)
    }else{
      return(newlist)
    }
  }

  else if(sum(as.character(loss) == type.list2) != 0){
    weight_list <- seq(0, 1, length=H+1)[2:(H)]
    bar.x <- apply(Psi.new, 2, mean)
    x.star <- cbind(t(t(Psi.new)-bar.x), -1)
    cov.x.star <- cov(x.star)
    if(as.character(loss) == "lssvm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      # estimate bases
      r.H <- matrix(0, ncol=p+1, nrow=H-1)
      for (s in 1:length(qprob)){
        y.tilde <- rep(1, nrow(Psi.new))
        y.tilde[y < qy[s]] <- -1  #s
        M <- solve(n * cov.x.star/lambda+t(x.star) %*% x.star)
        C <- t(x.star) %*%  y.tilde
        r.H[s,] <- M %*% C
      }
    }
    if(as.character(loss)=="wlssvm"){
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      #y.binary <- y
      for (pi in weight_list){
        W <- diag(c(ifelse(y==1, 1-pi, pi)))
        M <- solve(n*cov.x.star/lambda + t(x.star) %*% W %*% x.star)
        C <- t(x.star) %*% W %*% y
        r.H[which(weight_list==pi),] <- M %*% C
      }
    }
    Working_mat <- eigen(t(r.H[,1:p])%*%r.H[,1:p])
    v <- Working_mat$vectors
    u <- Working_mat$values
    newlist <- list("x" = x, "y" = y, "Mn"=Working_mat, "evalues" = u, "evectors" = v, "obj.psi" = psi.gen)

    class(newlist) <- c("npsdr", class(newlist))
    structure(class = "npsdr", newlist)
    if(plot == TRUE){
      plot(newlist, d=1, lowess=FALSE)
      return(newlist)
    }else{
      return(newlist)
    }
  }
}

#' @noRd
#' @export
print.npsdr <- function(x, ...) {
  obj <- x
  d <- list(evalues = obj$evalues, evectors= obj$evectors)
  print(d, ...)
  invisible(d)
}



#' @noRd
get.psi <- function(x, y, b=floor(length(y)/3)) {
  n <- nrow(x)
  x <- scale(x)
  bw <- 1/mean(as.numeric(stats::dist(x)))^2 # bw parameter for kernel
  Kn <- kernel.function(x, y = x, param.kernel = bw)
  #Kn <- svmpath::radial.kernel(x, x, bw)
  Qn <- diag(n) - matrix(1/n, n, n)

  eigen.psi <- eigen(Qn %*% Kn %*% Qn)
  Psi.new <- eigen.psi$vectors[,1:b, drop = F] # Psi
  l <- eigen.psi$values[1:b]
  tmp.obj <- list("w"=Psi.new, "l"=l, "scaled.x"= x, "bw" = bw, "b" = b)
  tmp.obj
  #class(tmp.obj) <- "npsdr"
}


#'@noRd
phix <- function(value, object) {
  psi.function <- psi.function
  x <- value
  v <- object$evector
  w <- object$obj.psi$w
  l <- object$obj.psi$l
  d <- p <- dim(x)[2]
  kernel.function <- kernel.function(x, y=x, param.kernel = 1/p)
  tau <- mean(as.numeric(dist(x)))
  kernel.param <- 1/tau^2
  p <- ncol(x)
  if (length(value) == p) {
    temp <- psi.function(value, x, v[,1:d, drop = F], w, l, kernel.function, kernel.param)
  } else if (ncol(value) == p) {
    temp <- t(apply(value, 1, psi.function, x, v[,1:d, drop = F], w, l, kernel.function, kernel.param))
  } else if (nrow(value) == p) {
    temp <- t(apply(value, 2, psi.function, x, v[,1:d, drop = F], w, l, kernel.function, kernel.param))
  } else stop("check `str(value)`")
  temp
}


#'@noRd
psi.function <- function(value, x, v, w, l, kernel.function, kernel.param){
  value <- matrix(value, 1, length(value))
  temp <- kernel.function(value, x, kernel.param)
  psi.value <- apply(w * c(temp - mean(temp)), 2, sum)/l
  rslt <- psi.value %*% v
  class(rslt) <- "npsdr"
  return(rslt)
}


#'@noRd
kernel.function <- function (x, y = x, param.kernel = 1/p) {
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  exp(-a * param.kernel)
}


