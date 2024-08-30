#' Unified linear principal sufficient dimension reduction methods
#'
#' A function for a linear principal sufficient dimension reduction.
#'
#' Two examples of the usage of user-defined losses are presented below (\code{u} represents a margin):
#'
#' \code{mylogit <- function(u, ...) log(1+exp(-u))},
#'
#' \code{myls <- function(u ...) u^2}.
#'
#' Argument \code{u} is a function variable  (any character is possible) and the argument \code{mtype} for \code{psdr()} determines a type of a margin, either (\code{type="m"}) or (\code{type="r"}) method. \code{type="m"} is a default.
#' Users have to change \code{type="r"}, when applying residual type loss.
#' Any additional parameters of the loss can be specified via \code{...} argument.
#'
#' @param x input matrix, of dimension \code{nobs} x \code{nvars}; each row is an observation vector. Requirement: \code{nvars}>1; in other words, \code{x} should have 2 or more columns.
#' @param y response variable, either can be continuous variable or (+1,-1) coded binary response vector.
#' @param loss pre-specified loss functions belongs to "svm", "logit","l2svm","wsvm", and etc., and user-defined loss function object also can be used formed by inside double (or single) quotation mark. Default is 'svm'.
#' @param h the number of slices and probabilities equally spaced in \eqn{(0,1)}. Default value is 10.
#' @param lambda the cost parameter for the svm loss function. The default value is 1.
#' @param eps the threshold for stopping iteration with respect to the magnitude of the change of the derivative. The default value is 1.0e-5.
#' @param max.iter maximum iteration number for the optimization process. default value is 100.
#' @param eta learning rate for the gradient descent algorithm. The default value is 0.1.
#' @param mtype a margin type, which is either margin ("m") or residual ("r") (See, Table 1 in manuscript). Only need when user-defined loss is used. Default is "m".
#' @param plot If \code{TRUE} then it produces scatter plots of \eqn{Y} versus \eqn{\hat{B^{\top}}_{j}\mathbf{X}}. \eqn{j} can be specified by the user with \eqn{j=1} as a default. The default is FALSE.
#' @return An object with S3 class "psdr". Details are listed below.
#' \item{\code{Mn}}{The estimated working matrix, which is obtained by the cumulative
#' outer product of the estimated parameters over the slices. It will not print out, unless it is called manually.}
#' \item{\code{evalues}}{Eigenvalues of the working matrix \eqn{Mn}}
#' \item{\code{evectors}}{Eigenvectors of the \eqn{Mn}, the first leading \eqn{d} eigenvectors consists
#' the basis of the central subspace}
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
#' @seealso \code{\link{psdr_bic}}, \code{\link{rtpsdr}}
#'@examples
#'## ----------------------------
#'## Linear PM
#'## ----------------------------
#' set.seed(1)
#' n <- 200; p <- 5;
#' x <- matrix(rnorm(n*p, 0, 2), n, p)
#' y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
#' y.tilde <- sign(y)
#' obj <- psdr(x, y)
#' print(obj)
#' plot(obj, d=2)
#'
#'## ----------------------------
#'## Kernel PM
#'## ----------------------------
#' obj_wsvm <- psdr(x, y.tilde, loss="wsvm")
#' plot(obj_wsvm)
#'
#'## ----------------------------
#'## User-defined loss function
#'## ----------------------------
#' mylogistic <- function(u) log(1+exp(-u))
#' psdr(x, y, loss="mylogistic")
#'
#'@import stats graphics
#'@export psdr
#'
#'
psdr <- function(x, y, loss="svm", h=10, lambda=1, eps=1.0e-5, max.iter=100, eta=0.1, mtype="m", plot=FALSE){
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
      #message("max.iter is set to 30 as a default.")
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

  n <- nrow(x)
  p <- ncol(x)
  H <- h

  init <- init.theta <- rnorm(dim(x)[2],0,1)

  bar.x <- apply(x, 2, mean)
  x.star <- cbind(t(t(x)-bar.x), -1)
  cov.x.star <- cov(x.star)
  cov.x <- cov(x)
  step <- 1/H
  pi.grid <- seq(step, 1-step, by = step)

  # generate y.tilde
  qprob <- (1:(H-1))/H
  qy <- stats::quantile(y, qprob)

  # standardization
  temp <- eigen(cov.x)
  D <- diag(sqrt(temp$values))
  V <- temp$vectors
  sd.x <-  V %*% D
  inv.sd.x <- diag(1/(sqrt(temp$values))) %*% t(V)

  centered.x <- t(x) - bar.x  # p x n
  z.new <- t((inv.sd.x %*% centered.x))
  eigen.mat <- diag(1,p,p)
  theta.new <- rep(0,p)
  type.list <- c("svm","logit","l2svm", "wsvm", "qr","asls", "wlogit","wl2svm")
  type.list2 <- c("lssvm","wlssvm")

  if(sum(as.character(loss) == type.list) != 0){
    w.init <- matrix(init, nrow=p, ncol=length(qprob))
    w.final <- matrix(0, nrow=p, ncol=length(qprob))

    if(as.character(loss) == "svm"){
      if(sum(unique(y)) == 0){
        stop("response variable should be continuous!")
      }
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          z <- z.new
          y.tilde <- y.tilde.new
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.tilde #s
            deriv <- -z[,k]*y.tilde*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          # if(max(abs(deriv)) < eps)
          #   break
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "l2svm"){
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      for (s in 1:length(qprob)) {
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          z <- z.new
          y.tilde <- y.tilde.new
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.tilde #s
            deriv <- -z[,k]*y.tilde*as.numeric(I((1-margin.v)>0))*2*(1-margin.v) #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          # if(max(abs(deriv)) < eps)
          #   break
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
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for (s in 1:length(qprob)) {
        w.init <- matrix(init, nrow=p, ncol=length(qprob))
        for(iter in 1:max.iter){
          z <- z.new
          y.bi <- y.new
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -z[,k]*y.bi*as.numeric(I((1-margin.v)>0))*2*(1-margin.v)*weight #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          # if(max(abs(deriv)) < eps)
          #   break
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
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        for(iter in 1:max.iter){
          z <- z.new
          y.tilde <- y.tilde.new
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.tilde #s
            deriv <- -z[,k]*y.tilde*(1/(1+exp(margin.v)))#k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          # if(max(abs(deriv)) < eps)
          #   break
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
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for (s in 1:length(qprob)){
        w.init <- matrix(init, nrow=p, ncol=length(qprob))
        for(iter in 1:max.iter){
          z <- z.new
          y.bi <- y.new
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- weight*(-z[,k])*y.bi*(1/(1+exp(margin.v)))   #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new
          # if(max(abs(deriv)) < eps)
          #   break
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
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for (s in 1:length(qprob)) {
        w.init <- matrix(init, nrow=p, ncol=length(qprob))
        for(iter in 1:max.iter){
          z <- z.new
          y.bi <- y.new
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            margin.v <- z %*% w[,s] * y.bi #s
            weight <- (1-pi.grid[s])*as.numeric(y.bi==1)+(pi.grid[s])*(as.numeric(y.bi==-1)) #s,s
            deriv <- -weight*z[,k]*y.bi*as.numeric(I((1-margin.v)>0)) #k
            derivative.j <- lambda*mean(deriv) + 2*w[k,s]  ##k,s
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          # if(max(abs(deriv)) < eps)
          #   break
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        }
        w.final[,s] <- w[,s]
      }
    }
    if(as.character(loss) == "qr"){
      #if(sum(unique(y)) == 0)
      #  stop("response variable should be continuous!")
      for (s in 1:length(pi.grid)) {
        w.init <- matrix(init, nrow=p, ncol=length(qprob))
        for(iter in 1:max.iter){
          z <- z.new
          y.new <- y
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            u <- y.new - z%*%w[,s]   #s
            derivative.j <- 2*w[k,s]+lambda*(1/length(y))*sum(-z[,k]*{pi.grid[s]*as.numeric(I(u>0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          # if(max(abs(deriv)) < eps)
          #   break
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }
    if(as.character(loss) == "asls"){
      for (s in 1:length(pi.grid)) {
        w.init <- matrix(init, nrow=p, ncol=length(qprob))
        for(iter in 1:max.iter){
          z <- z.new
          y.new <- y
          n <- nrow(z)
          w <- w.init
          for (k in 1:p){
            u <- y.new - z%*%w[,s]
            derivative.j <- 2*w[k,s]+lambda*(1/length(y))*sum((-z[,k]*2*u)*{pi.grid[s]*as.numeric(I(u>=0))+(1-pi.grid[s])*as.numeric(I(u<=0))})
            theta.new[k] <- w[k,s] -  eta*derivative.j  ##k,k,s
          }
          w[,s] <- theta.new  #s
          # if(max(abs(deriv)) < eps)
          #   break
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
          w.init <- matrix(theta.new, nrow=p, ncol = length(pi.grid))
        }
        w.final[,s] <- w[,s] #s,s
      }
    }

    psi <- t(inv.sd.x) %*% w.final
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])
    eigen.Mn <- eigen(Mn)
    newlist <- list("loss" = loss, "x" = x, "y" = y,"M"=Mn, "evalues" = eigen.Mn$values, "evectors" = eigen.Mn$vectors)

    class(newlist) <- c("psdr", class(newlist))
    #class(newlist) <- "psdr"
    structure(class = "psdr", newlist)

    if(plot == TRUE){
      plot.psdr(newlist, d=1)
      return(newlist)
    }else{
      return(newlist)
    }
  }
  else if(sum(as.character(loss) == type.list) == 0 & sum(as.character(loss)==type.list2)==0){
    ft <- E(loss)
    # grid.m <- seq(-2,2,length=100)
    w.init <- matrix(init, nrow=p, ncol=length(qprob))
    w.final <- matrix(0, nrow=p, ncol=length(qprob))

    ## arbitrary loss continuous ##
    if(sum(unique(y)) != 0){
      for(s in 1:length(qprob)){
        y.tilde.new <- rep(1, nrow(x))
        y.tilde.new[y < qy[s]] <- -1  #s
        pos.rate <- sum(y.tilde.new==1)/nrow(z.new)
        neg.rate <- sum(y.tilde.new==-1)/nrow(z.new)
        w <- w.init
        interval <- 1.0e-5
        for(iter in 1:max.iter){
          z <- z.new
          y.tilde <- y.tilde.new
          n <- nrow(z)
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            star <- (fn_arbitrary_loss(z, y.tilde,theta=w[,s]+interval*eigen.mat[,k],lambda,loss,prob=pi.grid[s], mtype="m")-
                       fn_arbitrary_loss(z, y.tilde,theta=w[,s]-interval*eigen.mat[,k],lambda,loss,prob=pi.grid[s], mtype="m"))
            derivative.vec[k] <- sign(star)*exp( log(abs(star))-log(2*interval))
            theta.new[k] <- w[k,s] - eta * derivative.vec[k] ###k,k,s,k
          }
          w[,s] <- theta.new   #s
          # if(max(abs(deriv)) < eps)
          #   break
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s] ##s,s
      }
    }

    ## arbitrary loss binary##
    if(sum(unique(y))==0){
      y.new <- y
      pos.rate <- sum(y.new==1)/nrow(z.new)
      neg.rate <- sum(y.new==-1)/nrow(z.new)
      for(s in 1:length(qprob)){
        w <- w.init
        interval <- 1.0e-5
        for(iter in 1:max.iter){
          z <- z.new
          y.bi <- y.new
          derivative.vec <- rep(0,p)
          for (k in 1:p){
            star <- (fn_arbitrary_binary_loss(z, y.bi, prob=pi.grid[s],theta=w[,s]+interval*eigen.mat[,k],lambda,loss, mtype="m")-
                       fn_arbitrary_binary_loss(z, y.bi,prob=pi.grid[s],theta=w[,s]-interval*eigen.mat[,k],lambda,loss, mtype="m"))
            derivative.vec[k] <- sign(star)*exp( log(abs(star))-log(2*interval))
            theta.new[k] <- w[k,s] - eta * derivative.vec[k]
          }
          w[,s] <- theta.new
          # if(max(abs(deriv)) < eps)
          #   break
          if(max(abs(theta.new-w.init[,s])) < eps)
            break
        }
        w.init <- matrix(theta.new, nrow=p, ncol = length(qprob))
        w.final[,s] <- w[,s]
      }
    }
    psi <- t(inv.sd.x) %*% w.final
    Mn <- matrix(0, p, p)
    for (h in 1:length(qprob)) Mn <- Mn + psi[,h, drop = F] %*% t(psi[,h, drop = F])
    eigen.Mn <- eigen(Mn)
    newlist <- list("loss" = loss, "x" = x, "y" = y, "M"=Mn, "evalues" = eigen.Mn$values, "evectors" = eigen.Mn$vectors)

    class(newlist) <- c("psdr", class(newlist))
    #class(newlist) <- "psdr"
    structure(class = "psdr", newlist)
    if(plot == TRUE){
      plot.psdr(newlist, d=1)
      return(newlist)
    }else{
      return(newlist)
    }
  }
  else if(sum(as.character(loss) == type.list2) != 0){
    #r.H <- matrix(0, ncol=p+1, nrow=H)
    #weight_list <- seq(0, 1, length=H+1)[2:(H)]
    weight_list <- seq(0, 1, length=H+2)[2:(H+1)]

    if(as.character(loss) == "lssvm"){
      A <- vector(mode = "list", length = length(qprob))
      if(sum(unique(y)) == 0)
        stop("response variable should be continuous!")
      # estimate bases
      r.H <- matrix(0, ncol=p+1, nrow=length(qprob))
      for (s in 1:length(qprob)){
        y.tilde <- rep(1, nrow(x))
        y.tilde[y < qy[s]] <- -1  #s
        A[[s]] <- M <- solve(n * cov.x.star/lambda+t(x.star) %*% x.star)
        C <- t(x.star) %*%  y.tilde
        r.H[s,] <- M %*% C
      }
    }
    if(as.character(loss)=="wlssvm"){
      A <- vector(mode = "list", length = H)
      r.H <- matrix(0, ncol=p+1, nrow=H)
      if(sum(unique(y)) != 0)
        stop("response variable should be a binary type!")
      for (i in 1:H){
        W <- diag(c(ifelse(y==1, 1-weight_list[i], weight_list[i])))
        A[[i]] <- solve(n*cov.x.star/lambda + t(x.star) %*% W %*% x.star)
        C <- t(x.star) %*% W %*% y
        r.H[i,] <- A[[i]] %*% C
      }
    }
    Working_mat <- t(r.H[,1:p])%*%r.H[,1:p]
    eigen.Mn <- eigen(Working_mat)
    newlist <- list("loss" = loss, "x" = x, "y" = y, "M"=Working_mat, "evectors" = eigen.Mn$vectors,"evalues" = eigen.Mn$values, "N"=n, "Xbar"=apply(x, 2, mean), "r"=r.H, "A"=A)

    class(newlist) <- c("psdr", class(newlist))
    #class(newlist) <- "psdr"
    structure(class = "psdr", newlist)

    if(plot == TRUE){
      plot.psdr(newlist, d=1)
      return(newlist)
    }else{
      return(newlist)
    }
  }
}



#' @noRd
#' @export
print.psdr <- function(x, ...) {
  obj <- x
  if(!(obj$loss %in% c("lssvm", "wlssvm"))){
  d <- list(evalues = obj$evalues, evectors = obj$evectors)
  }else{
  d <- list(evalues = obj$evalues, evectors = obj$evectors, r = obj$r, A = obj$A)
  }
  print(d, ...)
  invisible(d)
}













