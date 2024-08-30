#'Pre-embedded arbitrary loss
#' @noRd
#'
fn_arbitrary <- function(object){UseMethod("fn_arbitrary")}
fn_arbitrary_loss <- function(object){UseMethod("fn_arbitrary_loss")}
fn_arbitrary_binary_loss <- function(object){UseMethod("fn_arbitrary_binary_loss")}
fn_arbitrary_nonlinear_loss <- function(object){UseMethod("fn_arbitrary_nonlinear_loss")}
fn_arbitrary_nonlinear_binary_loss <- function(object){UseMethod("fn_arbitrary_nonlinear_binary_loss")}
wvec <- function(object){UseMethod("wvec")}
d2 <- function(object){UseMethod("d2")}
E <- function(object){UseMethod("E")}

E <- function (...) {eval(parse(text=paste(...,collapse=" ")))}

fn_arbitrary <- function(loss){
  E(loss)
}

mylogistic <- function(u) log(1+exp(-u))





############################
### Derivative functions ###
############################
##linear##
fn_arbitrary_loss <- function(x, y, theta, prob=0.5, lambda, loss, mtype){
  #type <- formals(loss)$type
  if (is.null(mtype) == T){     #check lengths
    mtype <- "m"
  }
  if(as.character(mtype) == "r"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- ft(u,prob)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- ft(m)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }
  return(loss.output)
}

fn_arbitrary_binary_loss <- function(x, y, prob, theta, lambda, loss, mtype){
  #type <- formals(loss)$type
  if (is.null(mtype) == T){     #check lengths
    mtype <- "m"
  }
  weight <- (1-prob)*(as.numeric(y==1)) + (prob)*(as.numeric(y==-1))
  if(as.character(mtype)=="r"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(u,prob)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(m)
    loss.output <- lambda*mean(losses) + t(theta)%*%theta
  }
  return(loss.output)
}

##nonlinear##
fn_arbitrary_nonlinear_loss <- function(x, y, theta, prob=NULL,lambda, loss, mtype){
  losses <- c()
  A <- t(x) %*% x
  #type <- formals(loss)$type
  if (is.null(mtype) == T){     #check lengths
    mtype <- "m"
  }
  if(as.character(mtype)=="r"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- ft(u,prob)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- ft(m)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  }
  return(loss.output)
}
fn_arbitrary_nonlinear_binary_loss <- function(x, y, theta, prob, lambda, loss, mtype){
  losses <- c()
  A <- t(x) %*% x
  #type <- formals(loss)$type
  if (is.null(mtype) == T){
    mtype <- "m"
  }
  weight <- (1-prob)*(as.numeric(y==1)) + (prob)*(as.numeric(y==-1))
  if(as.character(mtype)=="r"){
    u <- y - x%*%theta
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(u,prob)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  }else{
    m <- (x %*% theta)*y
    ft <- fn_arbitrary(loss)
    losses <- weight*ft(m)
    loss.output <- lambda*(sum(losses)/nrow(x)) + (t(theta) %*% A %*% theta)/nrow(x)
  }
  return(loss.output)
}


wvec <- function(x, y = y) {
  w <- NULL
  w[y == 1] <- 1-x
  w[y ==-1] <- x
  w
}

d2 <- function(Bhat, B) {
  # This function reports performance in terms of the Frobenius norm of the projection matrix difference

  return(norm(B%*%solve(t(B)%*%B)%*%t(B)-Bhat%*%solve(t(Bhat)%*%Bhat)%*%t(Bhat),"f"))
}






