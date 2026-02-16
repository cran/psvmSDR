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

#' @noRd
# --------------------------
# Input validation utilities
# --------------------------
.check_input <- function(x, y) {
  if (!is.matrix(x) && !is.data.frame(x))
    stop("x must be a numeric matrix or data.frame.")
  if (length(y) != nrow(x))
    stop("length(y) must match nrow(x).")
  invisible(TRUE)
}

# --------------------------
# Response-type detection
# --------------------------
#' @noRd
.response_type <- function(y) {
  y0 <- stats::na.omit(y)
  if (is.factor(y0) || is.ordered(y0)) {
    lev <- levels(y0)
    lev_num <- suppressWarnings(as.numeric(lev))
    is_numeric_levels <- !any(is.na(lev_num))
    if (is_numeric_levels) {
      uniq <- unique(lev_num)
      if (length(uniq) == 2L) return("binary")
      if (length(uniq) <= 5L) return("ordinal")
      return("categorical")
    }
    if (length(lev) == 2L) return("binary")
    return(if (is.ordered(y0)) "ordinal" else "categorical")
  }
  if (is.character(y0)) {
    u <- unique(trimws(tolower(y0)))
    return(if (length(u) == 2L) "binary" else "categorical")
  }
  if (is.logical(y0)) return("binary")
  if (is.numeric(y0) || is.integer(y0)) {
    u <- sort(unique(y0))

    if (length(u) == 2L) return("binary")
    if (length(u) <= 5L && all(u == round(u))) return("ordinal")
    return("continuous")
  }
  return("unknown")
}


# --------------------------
# Safe linear solver (no explicit inverse)
# --------------------------
.safe_solve <- function(A, b, eps = 1e-6, max_tries = 5) {
  p <- nrow(A)
  for (i in 0:max_tries) {
    eps_i <- eps * 10^i
    A_eps <- A + diag(eps_i, p)
    chol_try <- try(chol(A_eps), silent = TRUE)
    if (!inherits(chol_try, "try-error")) {
      return(backsolve(chol_try, forwardsolve(t(chol_try), b)))
    }
  }
  return(qr.solve(A, b))
}



# Internal helper: preprocessing for real-time updates
.rt_prepare <- function(x, Xbar = NULL) {
  n <- nrow(x); p <- ncol(x)
  if (is.null(Xbar)) Xbar <- colMeans(x)
  x.centered <- sweep(x, 2, Xbar, FUN = "-")
  x.star <- cbind(x.centered, -1)                 # (n x (p+1))
  cov.x.star <- stats::cov(x.star)                # (p+1) x (p+1)
  list(x.star = x.star, cov.x.star = cov.x.star, Xbar = Xbar)
}

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






