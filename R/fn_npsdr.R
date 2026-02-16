#'A unified Principal sufficient dimension reduction method via kernel trick
#'@description
#' This function extends principal SDR to nonlinear relationships between
#' predictors and the response using a kernel feature map. The kernel basis
#' is constructed internally using a data-driven number of basis functions,
#' and the working matrix is formed analogously to linear principal SDR but
#' in the transformed feature space.
#'
#' Users may choose from built-in loss functions or provide a custom loss
#' through the same interface as psdr(). The method supports both continuous
#' and binary responses and can visualize the nonlinear sufficient predictors.
#'
#' The output contains the kernel basis object, the working matrix M,
#' eigenvalues and eigenvectors, and detailed fitting metadata.
#'@param x data matrix
#'@param y either continuous or (+1,-1) typed binary response vector
#'@param loss pre-specified loss functions belongs to \code{"svm", "logit", "l2svm", "wsvm", "qr", "asls", "wlogit", "wl2svm", "lssvm", "wlssvm"}, and user-defined loss function object also can be used formed by inside double (or single) quotation mark. Default is 'svm'.
#'@param h unified control for slicing or weighting; accepts either an integer or a numeric vector.
#'@param lambda hyperparameter for the loss function. default value is 1
#'@param b number of basis functions for a kernel trick, floor(length(y)/3) is default
#'@param eps threshold for stopping iteration with respect to the magnitude of derivative, default value is 1.0e-4
#'@param max.iter maximum iteration number for the optimization process. default value is 30
#'@param eta learning rate for gradient descent method. default value is 0.1
#' @param mtype type of margin, either "m" or "r" refer margin and residual, respectively (See, Table 1 in the pacakge manuscript). When one use user-defined loss function this argument should be specified. Default is "m".
#' @param plot If \code{TRUE} then it produces scatter plots of \eqn{Y} versus the first sufficient predictor. The default is FALSE.
#' @return An object of class "npsdr" containing:
#' \itemize{
#'   \item \code{x}, \code{y}: input data
#'   \item \code{M}: working matrix
#'   \item \code{evalues}, \code{evectors}: eigen-decomposition of M
#'   \item \code{obj.psi}: kernel basis object from \code{get.psi()}
#'   \item \code{fit}: metadata (loss, h, lambda, eps, max.iter, eta, b, response.type, cutpoints, weight_cutpoints)
#' }
#' @author Jungmin Shin, \email{c16267@gmail.com}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
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
#'summary(obj_kernel)
#'plot(obj_kernel)
#'
#'}
#'@import stats graphics
#'@importFrom graphics abline grid par
#'@importFrom stats cov rnorm
#'@importFrom utils head
#'@export npsdr

npsdr <- function(x, y, loss = "svm", h = 10, lambda = 1,
                  b = floor(length(y) / 3), eps = 1e-5,
                  max.iter = 100, eta = 0.1, mtype = "m", plot = TRUE) {


  .check_input(x, y)
  if (!is.numeric(lambda) || lambda <= 0) stop("lambda must be positive.")
  if (!is.numeric(eta) || eta <= 0) stop("eta must be positive.")
  if (!is.numeric(max.iter) || max.iter <= 0) stop("max.iter must be positive.")
  if (!is.numeric(eps) || eps <= 0) stop("eps must be positive.")
  if (!is.numeric(b) || b < 2) stop("b must be at least 2.")

  ytype <- .response_type(y)
  if (!ytype %in% c("continuous", "binary"))
    stop("Unsupported response type for npsdr(): ", ytype)

  if (ytype == "binary") {
    # Ensure factor → numeric conversion
    if (is.factor(y)) {
      lev <- levels(y)
      y <- ifelse(y == lev[2], 1, -1)
    } else if (is.character(y)) {
      u <- sort(unique(trimws(tolower(y))))
      y <- ifelse(trimws(tolower(y)) == u[2], 1, -1)
    } else if (is.logical(y)) {
      y <- ifelse(y, 1, -1)
    } else if (is.numeric(y) || is.integer(y)) {
      # numeric binary: map min→-1, max→+1
      u <- sort(unique(y))
      y <- ifelse(y == u[2], 1, -1)
    } else {
      stop("Binary response type could not be safely normalized.")
    }
  }

  if (missing(loss) || loss == "svm") {
    if (ytype == "binary") {
      loss <- "wsvm"
    } else {
      loss <- "svm"
    }
  }


  psi.gen <- get.psi(x, y, b)
  Psi.new <- psi.gen$w
  n <- nrow(Psi.new)
  p <- ncol(Psi.new)


  if (ytype == "binary") {
    if (is.list(h)) {
      weight_cutpoints <- if (!is.null(h$weight)) h$weight else NULL
    } else if (is.numeric(h) && length(h) > 1) {
      if (any(h <= 0 | h >= 1)) {
        warning("Values 0 or 1 in h were replaced by 0.01 and 0.99 for stability.")
        h[h <= 0] <- 0.01
        h[h >= 1] <- 0.99
      }
      weight_cutpoints <- sort(unique(h))
    } else if (length(h) == 1L && is.numeric(h)) {
      weight_cutpoints <- seq(0, 1, length = as.integer(h) + 2)[2:(as.integer(h) + 1)]
    } else {
      stop("Invalid input for h.")
    }
    weight_list <- weight_cutpoints
    pi.grid <- weight_list
    H <- length(weight_list)
    qy <- NULL

  } else {
    if (is.list(h)) {
      cutpoints <- if (!is.null(h$slice)) h$slice else NULL
    } else if (is.numeric(h) && length(h) > 1) {
      if (any(h <= 0 | h >= 1)) {
        h[h <= 0] <- 0.01
        h[h >= 1] <- 0.99
      }
      cutpoints <- sort(unique(h))
    } else if (length(h) == 1L && is.numeric(h)) {
      step <- 1 / as.integer(h)
      cutpoints <- seq(step, 1 - step, by = step)
    } else {
      stop("Invalid input for h.")
    }
    pi.grid <- cutpoints
    H <- length(pi.grid) + 1L
    weight_list <- seq(0, 1, length = H + 2)[2:(H + 1)]
    qprob <- (1:(H - 1)) / H
    qy <- stats::quantile(y, qprob)
  }

  set.seed(1)
  init.theta <- rnorm(sd = 1, n = p)
  theta.new <- rep(0, p)
  qprob <- if (ytype == "continuous") (1:(H - 1)) / H else pi.grid
  w.init  <- matrix(init.theta, nrow = p, ncol = length(qprob))
  w.final <- matrix(0, nrow = p, ncol = length(qprob))
  eigen.mat <- diag(1, p, p)

  type.list  <- c("svm","logit","l2svm","wsvm","qr","asls","wlogit","wl2svm")
  type.list2 <- c("lssvm","wlssvm")


  if (as.character(loss) %in% type.list) {

    # === Continuous: SVM, L2SVM, LOGIT, etc. ===
    if (ytype == "continuous") {
      for (s in 1:length(pi.grid)) {
        y.tilde.new <- ifelse(y < stats::quantile(y, pi.grid[s]), -1, 1)
        w <- w.init
        for (iter in 1:max.iter) {
          A <- t(Psi.new) %*% Psi.new
          for (k in 1:p) {
            margin.v <- as.numeric(Psi.new %*% w[, s]) * y.tilde.new
            deriv <- -Psi.new[, k] * y.tilde.new * as.numeric((1 - margin.v) > 0)
            derivative.j <- lambda * mean(deriv) + 2 * (1 / nrow(Psi.new)) * (A[k, ] %*% w[, s])
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          if (max(abs(theta.new - w[, s])) < eps) break
          w[, s] <- theta.new
        }
        w.final[, s] <- w[, s]
      }
    }

    # === Binary: wsvm, wl2svm, wlogit ===
    if (ytype == "binary") {
      y.bin <- ifelse(y == max(y), 1, -1)
      for (s in 1:length(pi.grid)) {
        w <- w.init
        for (iter in 1:max.iter) {
          A <- t(Psi.new) %*% Psi.new
          for (k in 1:p) {
            margin.v <- as.numeric(Psi.new %*% w[, s]) * y.bin
            weight <- (1 - pi.grid[s]) * as.numeric(y.bin == 1) + pi.grid[s] * as.numeric(y.bin == -1)
            deriv <- -weight * Psi.new[, k] * y.bin * as.numeric((1 - margin.v) > 0)
            derivative.j <- lambda * mean(deriv) + 2 * (1 / nrow(Psi.new)) * (A[k, ] %*% w[, s])
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          if (max(abs(theta.new - w[, s])) < eps) break
          w[, s] <- theta.new
        }
        w.final[, s] <- w[, s]
      }
    }

    # Working matrix
    Mn <- matrix(0, p, p)
    for (hidx in 1:ncol(w.final))
      Mn <- Mn + w.final[, hidx, drop = FALSE] %*% t(w.final[, hidx, drop = FALSE])
    result <- eigen(Mn)
  }


  else if (as.character(loss) %in% type.list2) {
    bar.x <- colMeans(Psi.new)
    x.star <- cbind(scale(Psi.new, center = TRUE, scale = FALSE), -1)
    cov.x.star <- cov(x.star)

    if (as.character(loss) == "lssvm") {
      if (ytype != "continuous") stop("response variable should be continuous!")
      r.H <- matrix(0, ncol = p + 1, nrow = H)
      for (s in 1:(H - 1)) {
        y.tilde <- ifelse(y < stats::quantile(y, s / H), -1, 1)
        M <- .safe_solve(n * cov.x.star / lambda + t(x.star) %*% x.star,
                         crossprod(x.star, y.tilde))
        r.H[s, ] <- M
      }
    }

    if (as.character(loss) == "wlssvm") {
      if (ytype != "binary") stop("response variable should be binary!")
      y.bin <- ifelse(y == max(y), 1, -1)
      r.H <- matrix(0, ncol = p + 1, nrow = H)
      for (pi in weight_list) {
        W <- diag(ifelse(y.bin == 1, 1 - pi, pi), nrow = n)
        M <- .safe_solve(n * cov.x.star / lambda + t(x.star) %*% W %*% x.star,
                         crossprod(x.star, W %*% y.bin))
        r.H[which(weight_list == pi), ] <- M
      }
    }

    Mn <- t(r.H[, 1:p]) %*% r.H[, 1:p]
    result <- eigen(Mn)
  }

  else {
    stop("Unsupported or user-defined loss not implemented.")
  }

  out <- list(
    x = x,
    y = y,
    M = Mn,
    evalues = result$values,
    evectors = result$vectors,
    obj.psi = psi.gen,
    fit = list(
      loss = loss,
      lambda = lambda,
      eta = eta,
      eps = eps,
      max.iter = max.iter,
      h = h,
      b = b,
      mode = "kernel",
      response.type = ytype,
      cutpoints = if (ytype == "continuous") pi.grid else NULL,
      weight_cutpoints = if (ytype == "binary") pi.grid else NULL
    )
  )
  class(out) <- c("npsdr", "psdr")

  if (plot) {
    plot.npsdr(out, d = 1, lowess = FALSE)
    invisible(out)
  } else {
    invisible(out)
  }
}


# ------------------------------------------------------------------------------
# Print and Summary for npsdr (consistent with psdr / rtpsdr)
# ------------------------------------------------------------------------------

#' @noRd
#' @export
print.npsdr <- function(x, ...) {
  cat("--- Nonlinear Principal SDR (kernel) ---\n")
  cat("Loss:", x$fit$loss,
      " | n:", nrow(x$x),
      " | p:", ncol(x$x),
      " | Basis (b):", x$fit$b,
      " | Mode:", x$fit$mode, "\n")
  cat("Lambda:", x$fit$lambda,
      " | Eta:", x$fit$eta,
      " | Max.iter:", x$fit$max.iter, "\n")

  if (!is.null(x$fit$cutpoints)) {
    cat("Cutpoints:", paste(round(x$fit$cutpoints, 3), collapse = ", "), "\n")
  }
  if (!is.null(x$fit$weight_cutpoints) && x$fit$loss %in% c("wsvm", "wlssvm", "wlogit", "wl2svm")) {
    cat("Weight cutpoints:", paste(round(x$fit$weight_cutpoints, 3), collapse = ", "), "\n")
  }

  if (!is.null(x$evalues)) {
    topk <- head(round(x$evalues, 4), 5)
    cat("Eigenvalues (first 5):", paste(topk, collapse = ", "), "\n\n")
  }

  if (!is.null(x$evectors)) {
    cat("Eigenvectors (columns are SDR directions):\n")
    print(round(x$evectors, 4))
  }
  invisible(x)
}



#' @export
summary.npsdr <- function(object, digits = 4, ...) {
  cat("=== Summary of npsdr Object ===\n")

  n <- nrow(object$x)
  p <- ncol(object$x)
  fit <- object$fit

  cat("Loss:", fit$loss,
      " | n:", n,
      " | p:", p,
      " | Basis (b):", fit$b,
      " | Mode:", fit$mode, "\n")

  cat("Lambda:", fit$lambda,
      " | Eta:", fit$eta,
      " | Eps:", fit$eps,
      " | Max.iter:", fit$max.iter, "\n")

  cat("Response type:", fit$response.type, "\n")

  if (!is.null(fit$cutpoints)) {
    cat("Cutpoints:", paste(round(fit$cutpoints, 3), collapse = ", "), "\n")
  }
  if (!is.null(fit$weight_cutpoints) && fit$loss %in% c("wsvm", "wlssvm", "wlogit", "wl2svm")) {
    cat("Weight cutpoints:", paste(round(fit$weight_cutpoints, 3), collapse = ", "), "\n")
  }

  if (!is.null(object$obj.psi)) {
    psi <- object$obj.psi
    cat("\n--- Kernel Feature Information ---\n")
    cat("Bandwidth (bw):", round(psi$bw, digits),
        " | Basis functions (b):", psi$b, "\n")
    cat("Mean of scaled predictors:", round(mean(psi$scaled.x), digits),
        " | SD:", round(sd(psi$scaled.x), digits), "\n")
  }

  cat("\n--- Eigen Decomposition of Working Matrix (M) ---\n")
  if (!is.null(object$evalues)) {
    cat("Top eigenvalues (up to 10):\n")
    print(round(head(object$evalues, 10), digits))
  } else {
    cat("Eigenvalues not available.\n")
  }

  if (!is.null(object$evectors)) {
    cat("\nEstimated Eigenvectors (columns = central subspace basis):\n")
    print(round(object$evectors, digits))
  }

  invisible(object)
}


#' @noRd
get.psi <- function(x, y, b=floor(length(y)/3)) {
  n <- nrow(x)
  x <- scale(x)
  bw <- 1/mean(as.numeric(stats::dist(x)))^2
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
  original_scaled_x <- object$obj.psi$scaled.x
  original_bw <- object$obj.psi$bw

  v <- object$evectors
  w <- object$obj.psi$w
  l <- object$obj.psi$l
  d <- ncol(v)
  p <- ncol(original_scaled_x)

  m <- attr(original_scaled_x, "scaled:center")
  s <- attr(original_scaled_x, "scaled:scale")
  scaled_value <- t((t(value) - m)/s)

  if (is.vector(scaled_value)) {
    temp <- psi.function(scaled_value, original_scaled_x, v[,1:d, drop = F], w, l, original_bw)
  } else if (is.matrix(scaled_value) || is.data.frame(scaled_value)) {
    temp <- t(apply(scaled_value, 1, psi.function, original_scaled_x, v[,1:d, drop = F], w, l, original_bw))
  } else {
    stop("check `str(value)`")
  }
  temp
}



#'@noRd
psi.function <- function(value_scaled, original_scaled_x, v, w, l, kernel_param){

  value_scaled <- matrix(value_scaled, 1, length(value_scaled))

  temp_kernel <- kernel.function(value_scaled, original_scaled_x, kernel_param)
  temp_kernel_centered <- c(temp_kernel - mean(temp_kernel))
  psi_value <- colSums(w * temp_kernel_centered) / l

  rslt <- psi_value %*% v
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



