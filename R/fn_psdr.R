#' Unified linear principal sufficient dimension reduction methods
#' @description
#' This function implements a unified framework for linear principal SDR methods.
#' It provides a single interface that covers many existing principal-machine
#' approaches, such as principal SVM, weighted SVM, logistic, quantile, and
#' asymmetric least squares SDR. The method estimates the central subspace by
#' constructing a working matrix M derived from user-specified loss functions,
#' slicing or weighting schemes, and regularization.
#'
#' The function is designed for both continuous responses and binary
#' classification (with any two-level coding). Users may choose among several
#' built-in loss functions or supply a custom loss function.
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
#' The output includes the estimated eigenvalues and eigenvectors of M,
#' which form the basis of the estimated central subspace, as well as
#' detailed metadata used to summarize model fitting and diagnostics.
#'
#' @param x input matrix, of dimension \code{nobs} x \code{nvars}; each row is an observation vector.
#' @param y response variable, either continuous or binary (any 2-level coding; e.g., -1/1, 0/1, 1/2, TRUE/FALSE, factor/character).
#' @param loss pre-specified loss functions belongs to \code{"svm", "logit", "l2svm", "wsvm", "qr", "asls", "wlogit", "wl2svm", "lssvm", "wlssvm"}, and user-defined loss function object also can be used formed by inside double (or single) quotation mark. Default is 'svm'.
#' @param h unified control for slicing or weighting; accepts either an integer or a numeric vector.
#' @param lambda regularization parameter (default \code{1}).
#' @param eps convergence threshold on parameter change (default \code{1e-5}).
#' @param max.iter maximum number of iterations (default \code{100}).
#' @param eta learning rate for gradient descent (default \code{0.1}).
#' @param mtype a margin type, which is either margin ("m") or residual ("r") (See, Table 1 in the manuscript). Only need when user-defined loss is used. Default is "m".
#' @param plot logical; if TRUE, produces diagnostic plot.
#' @return An object of S3 class \code{"psdr"} containing
#'   \itemize{
#'     \item \code{M}: working matrix
#'     \item \code{evalues}, \code{evectors}: eigen decomposition of \code{M}
#'     \item \code{fit}: metadata (n, p, ytype, hyperparameters, per-slice iteration/convergence info)
#'   }
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
#' @seealso \code{\link{psdr_bic}}, \code{\link{rtpsdr}}
#'@examples
#'\donttest{
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
#'## --------------------------
#'## User defined cutoff points
#'## --------------------------
#'obj_cut <- psdr(x, y, h = c(0.1, 0.3, 0.5, 0.7))
#'print(obj_cut)
#'
#'## --------------------------------
#'## Linear PM (Binary classification)
#'## --------------------------------
#' obj_wsvm <- psdr(x, y.tilde, loss="wsvm")
#' plot(obj_wsvm)
#'
#'## ----------------------------
#'## User-defined loss function
#'## ----------------------------
#' mylogistic <- function(u) log(1+exp(-u))
#' psdr(x, y, loss="mylogistic")
#'
#'## ----------------------------
#'## Real-data example: iris (binary subset)
#'## ----------------------------
#'iris_binary <- droplevels(subset(iris, Species %in% c("setosa", "versicolor")))
#'psdr(x = iris_binary[, 1:4], y = iris_binary$Species, plot = TRUE)
#'}
#'
#'@import stats graphics
#'@importFrom utils head
#'@export psdr
#'
#'
psdr <- function(x, y, loss = "svm", h = 10, lambda = 1, eps = 1e-5, max.iter = 100, eta = 0.1, mtype = "m", plot = FALSE) {

  # ---- Input validation ---- #
  .check_input(x, y)
  ytype <- .response_type(y)

  if (!ytype %in% c("binary", "continuous"))
    stop("Unsupported response type: only continuous or binary responses are allowed. Detected: ", ytype)

  if (missing(loss) || loss == "svm") {
    if (ytype == "binary") {
      loss <- "wsvm"
    } else {
      loss <- "svm"
    }
  }

  n <- nrow(x)
  p <- ncol(x)

  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0) stop("lambda must be a positive scalar.")
  if (!is.numeric(eta)     || length(eta)     != 1L || eta <= 0) stop("eta must be a positive scalar.")
  if (!is.numeric(max.iter)|| length(max.iter)!= 1L || max.iter<= 0) stop("max.iter must be a positive integer.")
  if (!is.numeric(eps)     || length(eps)     != 1L || eps <= 0) stop("eps must be a positive scalar.")


  if (ytype == "binary") {
    # --- Binary response: interpret h as class-weight grid ---
    if (is.list(h)) {
      weight_cutpoints <- if (!is.null(h$weight)) h$weight else NULL
    } else if (is.numeric(h) && length(h) > 1) {
      # direct numeric vector input (e.g., h = c(0.1, 0.3, 0.6, 0.8))
      if (any(h < 0 | h > 1))
        stop("When y is binary, h must contain numeric values in [0, 1].")
      weight_cutpoints <- sort(unique(h))
    } else if (length(h) == 1L && is.numeric(h)) {
      # integer number of class-weight grid points
      weight_cutpoints <- seq(0, 1, length = as.integer(h) + 2)[2:(as.integer(h) + 1)]
    } else {
      stop("Invalid input for h.")
    }

    weight_list <- weight_cutpoints
    H <- length(weight_list)
    pi.grid <- weight_list  # for consistent downstream use

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
  }

  if (ytype == "continuous") {
    qprob <- (1:(H - 1)) / H
    qy <- stats::quantile(y, qprob)
  } else {
    qy <- NULL  # not needed for binary case
  }

  bar.x <- colMeans(x)
  x.star <- cbind(scale(x, center = TRUE, scale = FALSE), -1)  # (n x (p+1))
  cov.x.star <- stats::cov(x.star)
  cov.x <- stats::cov(x)

  eigx <- eigen(cov.x)
  D <- diag(sqrt(eigx$values), nrow = p, ncol = p)
  V <- eigx$vectors
  inv.sd.x <- diag(1 / sqrt(eigx$values), nrow = p, ncol = p) %*% t(V)

  centered.x <- t(x) - bar.x     # p x n
  z.new <- t(inv.sd.x %*% centered.x)  # n x p

  .to_pm_binary <- function(y) {
    y0 <- stats::na.omit(y)
    if (is.logical(y0)) return(ifelse(y, 1, -1))
    if (is.factor(y0) || is.character(y0)) {
      lev <- if (is.factor(y0)) levels(y0) else unique(as.character(y0))
      lev <- sort(lev)
      return(ifelse(as.character(y) == lev[2], 1, -1))
    }
    u <- sort(unique(y0))
    if (length(u) != 2L) stop("Binary coercion requires exactly 2 unique values.")
    return(ifelse(y == u[2], 1, -1))
  }
  y.bin <- if (ytype == "binary") .to_pm_binary(y) else NULL
  type.list  <- c("svm","logit","l2svm","wsvm","qr","asls","wlogit","wl2svm")
  type.list2 <- c("lssvm","wlssvm")

  init <- rep(1, p)
  eigen.mat <- diag(1, p, p)
  theta.new <- rep(0, p)

  slice_iters <- integer(length(pi.grid))
  slice_converged <- logical(length(pi.grid))
  slice_obj <- rep(NA_real_, length(pi.grid))


  if (as.character(loss) %in% type.list2) {

    if (loss == "lssvm"  && ytype == "binary")
      stop("lssvm requires a continuous response.")
    if (loss == "wlssvm" && ytype == "continuous")
      stop("wlssvm requires a binary response.")

    A <- vector("list", length = H)
    r.H <- matrix(0, ncol = p + 1, nrow = H)

    if (loss == "lssvm") {
      for (s in 1:H) {
        thr_idx <- min(s, length(qy))
        y.tilde <- ifelse(y < qy[thr_idx], -1, 1)
        A[[s]] <- .safe_solve(n * cov.x.star / lambda + t(x.star) %*% x.star, diag(p + 1))
        C <- t(x.star) %*% y.tilde
        r.H[s, ] <- A[[s]] %*% C
      }
    } else if (loss == "wlssvm") {
      for (i in 1:H) {
        wgt <- weight_list[i]
        W <- diag(ifelse(y.bin == 1, 1 - wgt, wgt), nrow = n, ncol = n)
        A[[i]] <- .safe_solve(n * cov.x.star / lambda + t(x.star) %*% W %*% x.star, diag(p + 1))
        C <- t(x.star) %*% W %*% y.bin
        r.H[i, ] <- A[[i]] %*% C
      }
    }

    Working_mat <- t(r.H[, 1:p]) %*% r.H[, 1:p]
    eigen.Mn <- eigen(Working_mat)

    meta <- list(
      n = n, p = p, ytype = ytype,
      lambda = lambda, eta = eta, eps = eps, max.iter = max.iter,
      converged = TRUE,
      n_iter = NA_integer_,
      per_slice = data.frame(slice = seq_along(pi.grid),
                             iter = NA_integer_, converged = TRUE, obj = NA_real_)
    )

    out <- list(
      loss = loss, x = x, y = y,
      M = Working_mat,
      evalues = eigen.Mn$values, evectors = eigen.Mn$vectors,
      N = n, Xbar = bar.x, r = r.H, A = A,
      fit = meta
    )
    class(out) <- c("psdr", class(out))
    if (plot) plot.psdr(out, d = 1)
    return(invisible(out))
  }

  if (as.character(loss) %in% type.list) {

    qlen <- length(pi.grid)
    w.init  <- matrix(init, nrow = p, ncol = qlen)
    w.final <- matrix(0,    nrow = p, ncol = qlen)

    .finalize_slice <- function(s, w_s, iter, obj_val, converged_flag) {
      w.final[, s] <<- w_s
      slice_iters[s] <<- iter
      slice_converged[s] <<- converged_flag
      slice_obj[s] <<- obj_val
    }

    ## ========== svm ==========
    if (as.character(loss) == "svm") {
      if (ytype == "binary") stop("svm requires a continuous response.")
      for (s in 1:qlen) {
        y.tilde.new <- ifelse(y < qy[s], -1, 1)
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.tilde <- y.tilde.new; nloc <- nrow(z)
          for (k in 1:p) {
            margin.v <- (z %*% w[, s]) * y.tilde
            deriv <- -z[, k] * y.tilde * as.numeric((1 - margin.v) > 0)
            derivative.j <- lambda * mean(deriv) + 2 * w[k, s]
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== l2svm ==========
    if (as.character(loss) == "l2svm") {
      if (ytype == "binary") stop("l2svm requires a continuous response.")
      for (s in 1:qlen) {
        y.tilde.new <- ifelse(y < qy[s], -1, 1)
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.tilde <- y.tilde.new; nloc <- nrow(z)
          for (k in 1:p) {
            margin.v <- (z %*% w[, s]) * y.tilde
            deriv <- -z[, k] * y.tilde * as.numeric((1 - margin.v) > 0) * 2 * (1 - margin.v)
            derivative.j <- lambda * mean(deriv) + 2 * w[k, s]
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== logit ==========
    if (as.character(loss) == "logit") {
      if (ytype == "binary") stop("logit requires a continuous response.")
      for (s in 1:qlen) {
        y.tilde.new <- ifelse(y < qy[s], -1, 1)
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.tilde <- y.tilde.new
          for (k in 1:p) {
            margin.v <- (z %*% w[, s]) * y.tilde
            deriv <- -z[, k] * y.tilde * (1 / (1 + exp(margin.v)))
            derivative.j <- lambda * mean(deriv) + 2 * w[k, s]
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== wsvm (binary) ==========
    if (as.character(loss) == "wsvm") {
      if (ytype == "continuous") stop("wsvm requires a binary response.")
      for (s in 1:qlen) {
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        wgt <- weight_list[s %% length(weight_list) + (length(weight_list) * (s %% length(weight_list) == 0))]
        for (iter in 1:max.iter) {
          z <- z.new; y.bi <- y.bin
          for (k in 1:p) {
            margin.v <- (z %*% w[, s]) * y.bi
            weight <- (1 - pi.grid[s]) * as.numeric(y.bi == 1) + pi.grid[s] * as.numeric(y.bi == -1)
            deriv <- -weight * z[, k] * y.bi * as.numeric((1 - margin.v) > 0)
            derivative.j <- lambda * mean(deriv) + 2 * w[k, s]
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== wlogit (binary) ==========
    if (as.character(loss) == "wlogit") {
      if (ytype == "continuous") stop("wlogit requires a binary response.")
      for (s in 1:qlen) {
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.bi <- y.bin
          for (k in 1:p) {
            margin.v <- (z %*% w[, s]) * y.bi
            weight <- (1 - pi.grid[s]) * as.numeric(y.bi == 1) + pi.grid[s] * as.numeric(y.bi == -1)
            deriv <- weight * (-z[, k]) * y.bi * (1 / (1 + exp(margin.v)))
            derivative.j <- lambda * mean(deriv) + 2 * w[k, s]
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== wl2svm (binary) ==========
    if (as.character(loss) == "wl2svm") {
      if (ytype == "continuous") stop("wl2svm requires a binary response.")
      for (s in 1:qlen) {
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.bi <- y.bin
          for (k in 1:p) {
            margin.v <- (z %*% w[, s]) * y.bi
            weight <- (1 - pi.grid[s]) * as.numeric(y.bi == 1) + pi.grid[s] * as.numeric(y.bi == -1)
            deriv <- -z[, k] * y.bi * as.numeric((1 - margin.v) > 0) * 2 * (1 - margin.v) * weight
            derivative.j <- lambda * mean(deriv) + 2 * w[k, s]
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== qr (quantile regression) ==========
    if (as.character(loss) == "qr") {
      for (s in 1:length(pi.grid)) {
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.new <- y
          for (k in 1:p) {
            u <- y.new - (z %*% w[, s])
            derivative.j <- 2 * w[k, s] +
              lambda * (1 / length(y)) * sum(-z[, k] * (pi.grid[s] * as.numeric(u > 0) + (1 - pi.grid[s]) * as.numeric(u <= 0)))
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## ========== asls ==========
    if (as.character(loss) == "asls") {
      for (s in 1:length(pi.grid)) {
        converged_flag <- FALSE
        obj_val <- NA_real_
        w <- w.init
        for (iter in 1:max.iter) {
          z <- z.new; y.new <- y
          for (k in 1:p) {
            u <- y.new - (z %*% w[, s])
            derivative.j <- 2 * w[k, s] +
              lambda * (1 / length(y)) *
              sum((-z[, k] * 2 * u) * (pi.grid[s] * as.numeric(u >= 0) + (1 - pi.grid[s]) * as.numeric(u <= 0)))
            theta.new[k] <- w[k, s] - eta * derivative.j
          }
          old <- w[, s]
          w[, s] <- theta.new
          obj_val <- mean(abs(theta.new - old))
          if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
        }
        .finalize_slice(s, w[, s], iter, obj_val, converged_flag)
        w.init[, s] <- w[, s]
      }
    }

    ## Build working matrix and return
    psi <- t(inv.sd.x) %*% w.final
    Mn <- matrix(0, p, p)
    for (hidx in 1:length(pi.grid)) {
      Mn <- Mn + psi[, hidx, drop = FALSE] %*% t(psi[, hidx, drop = FALSE])
    }
    eigen.Mn <- eigen(Mn)

    meta <- list(
      n = n, p = p, ytype = ytype,
      lambda = lambda, eta = eta, eps = eps, max.iter = max.iter,
      converged = all(slice_converged),
      n_iter = sum(slice_iters),
      per_slice = data.frame(slice = seq_along(pi.grid),
                             iter = slice_iters,
                             converged = slice_converged,
                             obj = slice_obj)
    )

    out <- list(
      loss = loss, x = x, y = y,
      M = Mn,
      evalues = eigen.Mn$values, evectors = eigen.Mn$vectors,
      fit = meta
    )
    class(out) <- c("psdr", class(out))
    if (plot) plot.psdr(out, d = 1)
    return(invisible(out))
  }

  ft <- E(loss)
  qlen <- length(pi.grid)
  w.init  <- matrix(init, nrow = p, ncol = qlen)
  w.final <- matrix(0,    nrow = p, ncol = qlen)

  # continuous case
  if (ytype == "continuous") {
    for (s in 1:qlen) {
      y.tilde.new <- ifelse(y < qy[s], -1, 1)
      converged_flag <- FALSE
      obj_val <- NA_real_
      w <- w.init
      interval <- 1.0e-5
      for (iter in 1:max.iter) {
        z <- z.new; y.tilde <- y.tilde.new
        derivative.vec <- rep(0, p)
        for (k in 1:p) {
          star <- (fn_arbitrary_loss(z, y.tilde, theta = w[, s] + interval * eigen.mat[, k],
                                     lambda = lambda, loss = loss, prob = pi.grid[s], mtype = "m") -
                     fn_arbitrary_loss(z, y.tilde, theta = w[, s] - interval * eigen.mat[, k],
                                       lambda = lambda, loss = loss, prob = pi.grid[s], mtype = "m"))
          derivative.vec[k] <- sign(star) * exp(log(abs(star)) - log(2 * interval))
          theta.new[k] <- w[k, s] - eta * derivative.vec[k]
        }
        old <- w[, s]
        w[, s] <- theta.new
        obj_val <- mean(abs(theta.new - old))
        if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
      }
      w.final[, s] <- w[, s]
      slice_iters[s] <- iter; slice_converged[s] <- converged_flag; slice_obj[s] <- obj_val
      w.init[, s] <- w[, s]
    }
  }

  # binary case
  if (ytype == "binary") {
    y.new <- y.bin
    for (s in 1:qlen) {
      converged_flag <- FALSE
      obj_val <- NA_real_
      w <- w.init
      interval <- 1.0e-5
      for (iter in 1:max.iter) {
        z <- z.new; y.bi <- y.new
        derivative.vec <- rep(0, p)
        for (k in 1:p) {
          star <- (fn_arbitrary_binary_loss(z, y.bi, prob = pi.grid[s],
                                            theta = w[, s] + interval * eigen.mat[, k],
                                            lambda = lambda, loss = loss, mtype = "m") -
                     fn_arbitrary_binary_loss(z, y.bi, prob = pi.grid[s],
                                              theta = w[, s] - interval * eigen.mat[, k],
                                              lambda = lambda, loss = loss, mtype = "m"))
          derivative.vec[k] <- sign(star) * exp(log(abs(star)) - log(2 * interval))
          theta.new[k] <- w[k, s] - eta * derivative.vec[k]
        }
        old <- w[, s]
        w[, s] <- theta.new
        obj_val <- mean(abs(theta.new - old))
        if (max(abs(theta.new - old)) < eps) { converged_flag <- TRUE; break }
      }
      w.final[, s] <- w[, s]
      slice_iters[s] <- iter; slice_converged[s] <- converged_flag; slice_obj[s] <- obj_val
      w.init[, s] <- w[, s]
    }
  }

  psi <- t(inv.sd.x) %*% w.final
  Mn <- matrix(0, p, p)

  for (hidx in 1:length(pi.grid)) {
    Mn <- Mn + psi[, hidx, drop = FALSE] %*% t(psi[, hidx, drop = FALSE])
  }

  eigen.Mn <- eigen(Mn)

  meta <- list(
    n = n, p = p, ytype = ytype,
    lambda = lambda, eta = eta, eps = eps, max.iter = max.iter,
    converged = all(slice_converged),
    n_iter = sum(slice_iters),
    per_slice = data.frame(slice = seq_along(pi.grid),
                           iter = slice_iters,
                           converged = slice_converged,
                           obj = slice_obj)
  )

  out <- list(
    loss = loss, x = x, y = y,
    M = Mn,
    evalues = eigen.Mn$values, evectors = eigen.Mn$vectors,
    fit = meta
  )
  class(out) <- c("psdr", class(out))
  if (plot) plot.psdr(out, d = 1)
  invisible(out)
}


#' @noRd
#' @export
print.psdr <- function(x, digits = 4, ...) {
  cat("\n--- Principal Sufficient Dimension Reduction (linear) ---\n")
  cat("Loss:", x$loss,
      " | n:", x$fit$n, " p:", x$fit$p,
      " | Response:", x$fit$ytype, "\n")
  cat("Lambda:", x$fit$lambda,
      " | Eta:", x$fit$eta,
      " | Max.iter:", x$fit$max.iter,
      " | Converged:", isTRUE(x$fit$converged), "\n")

  cat("\nEigenvalues (first 5): ",
      paste(round(head(x$evalues, 5), digits), collapse = ", "), "\n")

  cat("\nEigenvectors (columns are SDR directions):\n")
  print(round(x$evectors, digits))
  invisible(x)
}


#' @noRd
#' @export
summary.psdr <- function(object, digits = 4, ...) {
  cat("\n=== Summary of psdr Object ===\n")

  # --- Basic info ---
  cat("Loss function:", object$loss, "\n")
  cat("Sample size (n):", object$fit$n,
      " | Variables (p):", object$fit$p,
      " | Response type:", object$fit$ytype, "\n")
  cat("Lambda:", object$fit$lambda,
      " | Eta:", object$fit$eta,
      " | Eps:", object$fit$eps,
      " | Max.iter:", object$fit$max.iter, "\n")

  # --- Eigen decomposition ---
  cat("\n--- Eigen Decomposition of Working Matrix (M) ---\n")
  cat("Top eigenvalues (up to 10):\n")
  print(round(head(object$evalues, 10), digits))

  cat("\nEstimated Eigenvectors (columns = central subspace basis):\n")
  print(round(object$evectors, digits))

  # --- Convergence diagnostics ---
  if (!is.null(object$fit$per_slice)) {
    cat("\n--- Per-slice diagnostics ---\n")
    print(head(object$fit$per_slice, 10))
    if (nrow(object$fit$per_slice) > 10) cat("...\n")
  }

  cat("\nConvergence summary:\n")
  cat("  Total iterations:", object$fit$n_iter,
      " | All slices converged:", isTRUE(object$fit$converged), "\n")

  invisible(object)
}
