#' Plot sufficient predictors from an npsdr object
#'@description
#'
#' Creates diagnostic scatter plots of nonlinear sufficient predictors produced
#' by npsdr(). The function visualizes the estimated transformed directions
#' and optionally overlays a lowess smoothing curve for continuous responses.
#'
#' Additional graphical arguments can be provided. These plots help assess
#' nonlinear structure in the data and evaluate how effectively the kernel SDR
#' method reduces dimensionality.
#'
#' @param x object from \code{npsdr()}
#' @param d number of sufficient predictors to plot (default = 1)
#' @param lowess draw a lowess curve for continuous responses (default = TRUE)
#' @param col point color(s)
#' @param line.col color for lowess smoothing line (default = "red")
#' @param pch point character (default = 16)
#' @param lwd line width for smoothing (default = 1.2)
#' @param xlab x-axis label (default depends on predictor index)
#' @param ylab y-axis label (default = "Y" for continuous)
#' @param ... additional graphical parameters for \code{plot()}
#'
#' @return A scatter plot with sufficient predictors.
#' @author Jungmin Shin, \email{c16267@gmail.com}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#' @seealso \code{\link{npsdr_x}}, \code{\link{npsdr}}
#' @examples
#'\donttest{
#'set.seed(1)
#'n <- 200;
#'p <- 5;
#'x <- matrix(rnorm(n*p, 0, 2), n, p)
#'y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
#'obj_kernel <- npsdr(x, y, plot=FALSE)
#'plot(obj_kernel, d = 1)
#' }
#' @import stats
#' @importFrom graphics plot points lines par grid legend
#' @export


plot.npsdr <- function(x, ..., d = 1, lowess = TRUE,
                       col = NULL, line.col = "red",
                       pch = 16, lwd = 1.2,
                       xlab = NULL, ylab = NULL) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 1))

  if (!inherits(x, "npsdr"))
    stop("This method is only for objects of class 'npsdr'.")
  if (!is.numeric(d) || length(d) != 1L || d < 1)
    stop("d must be a positive integer.")

  object <- x
  y        <- object$y
  ytype    <- object$fit$response.type
  v        <- object$evectors
  phi_mat  <- phix(object$x, object)

  d <- min(d, ncol(v))

  dots <- list(...)

  if (ytype == "continuous") {

    if (d <= 2) {
      par(mfrow = c(1, d))
    } else {
      side <- ceiling(sqrt(d))
      par(mfrow = c(side, side))
    }

    for (j in 1:d) {
      xlab_j <- if (!is.null(xlab)) xlab else paste("Sufficient predictor", j)
      ylab_j <- if (!is.null(ylab)) ylab else "Y"
      col_j  <- if (!is.null(col)) col else "black"

      plot_args <- c(
        list(
          x = phi_mat[, j],
          y = y,
          xlab = xlab_j,
          ylab = ylab_j,
          col  = col_j,
          pch  = pch
        ),
        dots
      )

      do.call(graphics::plot, plot_args)
      grid(col = "gray", lty = 1)

      if (isTRUE(lowess)) {
        sm <- stats::lowess(phi_mat[, j], y, f = 0.7)
        graphics::lines(sm, col = line.col, lwd = lwd)
      }
    }

  } else if (ytype == "binary") {

    y_bin <- as.numeric(factor(y))
    if (ncol(phi_mat) < 2)
      stop("At least two sufficient predictors are needed for binary plotting.")

    col2 <- if (!is.null(col)) col else c("steelblue", "tomato2")
    xlab_main <- if (!is.null(xlab)) xlab else "Sufficient predictor 1"
    ylab_main <- if (!is.null(ylab)) ylab else "Sufficient predictor 2"

    plot_args <- c(
      list(
        x = phi_mat[, 1],
        y = phi_mat[, 2],
        xlab = xlab_main,
        ylab = ylab_main,
        type = "n"
      ),
      dots
    )
    do.call(graphics::plot, plot_args)

    graphics::points(phi_mat[y_bin == 1, 1], phi_mat[y_bin == 1, 2],
                     col = col2[1], pch = pch)
    graphics::points(phi_mat[y_bin == 2, 1], phi_mat[y_bin == 2, 2],
                     col = col2[2], pch = pch)
    grid(col = "gray", lty = 1)

    legend("topright",
           legend = levels(factor(y)),
           col = col2,
           pch = pch, bty = "n")

  } else {
    warning("Plotting supported only for continuous or binary responses.")
  }

  invisible(object)
}
