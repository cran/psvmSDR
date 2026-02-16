#' Plot sufficient predictors from a psdr object
#' @description
#' Produces scatter plots of the sufficient predictors obtained from psdr().
#' For continuous responses, the function plots Y versus each selected
#' sufficient predictor along with an optional lowess curve. For binary
#' responses, a two-dimensional scatter plot of the first two sufficient
#' predictors is produced with class-specific point colors.
#'
#' Additional graphical parameters may be passed to the underlying plot()
#' function. The plot is intended as a diagnostic tool to visualize the
#' estimated central subspace and assess how well the sufficient predictors
#' capture the relationship between X and Y.
#'
#' @param x object from the function \code{psdr()}
#' @param d number of sufficient predictors. Default is 1.
#' @param lowess draw a locally weighted scatterplot smoothing curve. Default is TRUE.
#' @param col color vector for points (optional; defaults depend on response type)
#' @param line.col color for lowess smoothing line (default = "red")
#' @param pch plotting character (default = 16)
#' @param lwd line width for smoothing curve (default = 1.2)
#' @param xlab label for x-axis (default depends on d)
#' @param ylab label for y-axis (default depends on response type)
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @return A scatter plot with sufficient predictors and the lowess curve is overlayed as default.
#' @author Jungmin Shin, \email{c16267@gmail.com}, Seung Jun Shin, \email{sjshin@korea.ac.kr}, Andreas Artemiou \email{artemiou@uol.ac.cy}
#' @seealso \code{\link{psdr_bic}}, \code{\link{psdr}}
#' @examples
#'\donttest{
#' set.seed(1)
#' n <- 200; p <- 5;
#' x <- matrix(rnorm(n*p, 0, 2), n, p)
#' y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
#' obj <- psdr(x, y)
#' plot(obj)
#' }
#' @import stats
#' @importFrom graphics plot points lines par grid
#' @method plot psdr
#' @export

plot.psdr <- function(x, ..., d = 1, lowess = TRUE, col = NULL,
                      line.col = "red", pch = 16, lwd = 1.2,
                      xlab = NULL, ylab = NULL) {

  dots <- list(...)

  # graphical layout handling
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 1))

  if (!inherits(x, "psdr"))
    stop("This method is only for objects of class 'psdr'.")
  if (!is.numeric(d) || length(d) != 1L || d < 1)
    stop("d must be a positive integer.")

  object  <- x
  ytype   <- object$fit$ytype
  evectors <- object$evectors
  x_proj  <- as.matrix(object$x) %*% evectors

  if (ytype == "continuous") {

    y <- object$y
    d <- min(d, ncol(evectors))

    # auto layout
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
          x = x_proj[, j],
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
        sm <- stats::lowess(x_proj[, j], y, f = 0.7)
        graphics::lines(sm, col = line.col, lwd = lwd)
      }
    }

  } else if (ytype == "binary") {

    y <- as.numeric(factor(object$y))
    d <- min(max(2, d), ncol(evectors))
    if (ncol(x_proj) < 2)
      stop("At least two sufficient predictors are required for binary visualization.")

    col2 <- if (!is.null(col)) col else c("steelblue", "tomato2")
    xlab_main <- if (!is.null(xlab)) xlab else "Sufficient predictor 1"
    ylab_main <- if (!is.null(ylab)) ylab else "Sufficient predictor 2"

    plot_args <- c(
      list(
        x = x_proj[, 1],
        y = x_proj[, 2],
        xlab = xlab_main,
        ylab = ylab_main,
        type = "n"
      ),
      dots
    )
    do.call(graphics::plot, plot_args)

    graphics::points(x_proj[y == 1, 1], x_proj[y == 1, 2],
                     col = col2[1], pch = pch)
    graphics::points(x_proj[y == 2, 1], x_proj[y == 2, 2],
                     col = col2[2], pch = pch)
    grid(col = "gray", lty = 1)
    legend("topright", legend = levels(factor(object$y)),
           col = col2, pch = pch, bty = "n")

  } else {
    stop("Plotting is supported only for continuous or binary responses.")
  }

  invisible(object)
}


