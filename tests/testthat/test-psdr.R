## ----------------------------
## Linear PM
## ----------------------------
set.seed(1)
n <- 200; p <- 5;
x <- matrix(rnorm(n*p, 0, 2), n, p)
y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
y.tilde <- sign(y)
obj <- psdr(x, y)
print(obj)
plot(obj, d=2)

## --------------------------
## User defined cutoff points
## --------------------------
obj_cut <- psdr(x, y, h = c(0.1, 0.3, 0.5, 0.7))
print(obj_cut)

## --------------------------------
## Linear PM (Binary classification)
## --------------------------------
obj_wsvm <- psdr(x, y.tilde, loss="wsvm")
plot(obj_wsvm)

## ----------------------------
## User-defined loss function
## ----------------------------
mylogistic <- function(u) log(1+exp(-u))
psdr(x, y, loss="mylogistic")

## ----------------------------
## Real-data example: iris (binary subset)
## ----------------------------
iris_binary <- droplevels(subset(iris, Species %in% c("setosa", "versicolor")))
psdr(x = iris_binary[, 1:4], y = iris_binary$Species, plot = TRUE)
