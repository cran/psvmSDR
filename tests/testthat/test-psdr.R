## ----------------------------
## Linear PM
## ----------------------------
set.seed(1)
n <- 300; p <- 5;
x <- matrix(rnorm(n*p, 0, 2), n, p)
y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2*rnorm(n)
y.tilde <- sign(y)
obj <- psdr(x, y)
print(obj)
plot(obj, d=4)

## --------------------------------
## Linear PM (Binary classification)
## --------------------------------
obj_wsvm <- psdr(x, y.tilde, loss="wsvm")
print(obj_wsvm)
plot(obj_wsvm)

## ----------------------------
## User-defined loss function
## ----------------------------
mylogistic <- function(u) log(1+exp(-u))
psdr(x, y, loss="mylogistic")
