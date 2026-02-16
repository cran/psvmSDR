set.seed(1)
p <- 5; m <- 300; B <- 3
obj <- NULL
for (b in 1:B) {
  x <- matrix(rnorm(m*p), m, p)
  y <- x[,1]/(0.5+(x[,2]+1)^2) + 0.2*rnorm(m)
  obj <- rtpsdr(x, y, obj=obj, h=8, lambda=1)
}
print(obj)
summary(obj)
