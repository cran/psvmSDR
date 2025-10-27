p <- 5
m <- 500 # batch size
N <- 10  # number of batches
obj <- NULL
for (iter in 1:N){
  set.seed(iter)
  x <- matrix(rnorm(m*p), m, p)
  y <-  x[,1]/(0.5 + (x[,2] + 1)^2) + 0.2 * rnorm(m)
  obj <- rtpsdr(x = x, y = y, obj=obj)
}
print(obj)
