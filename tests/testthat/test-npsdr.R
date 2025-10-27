set.seed(1)
n <- 200;
p <- 5;
x <- matrix(rnorm(n*p, 0, 2), n, p)
y <- 0.5*sqrt((x[,1]^2+x[,2]^2))*(log(x[,1]^2+x[,2]^2))+ 0.2*rnorm(n)
obj_kernel <- npsdr(x, y, plot=FALSE)
print(obj_kernel)
plot(obj_kernel)
