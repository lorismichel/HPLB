# helpers functions

generateWitnessFunctionGaussianKernel <- function(sample1, sample2, rho) {
  f1 <- function(x) mean(apply(sample1,1,function(y) exp(-sum((y-x)^2)/rho)))
  f2 <- function(x) mean(apply(sample2,1,function(y) exp(-sum((y-x)^2)/rho)))

  return(function(x) f1(x)-f2(x))
}

# small example
sample1 <- matrix(rnorm(300),ncol=1)
sample2 <- matrix(rnorm(300,sd=2),ncol=1)
f <- generateWitnessFunctionGaussianKernel(sample1, sample2, rho = 2)
plot(seq(-4,4,0.1), sapply(seq(-4,4,0.1), function(x) f(x)),type="b",pch=19)

