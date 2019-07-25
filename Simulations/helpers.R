# helpers functions
require(mvtnorm)

generateWitnessFunctionGaussianKernel <- function(sample1, sample2, rho=NULL) {

    # Get median heuristic


  if (is.null(rho)){
  aggsamp<-as.matrix(rbind(sample1,sample2))
  G <- apply( aggsamp^2,1,function(y) sum(y))
  Q <- matrix(rep(G,dim(aggsamp)[1]), ncol=dim(aggsamp)[1])
  dists <- Q + t(Q) - 2*aggsamp%*%t(aggsamp)
  dists <- dists[upper.tri(dists, diag = FALSE)]
  rho <- sqrt(0.5*median(dists[dists > 0]))
  }


  f1 <- function(x) mean(apply(sample1,1,function(y) exp(-sum((y-x)^2)/rho)))
  f2 <- function(x) mean(apply(sample2,1,function(y) exp(-sum((y-x)^2)/rho)))

  return(function(x) f1(x)-f2(x))
}

# small example
#sample1 <- matrix(rnorm(300, sd = 2),ncol=1)
#sample2 <- matrix(rnorm(300),ncol=1)
#f <- generateWitnessFunctionGaussianKernel(sample1, sample2, rho = 2)
#plot(seq(-4,4,0.1), sapply(seq(-4,4,0.1), function(x) f(x)),type="b",pch=19)


rtcopula <- function(n,R,df){

  p <- ncol(R)
  T <- matrix(NA, n, p)

  for(i in 1:n){
    r <- rmvt(1,R,df)
    for(j in 1:p){
      term1 <- pt(r[j], df)
      term2 <- qnorm(term1)

      T[i,j] <- term2
    }
  }

  return(T)

}
