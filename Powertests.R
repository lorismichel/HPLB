## requirements
require(dWit)
require(parallel)
require(data.table)

## data generation & bayes ratio
sampleExpContamination <- function(n, lambda) {
  return(ifelse(sample(c(0,1), prob = c(lambda,1-lambda), size = n, replace = TRUE)==0, -rexp(n), rexp(n)))
}
bayesRatio <- function(x, lambda) { dexp(x) / (lambda*dexp(-x)+(1-lambda)*dexp(x) + dexp(x))}

generateInputData <- function(n, lambda) {
  rho1 <- bayesRatio(sampleExpContamination(n = n, lambda = 0), lambda = lambda)
  rho2 <- bayesRatio(sampleExpContamination(n = n, lambda = lambda), lambda = lambda)
  t = c(rep(1, length(rho1)), rep(0, length(rho2)))
  rho <- c(rho1, rho2)
  return(list(t = t, rho = rho))
}

n<-1e2
gamma<- -1/2
for (i in 1:100) {
  lambda=n^(gamma)
  d <- generateInputData(n = n, lambda = lambda)
  tvsearchest<-dWit(t = d$t, rho = d$rho, s = 0.5, estimator.type = "empirical-tv-search", tv.seq = seq(0,1,lambda/10))$tvhat
  print(paste("TVSearch:", tvsearchest))
  print(paste("Binomial:",dWit(t = d$t, rho = d$rho, s = 0.5, estimator.type = "binomial")$tvhat))
  #print(lambda/tvsearchest)

  }
