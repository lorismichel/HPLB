# SIMUlATION_1
# type: checking power by a loglog plot of sample size against contamination/tv

## saving paths
PATH.SAVE = ""
## requirements
require(dWit)


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
inputs <- generateInputData(1000, lambda = 0.1)

grid.n <- 10^c(1,2,3,4,5,6,8,9,10)
gamma <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2)
res <- matrix(0,nrow=length(gamma),ncol=length(grid.n))
B <- 5
for (k in 1:B) {
for (i in 1:length(gamma)) {
  for (j in 1:length(grid.n)) {
    lambda <- grid.n[j]^{-gamma[i]}
    inputs <- generateInputData(grid.n[j], lambda = lambda)
    tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
    res[i,j] <- res[i,j] + as.numeric(tvhat > 0)/B
  }
}
  print(k)
}
mapply(FUN = function(n, gamma) {
  lambda <- n^{-gamma}
  inputs <- generateInputData(n, lambda = n^{-gamma})
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
  as.numeric(tvhat > 0)}, c(1000,1000,1000), c(0.5, 0.6, 0.6))


# running simulations on clusters
# params sims
nrep <- 5
grid.gamma <- rep(seq(0.1, 1, by = 0.1), nrep)
grid.n <- rep(10^{2:5},B)
grid <- expand.grid(grid.n, grid.gamma)

# running simulations
res <- mcmapply(FUN = function(n, gamma) {
  lambda <- n^{-gamma}
  inputs <- generateInputData(n, lambda = n^{-gamma})
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
  as.numeric(tvhat > 0)}, grid[,1], grid[,2])

# gathering data
power.data <- data.table(logn = log(grid[,1], base = 10), loglambda = -log(grid[,1], base = 10)*grid[,2], reject = res)
power.table <- power.data[,.(power = mean(reject)),by=c("logn","loglambda")]
plot(power.table$logn,power.table$loglambda,cex = power.table$power*2,pch=19)
abline(a = 0, b = -1,col="blue")

save(d, file = paste0(PATH.SAVE, "powerStudy.Rdata"))
