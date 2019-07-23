# SIMUlATION_1
# type: checking power by a loglog plot of sample size against contamination/tv
#       comparison with missclassification error rate test.

## saving paths
PATH.SAVE = "../Data/"

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



# running simulations on clusters
# params sims
nrep <- 10
grid.gamma <- rep(seq(0.1,  1.2, by = 0.1), nrep)
grid.n <- rep(10^{2:5}, nrep)
grid <- expand.grid(grid.n, grid.gamma)

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(n, gamma) {
  lambda <- n^{-gamma}
  inputs <- generateInputData(n, lambda = lambda)
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  as.numeric(tvhat > 0)}, grid[,1], grid[,2], mc.cores = 25)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(n, gamma) {
  lambda <- n^{-gamma}
  inputs <- generateInputData(n, lambda = lambda)
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  as.numeric(tvhat > 0)}, grid[,1], grid[,2], mc.cores = 25)

# gathering data
power.data <- data.table(logn = log(grid[,1], base = 10),
			 loglambda = -log(grid[,1], base = 10)*grid[,2],
			 reject_search = res1, reject_binomial = res2)
power.table <- power.data[,.(power_search = mean(reject_search), power_binomial = mean(reject_binomial)),by=c("logn","loglambda")]

# saving results of simulations
save(power.data, power.table, file = paste0(PATH.SAVE, "DATA_SIMULATION_1_new.Rdata"))
