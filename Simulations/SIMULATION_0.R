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

generateInputData <- function(n, lambda, imbalance.level = 0) {
  if (imbalance.level == 0) {
    rho1 <- bayesRatio(sampleExpContamination(n = n, lambda = 0), lambda = lambda)
    rho2 <- bayesRatio(sampleExpContamination(n = n, lambda = lambda), lambda = lambda)
  } else if (imbalance.level == 1) {
   rho1 <- bayesRatio(sampleExpContamination(n = n, lambda = 0), lambda = lambda)
   rho2 <- bayesRatio(sampleExpContamination(n = n-sqrt(n), lambda = lambda), lambda = lambda)
  } else {
   rho1 <- bayesRatio(sampleExpContamination(n = n, lambda = 0), lambda = lambda)
   rho2 <- bayesRatio(sampleExpContamination(n = n-(n/2), lambda = lambda), lambda = lambda)
  }
  t = c(rep(1, length(rho1)), rep(0, length(rho2)))
  rho <- c(rho1, rho2)
  return(list(t = t, rho = rho))
}



# running simulations on clusters
# params sims
nrep <- 50
#grid.gamma <- rep(seq(0.1,  1.2, by = 0.1), nrep)
grid <- rep(10^{seq(2,5,length.out = 10)}, nrep)
grid.imbalance.level <- c(0,1,2)
grid <- expand.grid(grid,grid.imbalance.level)
#grid <- expand.grid(grid.n, grid.gamma)

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(n, gamma) {
  inputs <- generateInputData(n, lambda = 0)
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "asymptotic-tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid[,1], grid[,2], mc.cores = 25)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(n, gamma) {
  inputs <- generateInputData(n, lambda = 0)
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid[,1], grid[,2], mc.cores = 25)

set.seed(123, "L'Ecuyer")
res3 <- mcmapply(FUN = function(n, gamma) {
  inputs <- generateInputData(n, lambda = 0)
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "tv-search")$tvhat
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "hyper-tv-search", z = n)$tvhat
  tvhat}, grid[,1], grid[,2], mc.cores = 25)

# gathering data
power.data <- data.table(logn = log(grid[,1], base = 10),
			 imbalance = grid[,2],
                         tv_search = res1, tv_binomial = res2, tv_hyper = res3)
power.table <- power.data[,.(power_search = mean(tv_search>0), power_binomial = mean(tv_binomial>0)), power_hyper = mean(tv_hyper>0),by=c("logn","imbalance")]

# saving results of simulations
save(power.data, power.table, file = paste0(PATH.SAVE, "DATA_SIMULATION_0.Rdata"))
