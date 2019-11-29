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
generateUnifMixturesData <- function(n, boundaries = matrix(c(-10,-9, -1, 0, 0, 1),ncol=2, byrow = TRUE), weights = c(0.1, 0.4, 0.5)) {

  force(boundaries)
  force(weights)

  # generate data
  ind.mix <- sample(c(1,2,3),size = n, prob = weights, replace = TRUE)
  x <- runif(n = n, min = boundaries[ind.mix, 1], max = boundaries[ind.mix, 2])

  return(x)
}

generateUnifMixturesDensities<- function(boundaries = matrix(c(-10,-9, -1, 0, 0, 1),ncol=2, byrow = TRUE), weights = c(0.1, 0.4, 0.5)) {

  force(boundaries)
  force(weights)

  # generate density function
  d <- function(y) sum(weights * sapply(1:length(weights), function(i) dunif(x = y, min = boundaries[i,1], max = boundaries[i,2])))

  return(d)
}

# running simulations on clusters


# params sims
nrep <- 50
grid.gamma <- rep(seq(0.1,  1.2, by = 0.1), nrep)
grid.n <- rep(10^{seq(2,5,length.out = 10)}, nrep)
grid <- expand.grid(grid.n, grid.gamma)


# scenario 1: bad for both methods

# running simulations
set.seed(123, "L'Ecuyer")
res_tv_search_1 <- mcmapply(FUN = function(n, gamma) {

  p1 <- n^{-gamma}
  p2 <- 0.5
  weights <- c(p1, p2)

  x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-2, -1, 1, 2) , ncol=2, byrow = TRUE), weights = weights)
  x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(9,10, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-10,-9, -2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
  d.right <- generateUnifMixturesDensities(boundaries = matrix(c(9,10, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))

  dWit(t = rep(0:1, each = n), rho = sapply(c(x1,x2), function(x) bayesRatio(x)), estimator.type = "asymptotic-tv-search")

  tvhat}, grid[,1], grid[,2], mc.cores = 25)

set.seed(123, "L'Ecuyer")
res_binomial_1 <- mcmapply(FUN = function(n, gamma) {

  p1 <- n^{-gamma}
  p2 <- 0.5
  weights <- c(p1, p2)

  x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-2, -1, 1, 2) , ncol=2, byrow = TRUE), weights = weights)
  x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(9,10, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-10,-9, -2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
  d.right <- generateUnifMixturesDensities(boundaries = matrix(c(9,10, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  dWit(t = rep(0:1, each = n), rho = sapply(c(x1,x2), function(x) bayesRatio(x)), estimator.type = "binomial-test")

  tvhat}, grid[,1], grid[,2], mc.cores = 25)


# gathering data
power.data <- data.table(logn = log(grid[,1], base = 10),
                         loglambda = -log(grid[,1], base = 10)*grid[,2],
                         tv_search = res_tv_search_1, tv_binomial = res_binomial_1)
power.table <- power.data[,.(power_search = mean(tv_search>0), power_binomial = mean(tv_binomial>0)),by=c("logn","loglambda")]

# saving results of simulations
save(power.data, power.table, file = paste0(PATH.SAVE, "DATA_SIMULATION_1_sc1.Rdata"))

