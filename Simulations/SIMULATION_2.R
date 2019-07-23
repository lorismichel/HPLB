# SIMUlATION_1
# type: checking power by a loglog plot of sample size against contamination/tv
#       comparison with missclassification error rate test.

## saving paths
PATH.SAVE = "../Data/"

## requirements
require(dWit)
require(parallel)
require(data.table)
require(stats)
require(ranger)


## data generation & bayes ratio
dataContamination <- function(y, p = 0) {
  y <- as.numeric(levels(y))[y]
  return(factor(ifelse(rbinom(n = length(y), size = 1, p = p) == 0, y, 1-y)))
}



# running simulations on clusters
# params sims
grid.p <- seq(0, 0.5, by = 0.025)
grid.mean <- c(0.5, 2, 4, 10)
grid <- expand.grid(grid.mean, grid.p)
n <- 10000

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(m, p) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- ifelse(y.train == 0, rnorm(n, 0), rnorm(n, m))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- ifelse(y.test == 0, rnorm(n, 0), rnorm(n, m))
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~x, data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid[,1], grid[,2], mc.cores = 20)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(m, p) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- ifelse(y.train == 0, rnorm(n, 0), rnorm(n, m))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- ifelse(y.test == 0, rnorm(n, 0), rnorm(n, m))
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~x, data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid[,1], grid[,2], mc.cores = 20)

# gathering data
power.table <- data.table(m     = grid[,1], 
			  p     = grid[,2], 
			  tvhat_search = res1,
                          tvhat_binomial = res2)

# saving results of simulations
save(power.table, file = paste0(PATH.SAVE, "DATA_SIMULATION_2.Rdata"))
