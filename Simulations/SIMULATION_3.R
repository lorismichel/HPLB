# SIMUlATION_3
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


# running simulations on clusters
# params sim
nrep <- 50
grid.mean <- c(rep(0.5, nrep), rep(2, nrep), rep(4, nrep), rep(10, nrep))
n <- 10000

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(m, p) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rnorm(n, m)
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rnorm(n, m)
  rf <- ranger(y~x, data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid.mean, mc.cores = 20)

# running simulations
set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(m, p) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rnorm(n, m)
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rnorm(n, m)
  rf <- ranger(y~x, data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid.mean, mc.cores = 20)

# gathering data
power.table <- data.table(m     = grid.mean, 
			  tvhat_search = res1,
                          tvhat_binomial = res2)

# saving results of simulations
save(power.table, file = paste0(PATH.SAVE, "DATA_SIMULATION_3.Rdata"))
