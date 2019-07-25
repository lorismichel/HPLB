# SIMUlATION_4
# type: copula shift

## saving paths
PATH.SAVE = "../Data/"

## requirements
require(dWit)
require(parallel)
require(data.table)
require(stats)
require(ranger)
require(MASS)
require(reticulate)

# Source mmd
source('./helpers.R')


## data generation & bayes ratio
genMultiVar <- function(n, Sigma = diag(1, 3)) {
  return(mvrnorm(n = n, mu = rep(0, nrow(Sigma)), Sigma = Sigma))
}
genMultiVarT <- function(n, R = diag(1, 3), df = 3) {
  return(rtcopula(n = n, R = R, df = 3))
}



# running simulations on clusters
# params sims
cov1 <- diag(1, 3)
cov2 <- diag(0.5, 3) + matrix(0.5, 3, 3)
R <- diag(1, 10)
n <- 700
nrep <- 100
grid <- c(1:nrep)


# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(r) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
  rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(r) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
  rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res3 <- mcmapply(FUN = function(r) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
  #rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  #rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  f <- generateWitnessFunctionGaussianKernel(sample2 = x.train[1:n,], sample1 = x.train[-c(1:n),])
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = apply(x.test, 1, function(x) f(x)), s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res4 <- mcmapply(FUN = function(r) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rbind(genMultiVar(n = n, R), genMultiVarT(n = n, R))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rbind(genMultiVar(n = n, R), genMultiVarT(n = n, R))
  #rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  #rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  f <- generateWitnessFunctionGaussianKernel(sample2 = x.train[1:n,], sample1 = x.train[-c(1:n),])
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = apply(x.test, 1, function(x) f(x)), s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res5 <- mcmapply(FUN = function(r) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rbind(genMultiVar(n = n, R), genMultiVarT(n = n, R))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rbind(genMultiVar(n = n, R), genMultiVarT(n = n, R))
  #rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  #rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res6 <- mcmapply(FUN = function(r) {
  y.train <- factor(c(rep(0,n), rep(1,n)))
  x.train <- rbind(genMultiVar(n = n, R), genMultiVarT(n = n, R))
  y.test <- factor(c(rep(0,n), rep(1,n)))
  x.test <- rbind(genMultiVar(n = n, R), genMultiVarT(n = n, R))
  #rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  #rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

#set.seed(123, "L'Ecuyer")
#res4 <- mcmapply(FUN = function(r) {
#  y.train <- factor(c(rep(0,n), rep(1,n)))
#  x.train <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
#  y.test <- factor(c(rep(0,n), rep(1,n)))
#  x.test <- rbind(genMultiVar(n = n, cov1), genMultiVar(n = n, cov2))
#  #rf <- ranger(y~., data = data.frame(y = y.train, x = x.train),classification = TRUE, probability = TRUE)
3  #rho <- predict(rf, data = data.frame(x = x.test))$predictions[,"1"]
#  f <- generateWitnessFunctionGaussianKernel(sample2 = x.train[1:n,], sample1 = x.train[-c(1:n),])
#  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = apply(x.test, 1, function(x) f(x)), s = 0.5, estimator.type = "binomial")$tvhat
#  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
#  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.cor <- data.table(rep     = grid,
                          tvhat_search_rf = res1,
                          tvhat_search_mmd = res3,
                          tvhat_binomial_rf = res2)
power.table.r <- data.table(rep     = grid,
                          tvhat_search_rf = res5,
                          tvhat_search_mmd = res4,
                          tvhat_binomial_rf = res6)

# saving results of simulations
save(power.table.cor, power.table.r, file = paste0(PATH.SAVE, "DATA_SIMULATION_4.Rdata"))
