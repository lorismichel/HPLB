# SIMULATION 0

# different types of level checks
# 0) RF
# 1) LDA
# 2) most likely class
# 3) prior prob
# 4) 0.5

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

generateInputData <- function(n, lambda, imbalance.level = 0, method = 0) {

  if (imbalance.level == 0) {

    if (method == 0) {
      # RF
      x1.train <- sampleExpContamination(n = n, lambda = 0)
      x0.train <- sampleExpContamination(n = n, lambda = lambda)
      rf <- ranger::ranger(factor(t)~., data.frame(t = c(rep(1, n),
                                                         rep(0, n)),
                                                   x = c(x1.train,x0.train)), probability = TRUE)
      rho1 <- predict(rf, data.frame(x=sampleExpContamination(n = n, lambda = 0)))$predictions[,"1"]
      rho2 <- predict(rf, data.frame(x=sampleExpContamination(n = n, lambda = lambda)))$predictions[,"1"]
    } else if (method == 1) {
      # LDA
      x1.train <- sampleExpContamination(n = n, lambda = 0)
      x0.train <- sampleExpContamination(n = n, lambda = lambda)
      ld <- MASS::lda(factor(t)~., data.frame(t = c(rep(1, n),
                                                         rep(0, n)),
                                                   x = c(x1.train,x0.train)))
      rho1 <- predict(ld, data.frame(x=sampleExpContamination(n = n, lambda = 0)))$posterior[,"1"]
      rho2 <- predict(ld, data.frame(x=sampleExpContamination(n = n, lambda = lambda)))$posterior[,"1"]
    } else if (method == 2) {
      # most proeminent
      rho1 <- rep(0, n)
      rho2 <- rep(0, n)
    } else if (method == 3) {
      # prior
      rho1 <- rep(0.5, n)
      rho2 <- rep(0.5, n)
    } else {
      # 1/2
      rho1 <- rep(0.5, n)
      rho2 <- rep(0.5, n)
    }
  } else if (imbalance.level == 1) {
    if (method == 0) {
      # RF
      x1.train <- sampleExpContamination(n = n, lambda = 0)
      x0.train <- sampleExpContamination(n = n-sqrt(n), lambda = lambda)
      rf <- ranger::ranger(factor(t)~., data.frame(t = c(rep(1, n),
                                                         rep(0, n-sqrt(n))),
                                                   x = c(x1.train,x0.train)), probability = TRUE)
      rho1 <- predict(rf, data.frame(x=sampleExpContamination(n = n, lambda = 0)))$predictions[,"1"]
      rho2 <- predict(rf, data.frame(x=sampleExpContamination(n = n-sqrt(n), lambda = lambda)))$predictions[,"1"]
    } else if (method == 1) {
      # LDA
      x1.train <- sampleExpContamination(n = n, lambda = 0)
      x0.train <- sampleExpContamination(n = n-sqrt(n), lambda = lambda)
      ld <- MASS::lda(factor(t)~., data.frame(t = c(rep(1, n),
                                                    rep(0, n-sqrt(n))),
                                              x = c(x1.train,x0.train)))
      rho1 <- predict(ld, data.frame(x=sampleExpContamination(n = n, lambda = 0)))$posterior[,"1"]
      rho2 <- predict(ld, data.frame(x=sampleExpContamination(n = n-sqrt(n), lambda = lambda)))$posterior[,"1"]
    } else if (method == 2) {
      # most proeminent
      rho1 <- rep(0, n)
      rho2 <- rep(0, n-sqrt(n))
    } else if (method == 3) {
      # prior
      rho1 <- rep(n / (2*n-sqrt(n)), n)
      rho2 <- rep(n / (2*n-sqrt(n)), n-sqrt(n))
    } else {
      # 1/2
      rho1 <- rep(0.5, n)
      rho2 <- rep(0.5, n-sqrt(n))
    }
  } else {
    if (method == 0) {
      # RF
      x1.train <- sampleExpContamination(n = n, lambda = 0)
      x0.train <- sampleExpContamination(n = n-n/2, lambda = lambda)
      rf <- ranger::ranger(factor(t)~., data.frame(t = c(rep(1, n),
                                                         rep(0, n-n/2)),
                                                   x = c(x1.train,x0.train)), probability = TRUE)
      rho1 <- predict(rf, data.frame(x=sampleExpContamination(n = n, lambda = 0)))$predictions[,"1"]
      rho2 <- predict(rf, data.frame(x=sampleExpContamination(n = n-n/2, lambda = lambda)))$predictions[,"1"]
    } else if (method == 1) {
      # LDA
      x1.train <- sampleExpContamination(n = n, lambda = 0)
      x0.train <- sampleExpContamination(n = n-n/2, lambda = lambda)
      ld <- MASS::lda(factor(t)~., data.frame(t = c(rep(1, n),
                                                    rep(0, n-n/2)),
                                              x = c(x1.train,x0.train)))
      rho1 <- predict(ld, data.frame(x=sampleExpContamination(n = n, lambda = 0)))$posterior[,"1"]
      rho2 <- predict(ld, data.frame(x=sampleExpContamination(n = n-n/2, lambda = lambda)))$posterior[,"1"]
    } else if (method == 2) {
      # most proeminent
      rho1 <- rep(0, n)
      rho2 <- rep(0, n-n/2)
    } else if (method == 3) {
      # prior
      rho1 <- rep(n / (2*n-n/2), n)
      rho2 <- rep(n / (2*n-n/2), n-n/2)
    } else {
      # 1/2
      rho1 <- rep(0.5, n)
      rho2 <- rep(0.5, n-n/2)
    }
  }
  #
  t = c(rep(1, length(rho1)), rep(0, length(rho2)))
  rho <- c(rho1, rho2)
  return(list(t = t, rho = rho))
}



# running simulations on clusters
# params sims
nrep <- 50
grid <- rep(10^{seq(2,5,length.out = 10)}, nrep)
grid.imbalance.level <- c(0:2)
grid.method <- c(0:4)
grid <- expand.grid(grid,grid.imbalance.level, grid.method)


# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(n, imb, meth) {
  inputs <- generateInputData(n, lambda = 0, imb, meth)
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "asymptotic-tv-search")$tvhat
  tvhat}, grid[,1], grid[,2], grid[,3], mc.cores = 25)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(n, imb, meth) {
  inputs <- generateInputData(n, lambda = 0, imb, meth)
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial-test")$tvhat
  tvhat}, grid[,1], grid[,2], grid[,3], mc.cores = 25)

set.seed(123, "L'Ecuyer")
res3 <- mcmapply(FUN = function(n, imb, meth) {
  inputs <- generateInputData(n, lambda = 0, imb, meth)
  tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "hypergeometric-test", z = n)$tvhat
  tvhat}, grid[,1], grid[,2], grid[,3], mc.cores = 25)

# gathering data
power.data <- data.table(logn = log(grid[,1], base = 10),
			 imbalance = grid[,2],
			 method = grid[,3],
                         tv_search = res1, tv_binomial = res2, tv_hyper = res3)
power.table <- power.data[,.(power_search = mean(tv_search>0), power_binomial = mean(tv_binomial>0), power_hyper = mean(tv_hyper>0)),by=c("logn","imbalance", "method")]

# saving results of simulations
save(power.data, power.table, file = paste0(PATH.SAVE, "DATA_SIMULATION_0.Rdata"))
