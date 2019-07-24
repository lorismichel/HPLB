# SIMUlATION_5
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
require(MASS)
require(mlbench)


## data generation & bayes ratio
dataContamination <- function(y, p = 0) {
  y <- as.numeric(levels(y))[y]
  return(factor(ifelse(rbinom(n = length(y), size = 1, p = p) == 0, y, 1-y)))
}


# Prepro
crime01 = rep(0, length(Boston$crim))
crime01[Boston$crim > median(Boston$crim)] = 1
Boston = data.frame(Boston)

train = 1:(dim(Boston)[1]/2)
test = (dim(Boston)[1]/2 + 1):dim(Boston)[1]
Boston.train = Boston[train, -c(1)]
Boston.test = Boston[test, -c(1)]
crim01.test = factor(crime01[test])
crim01.train = factor(crime01[train])
# running simulations on clusters
# params sims
grid.p <- seq(0, 0.5, by = 0.025)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- crim01.train
  x.train <- Boston.train
  y.test <- crim01.test
  x.test <- Boston.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, Boston.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(Boston.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- crim01.train
  x.train <- Boston.train
  y.test <- crim01.test
  x.test <- Boston.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, Boston.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(Boston.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.boston <- data.table(p     = grid,
                          tvhat_search = res1,
                          tvhat_binomial = res2)



####################### titanic dataset

# PrePro
titanic <- read.csv(url("https://web.stanford.edu/class/archive/cs/cs109/cs109.1166/stuff/titanic.csv"),
         header = T)

survived = titanic$Survived
titanic = data.frame(titanic)

train = 1:(dim(titanic)[1]/2)
test = (dim(titanic)[1]/2 + 1):dim(titanic)[1]
titanic.train = titanic[train, -c(1)]
titanic.test = titanic[test, -c(1)]
survived.test = factor(survived[test])
survived.train = factor(survived[train])
#


# running simulations on clusters
# params sims
grid.p <- seq(0, 0.5, by = 0.01)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- survived.train
  x.train <- titanic.train
  y.test <- survived.test
  x.test <- titanic.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- survived.train
  x.train <- titanic.train
  y.test <- survived.test
  x.test <- titanic.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.titanic <- data.table(p     = grid,
                          tvhat_search = res1,
                          tvhat_binomial = res2)





####################### breast cancer dataset

# PrePro
data("BreastCancer")
breast <- na.omit(BreastCancer)

response <- factor(breast$Class)
levels(response) <- c(0,1)

train = 1:(dim(breast)[1]/2)
test = (dim(breast)[1]/2 + 1):dim(breast)[1]
breast.train = breast[train, -c(1,10)]
breast.test = breast[test, -c(1,10)]
response.test = factor(response[test])
response.train = factor(response[train])
#


# running simulations on clusters
# params sims
grid.p <- seq(0, 0.5, by = 0.01)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- breast.train
  y.test <- response.test
  x.test <- breast.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- breast.train
  y.test <- response.test
  x.test <- breast.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.breast <- data.table(p     = grid,
                                  tvhat_search = res1,
                                  tvhat_binomial = res2)




####################### ionosphere dataset

# PrePro
data("Ionosphere")
iono <- na.omit(Ionosphere)

response <- factor(iono$Class)
levels(response) <- c(0,1)

train = 1:(dim(iono)[1]/2)
test = (dim(iono)[1]/2 + 1):dim(iono)[1]
iono.train = iono[train, -c(35)]
iono.test = iono[test, -c(35)]
response.test = factor(response[test])
response.train = factor(response[train])
#


# running simulations on clusters
# params sims
grid.p <- seq(0, 0.5, by = 0.01)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- iono.train
  y.test <- response.test
  x.test <- iono.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- iono.train
  y.test <- response.test
  x.test <- iono.test
  y.train <- dataContamination(y.train, p)
  y.test  <- dataContamination(y.test, p)
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho, s = 0.5, estimator.type = "binomial")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.iono <- data.table(p     = grid,
                                 tvhat_search = res1,
                                 tvhat_binomial = res2)


# saving results of simulations
save(power.table.boston, power.table.titanic, power.table.breast, power.table.iono, file = paste0(PATH.SAVE, "DATA_SIMULATION_5.Rdata"))
