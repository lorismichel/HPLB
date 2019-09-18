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
dataProbPermutation <- function(y, yhat, up, down,
                                direction = "in", use.perm=FALSE) {
  y <- as.numeric(levels(y))[y]

  if (!use.perm) {
    if (direction == "in") {
      return(factor(ifelse(yhat >= down & yhat <= up,rbinom(n = length(y), size = 1, p = 1/2),
                           y)))
    } else {
      return(factor(ifelse(yhat < down | yhat > up, rbinom(n = length(y), size = 1, p = 1/2),
                           y)))
    }
  } else {
    if (direction == "in") {
      return(factor(ifelse(yhat >= down & yhat <= up, sample(y),
                           y)))
    } else {
      return(factor(ifelse(yhat < down | yhat > up, sample(y),
                           y)))
    }
  }


}

permuteRho <- function(rho, u, d) {
  ind <- which(rho <= u & rho >= d)
  rho.perm <- rho
  rho.perm[ind] <- sample(rho.perm[ind])
  return(rho.perm)
}



# Prepro
crime01 = rep(0, length(Boston$crim))
crime01[Boston$crim > median(Boston$crim)] = 1
Boston = data.frame(Boston)

train = 1:(dim(Boston)[1]*(1/3))
test = setdiff(1:dim(Boston)[1], train)
Boston.train = Boston[train, -c(1)]
Boston.test = Boston[test, -c(1)]
crim01.test = crime01[test]
crim01.train = crime01[train]

Boston.train <- Boston.train[rep(1:nrow(Boston.train), each=1),]
Boston.test <- Boston.test[rep(1:nrow(Boston.test), each=1),]
#response.test = factor(response[test])
crim01.test = factor(rep(crim01.test,each=1))
#response.test = factor(response[test])
crim01.train = factor(rep(crim01.train,each=1))


# running simulations on clusters
# params sims
grid.p <- seq(0, 0.5, by = 0.01)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- crim01.train
  x.train <- Boston.train
  y.test <- crim01.test
  x.test <- Boston.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- crim01.train
  x.train <- Boston.train
  y.test <- crim01.test
  x.test <- Boston.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
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

train = 1:(dim(titanic)[1]*(1/3))
test = setdiff(1:dim(titanic)[1], train)
titanic.train = titanic[train, -c(1)]
titanic.test = titanic[test, -c(1)]
survived.test = factor(survived[test])
survived.train = factor(survived[train])


#titanic.train <- titanic.train[rep(1:nrow(titanic.train), each=10),]
#titanic.test <- titanic.test[rep(1:nrow(titanic.test), each=10),]
#response.test = factor(response[test])
#survived.test = factor(rep(survived[test],each=10))
#response.test = factor(response[test])
#survived.train = factor(rep(survived[train],each=10))

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

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- survived.train
  x.train <- titanic.train
  y.test <- survived.test
  x.test <- titanic.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
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

train = 1:(dim(breast)[1]*(1/3))
test = setdiff(1:dim(breast)[1], train)
breast.train = breast[train, -c(1,10)]
breast.test = breast[test, -c(1,10)]
#response.test = factor(response[test])
#response.train = factor(response[train])

#breast.train <- breast.train[rep(1:nrow(breast.train), each=10),]
#breast.test <- breast.test[rep(1:nrow(breast.test), each=10),]
response.train = factor(response[train])
#response.test = factor(rep(response[test],each=10))
response.test = factor(response[test])
#response.train = factor(rep(response[train],each=10))
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

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search")$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- breast.train
  y.test <- response.test
  x.test <- breast.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
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

train = 1:(dim(iono)[1]*(1/3))
test = setdiff(1:dim(iono)[1], train)
iono.train = iono[train, -c(35)]
iono.test = iono[test, -c(35)]
#iono.train <- iono.train[rep(1:nrow(iono.train), each=10),]
#iono.test <- iono.test[rep(1:nrow(iono.test), each=10),]
#response.test = factor(response[test])
#response.test = factor(rep(response[test],each=10))
response.test = factor(response[test])
#response.train = factor(rep(response[train],each=10))
response.train = factor(response[train])
#


# running simulations on clusters
# params sims
grid.p <- rep(seq(0, 0.5, by = 0.01),each=1)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- iono.train
  y.test <- response.test
  x.test <- iono.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- iono.train
  y.test <- response.test
  x.test <- iono.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.iono <- data.table(p     = grid,
                                 tvhat_search = res1,
                                 tvhat_binomial = res2)

#power.table.i <- power.table.iono[,.(meanB = mean(tvhat_binomial), meanS = mean(tvhat_search), sdS = sd(tvhat_search),sdB = sd(tvhat_binomial)),by=c("p")]


# saving results of simulations
#save(power.table.boston, power.table.titanic, power.table.breast, power.table.iono, file = paste0(PATH.SAVE, "DATA_SIMULATION_5.Rdata"))




abalone.cols = c("sex", "length", "diameter", "height", "whole.wt",
                "shucked.wt", "viscera.wt", "shell.wt", "rings")

url <- 'http://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data'
abalone <- read.table(url, sep=",", row.names=NULL, col.names=abalone.cols,
                      nrows=4177)
colnames(abalone) <- abalone.cols


train = 1:(dim(abalone)[1]*(1/3))
test = setdiff(1:dim(abalone)[1], train)
abalone.train = abalone[train, -c(9)]
abalone.test = abalone[test, -c(9)]

# young and adult vs old
response <- ifelse(abalone$rings <=5, 0, ifelse(abalone$rings <= 13, 0, 1))
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

#


# running simulations on clusters
# params sims
grid.p <- rep(seq(0, 0.5, by = 0.01),each=1)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- abalone.train
  y.test <- response.test
  x.test <- abalone.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- abalone.train
  y.test <- response.test
  x.test <- abalone.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train), classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.abalone <- data.table(p     = grid,
                               tvhat_search = res1,
                               tvhat_binomial = res2)

# adult
adult <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data',
                    sep = ',', fill = F, strip.white = T)
colnames(adult) <- c('age', 'workclass', 'fnlwgt', 'educatoin',
                     'educatoin_num', 'marital_status', 'occupation', 'relationship', 'race', 'sex',
                     'capital_gain', 'capital_loss', 'hours_per_week', 'native_country', 'income')
levels(adult$income) <- c(0,1)

train = 1:(dim(adult)[1]*(1/3))
test = setdiff(1:dim(adult)[1], train)
adult.train = adult[train, -c(15)]
adult.test = adult[test, -c(15)]

# young and adult vs old
response <- adult$income
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

#


# running simulations on clusters
# params sims
grid.p <- rep(seq(0, 0.5, by = 0.01),each=1)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- adult.train
  y.test <- response.test
  x.test <- adult.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- adult.train
  y.test <- response.test
  x.test <- adult.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.adult <- data.table(p     = grid,
                                  tvhat_search = res1,
                                  tvhat_binomial = res2)



# bank notes
banknotes <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/00267/data_banknote_authentication.txt',
                    sep = ',', fill = F, strip.white = T)
banknotes <- banknotes[sample(1:nrow(banknotes)),]
colnames(banknotes) <- c("X1","X2","X3","X4","y")



train = 1:(dim(banknotes)[1]*(1/3))
test = setdiff(1:dim(banknotes)[1], train)
banknotes.train = banknotes[train, -c(5)]
banknotes.test = banknotes[test, -c(5)]

# young and adult vs old
response <- banknotes$y
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

#


# running simulations on clusters
# params sims
grid.p <- rep(seq(0, 0.5, by = 0.01),each=1)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- banknotes.train
  y.test <- response.test
  x.test <- banknotes.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- banknotes.train
  y.test <- response.test
  x.test <- banknotes.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial",  threshold = mean(rho.perm))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.banknotes <- data.table(p     = grid,
                                tvhat_search = res1,
                                tvhat_binomial = res2)


# Default
default <- ISLR::Default


train = 1:(dim(default)[1]*(1/3))
test = setdiff(1:dim(default)[1], train)
default.train = default[train, -c(1)]
default.test = default[test, -c(1)]

# young and adult vs old
response <- ifelse(default$default=="Yes",1,0)
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

#


# running simulations on clusters
# params sims
grid.p <- rep(seq(0, 0.5, by = 0.01),each=1)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- default.train
  y.test <- response.test
  x.test <- default.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,  threshold = mean(rho.perm), estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- default.train
  y.test <- response.test
  x.test <- default.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial", threshold = mean(rho.perm))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.default <- data.table(p    = grid,
                                    tvhat_search = res1,
                                    tvhat_binomial = res2)


# credit
url="http://freakonometrics.free.fr/german_credit.csv"
credit=read.csv(url, header = TRUE, sep = ",")
credit <- credit[sample(1:nrow(credit)),]

train = 1:(dim(credit)[1]*(1/3))
test = setdiff(1:dim(credit)[1], train)
credit.train = credit[train, -c(1)]
credit.test = credit[test, -c(1)]

# young and adult vs old
response <- credit$Creditability
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

# running simulations on clusters
# params sims
grid.p <- rep(seq(0, 0.5, by = 0.01),each=1)
grid <- grid.p

# running simulations
set.seed(123, "L'Ecuyer")
res1 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- credit.train
  y.test <- response.test
  x.test <- credit.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

set.seed(123, "L'Ecuyer")
res2 <- mcmapply(FUN = function(p) {
  y.train <- response.train
  x.train <- credit.train
  y.test <- response.test
  x.test <- credit.test

  # 1) fit a forest RF1
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  # 2) set boundaries
  rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  dq <- quantile(rho, p)
  uq <- quantile(rho, 1-p)
  # 3) permute labels based on RF1 (oob on train, pred on test)
  #y.train <- dataProbPermutation(y = y.train, yhat = rf$predictions[,"1"],  up = uq, down = dq)
  #y.test <- dataProbPermutation(y = y.test, yhat = predict(rf, data = data.frame(x.test))$predictions[,"1"], up = uq, down = dq)
  # 4) fit a new forest RF2 to get an ordering
  #rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  rho.perm <- permuteRho(rho = rho, u = uq, d = dq)
  # 5) evaluate
  tvhat <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "binomial", threshold = mean(rho.perm))$tvhat
  #tvhat <- dWit(t = inputs$t, rho = inputs$rho, s = 0.5, estimator.type = "binomial")$tvhat
  tvhat}, grid, mc.cores = 10)

# gathering data
power.table.credit <- data.table(p    = grid,
                                  tvhat_search = res1,
                                  tvhat_binomial = res2)

save(power.table.boston,
     power.table.titanic,
     power.table.breast,
     power.table.iono,
     power.table.abalone,
     power.table.adult,
     power.table.banknotes,
     power.table.default,
     power.table.credit,
     file = paste0(PATH.SAVE, "DATA_SIMULATION_5.Rdata"))
