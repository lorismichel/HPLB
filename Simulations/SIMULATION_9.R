# SIMUlATION_9
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
require(pROC)



# summary function
rhoSummaries <- function(y.train, x.train, y.test, x.test) {

  # 1) fit models
  rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  ld <- MASS::lda(y~., data.frame(y = y.train, x.train))

  # summaries
  mat.metrics <- matrix(ncol=4,nrow=3)

  # RF
  rho_rf <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
  mat.metrics[1,1] <- auc(y.test, rho_rf)
  mat.metrics[1,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_rf, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  mat.metrics[1,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_rf>0.5))
  mat.metrics[1,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_rf>0.5))))/sum(tab)

  # H0
  rho_h0 <- rep(0, length(y.test))
  mat.metrics[2,1] <- auc(y.test, rho_h0)
  mat.metrics[2,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_h0, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  mat.metrics[2,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_h0>0.5))
  mat.metrics[2,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_h0>0.5))))/sum(tab)

  # LDA
  rho_lda <- predict(ld, newdata = data.frame(x.test))$posterior[,"1"]
  mat.metrics[3,1] <- auc(y.test, rho_lda)
  mat.metrics[3,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_lda, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  mat.metrics[3,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_lda>0.5))
  mat.metrics[3,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_lda>0.5))))/sum(tab)

  colnames(mat.metrics) <- c("AUC","TV","Fscore","Accuracy")
  rownames(mat.metrics) <- c("RF", "H0","LDA")

  return(mat.metrics)
}


# repro
set.seed(1)

# Prepro
crime01 = rep(0, length(Boston$crim))
crime01[Boston$crim > median(Boston$crim)] = 1
Boston = data.frame(Boston)
Boston <- Boston[sample(1:nrow(Boston), size = nrow(Boston), replace = FALSE),]

train = 1:(dim(Boston)[1]*(1/2))
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


mat.metrics.Boston <- rhoSummaries(y.train = crim01.train, x.train = Boston.train,
                                    y.test = crim01.test, x.test = Boston.test)

####################### titanic dataset

# PrePro
titanic <- read.csv(url("https://web.stanford.edu/class/archive/cs/cs109/cs109.1166/stuff/titanic.csv"),
                    header = T)

survived = titanic$Survived
titanic = data.frame(titanic)
titanic <- titanic[sample(1:nrow(titanic), replace = FALSE, size = nrow(titanic)),]
titanic$Sex <- ifelse(titanic$Sex=="male",0,1)
titanic <- titanic[,-c(3)]

train = 1:(dim(titanic)[1]*(1/2))
test = setdiff(1:dim(titanic)[1], train)
titanic.train = titanic[train, -c(1)]
titanic.test = titanic[test, -c(1)]
survived.test = factor(survived[test])
survived.train = factor(survived[train])


mat.metrics.titanic <- rhoSummaries(y.train = survived.train, x.train = titanic.train,
                                    y.test = survived.test, x.test = titanic.test)


####################### breast cancer dataset

# PrePro
data("BreastCancer")
breast <- na.omit(BreastCancer)
breast <- breast[sample(1:nrow(breast), replace = FALSE, size = nrow(breast)),]

response <- factor(breast$Class)
levels(response) <- c(0,1)

train = 1:(dim(breast)[1]*(1/2))
test = setdiff(1:dim(breast)[1], train)
breast.train = breast[train, -c(1,11)]
breast.test = breast[test, -c(1,11)]
#response.test = factor(response[test])
#response.train = factor(response[train])

#breast.train <- breast.train[rep(1:nrow(breast.train), each=10),]
#breast.test <- breast.test[rep(1:nrow(breast.test), each=10),]
response.train = factor(response[train])
#response.test = factor(rep(response[test],each=10))
response.test = factor(response[test])
#response.train = factor(rep(response[train],each=10))
#

mat.metrics.breast <- rhoSummaries(y.train = response.train, x.train = breast.train,
                                     y.test = response.test, x.test = breast.test)



####################### ionosphere dataset

# PrePro
data("Ionosphere")
iono <- na.omit(Ionosphere)

response <- factor(iono$Class)
levels(response) <- c(0,1)
iono <- iono[sample(1:nrow(iono), replace = FALSE, size = nrow(iono)),]


train = 1:(dim(iono)[1]*(1/2))
test = setdiff(1:dim(iono)[1], train)
iono.train = iono[train, -c(2,35)]
iono.test = iono[test, -c(2,35)]
#iono.train <- iono.train[rep(1:nrow(iono.train), each=10),]
#iono.test <- iono.test[rep(1:nrow(iono.test), each=10),]
#response.test = factor(response[test])
#response.test = factor(rep(response[test],each=10))
response.test = factor(response[test])
#response.train = factor(rep(response[train],each=10))
response.train = factor(response[train])

mat.metrics.iono <- rhoSummaries(y.train = response.train, x.train = iono.train,
                                    y.test = response.test, x.test = iono.test)




abalone.cols = c("sex", "length", "diameter", "height", "whole.wt",
                 "shucked.wt", "viscera.wt", "shell.wt", "rings")

url <- 'http://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data'
abalone <- read.table(url, sep=",", row.names=NULL, col.names=abalone.cols,
                      nrows=4177)
colnames(abalone) <- abalone.cols
abalone <- abalone[sample(1:nrow(abalone), replace = FALSE, size = nrow(abalone)),]


train = 1:(dim(abalone)[1]*(1/2))
test = setdiff(1:dim(abalone)[1], train)
abalone.train = abalone[train, -c(9)]
abalone.test = abalone[test, -c(9)]

# young and adult vs old
response <- ifelse(abalone$rings <=5, 0, ifelse(abalone$rings <= 13, 0, 1))
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

mat.metrics.abalone <- rhoSummaries(y.train = response.train, x.train = abalone.train,
                                    y.test = response.test, x.test = abalone.test)

# # adult
# adult <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data',
#                     sep = ',', fill = F, strip.white = T)
# colnames(adult) <- c('age', 'workclass', 'fnlwgt', 'educatoin',
#                      'educatoin_num', 'marital_status', 'occupation', 'relationship', 'race', 'sex',
#                      'capital_gain', 'capital_loss', 'hours_per_week', 'native_country', 'income')
# levels(adult$income) <- c(0,1)
# adult <- adult[sample(1:nrow(adult), replace = FALSE, size = nrow(adult)),]
#
#
# train = 1:(dim(adult)[1]*(1/2))
# test = setdiff(1:dim(adult)[1], train)
# adult.train = adult[train, -c(15)]
# adult.test = adult[test, -c(15)]
#
# # young and adult vs old
# response <- adult$income
# #response.test = factor(response[test])
# response.test = factor(response[test])
# response.train = factor(response[train])
#
# mat.metrics.adult <- rhoSummaries(y.train = response.train, x.train = adult.train,
#                                      y.test = response.test, x.test = adult.test)


# bank notes
banknotes <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/00267/data_banknote_authentication.txt',
                        sep = ',', fill = F, strip.white = T)
banknotes <- banknotes[sample(1:nrow(banknotes)),]
colnames(banknotes) <- c("X1","X2","X3","X4","y")
banknotes <- banknotes[sample(1:nrow(banknotes), replace = FALSE, size = nrow(banknotes)),]


train = 1:(dim(banknotes)[1]*(1/2))
test = setdiff(1:dim(banknotes)[1], train)
banknotes.train = banknotes[train, -c(5)]
banknotes.test = banknotes[test, -c(5)]

# young and adult vs old
response <- banknotes$y
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

mat.metrics.banknotes <- rhoSummaries(y.train = response.train, x.train = banknotes.train,
                                     y.test = response.test, x.test = banknotes.test)

# Default
default <- ISLR::Default
default <- default[sample(1:nrow(default), replace = FALSE, size = nrow(default)),]

train = 1:(dim(default)[1]*(1/2))
test = setdiff(1:dim(default)[1], train)
default.train = default[train, -c(1)]
default.test = default[test, -c(1)]

# young and adult vs old
response <- ifelse(default$default=="Yes",1,0)
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

y.train <- response.train
x.train <- default.train
y.test <- response.test
x.test <- default.test


mat.metrics.default <- rhoSummaries(y.train = response.train, x.train = default.train,
                                       y.test = response.test, x.test = default.test)



# credit
url="http://freakonometrics.free.fr/german_credit.csv"
credit=read.csv(url, header = TRUE, sep = ",")
credit <- credit[sample(1:nrow(credit), replace = FALSE, size = nrow(credit)),]

train = 1:(dim(credit)[1]*(1/2))
test = setdiff(1:dim(credit)[1], train)
credit.train = credit[train, -c(1)]
credit.test = credit[test, -c(1)]

# young and adult vs old
response <- credit$Creditability
#response.test = factor(response[test])
response.test = factor(response[test])
response.train = factor(response[train])

y.train <- response.train
x.train <- credit.train
y.test <- response.test
x.test <- credit.test

mat.metrics.credit <- rhoSummaries(y.train = response.train, x.train = credit.train,
                                     y.test = response.test, x.test = credit.test)


save(mat.metrics.Boston,
     mat.metrics.titanic,
     mat.metrics.breast,
     mat.metrics.iono,
     mat.metrics.abalone,
   #  mat.metrics.adult,
     mat.metrics.banknotes,
     mat.metrics.default,
     mat.metrics.credit,
     file = paste0(PATH.SAVE, "DATA_SIMULATION_9.Rdata"))


