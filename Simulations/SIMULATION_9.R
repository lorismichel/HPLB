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
require(e1071)
require(rpart)
require(class)


# summary function
rhoSummaries <- function(dataset = "Boston") {

  d <- getDataset(dataset = dataset)
  y.train <- d$y.train
  x.train <- d$x.train
  y.test  <- d$y.test
  x.test  <- d$x.test

  # fit models
  rf <- ranger::ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
  ld <- MASS::lda(y~., data.frame(y = y.train, x.train))
  gl <- glm(formula = y~., data =  data.frame(y = y.train, x.train), family = "binomial")
  #sv <- svm(formula = factor(y)~.,  data = data.frame(y = y.train, x.train), probability = TRUE)
  rp <- rpart(formula = y~., data =  data.frame(y = y.train, x.train))


  # summaries
  mat.metrics <- matrix(ncol=4,nrow=6)


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

  # GLM
  rho_glm <- predict(gl, newdata = data.frame(x.test), type = "response")
  mat.metrics[4,1] <- auc(y.test, rho_glm)
  mat.metrics[4,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_glm, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  mat.metrics[4,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_glm>0.5))
  mat.metrics[4,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_glm>0.5))))/sum(tab)

  # # SVM
  # rho_svm <- predict(sv, data = data.frame(x.test), probability = TRUE)
  # mat.metrics[4,1] <- auc(y.test, rho_sv)
  # mat.metrics[4,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_sv, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  # mat.metrics[4,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_sv>0.5))
  # mat.metrics[4,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_sv>0.5))))/sum(tab)

  # RPART
  rho_rp <- predict(rp, newdata = data.frame(x.test))[,"1"]
  mat.metrics[5,1] <- auc(y.test, rho_rp)
  mat.metrics[5,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_rp, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  mat.metrics[5,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_rp>0.5))
  mat.metrics[5,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_rp>0.5))))/sum(tab)

  # KNN
  kn <- knn(train = x.train, cl = y.train, test = x.test, k = 5, prob = TRUE)
  rho_knn <- ifelse(kn == 1, attr(kn, "prob"), 1-attr(kn, "prob"))
  mat.metrics[6,1] <- auc(y.test, rho_knn)
  mat.metrics[6,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho_knn, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
  mat.metrics[6,3] <- ModelMetrics::fScore(actual = as.numeric(y.test)-1, predicted = as.numeric(rho_knn>0.5))
  mat.metrics[6,4] <- sum(diag(tab <- table(as.numeric(y.test)-1, as.numeric(rho_knn>0.5))))/sum(tab)


  colnames(mat.metrics) <- c("AUC", "TV", "Fscore", "Accuracy")
  rownames(mat.metrics) <- c("RF", "H0", "LDA", "GLM", "RPART", "KNN")

  return(mat.metrics)
}

# dataset names
dataset.names <- c("Boston", "titanic", "BreastCancer", "Ionosphere", "abalone",
                   "banknotes", "Default", "credit")

# running the sims
res <- list()
for (n in dataset.names) {
  res$n <- rhoSummaries(dataset = n)
  print(n)
}


save(res,
     file = paste0(PATH.SAVE, "DATA_SIMULATION_9.Rdata"))


