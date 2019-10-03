# SIMULATION 8


# Default
default <- ISLR::Default

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

set.seed(0)
rf <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)
rf4 <- ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE, class.weights = c(0.965 , 0.035 ))
ld <- MASS::lda(y~., data.frame(y = y.train, x.train))

# different types of rho
rho1 <- predict(rf, data = data.frame(x.test))$predictions[,"1"]
rho2 <- rep(0, length(y.test))
rho3 <- predict(ld, data = data.frame(x.test))$posterior[,"1"]
rho4 <- predict(rf4, data = data.frame(x.test))$predictions[,"1"]

mat.metrics <- matrix(ncol=4,nrow=4)

# relevant measures
require(pROC)
# rho1
mat.metrics[1,1] <- auc(response.test, rho1)
mat.metrics[1,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho1, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
mat.metrics[1,3] <- ModelMetrics::fScore(actual = as.numeric(response.test)-1, predicted = as.numeric(rho1>0.5))
mat.metrics[1,4] <- sum(diag(tab <- table(as.numeric(response.test)-1, as.numeric(rho1>0.5))))/sum(tab)

# rho2
mat.metrics[2,1] <-auc(response.test, rho2)
mat.metrics[2,2] <-dWit(t = as.numeric(levels(y.test))[y.test], rho = rho2, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
mat.metrics[2,3] <-ModelMetrics::fScore(actual = as.numeric(response.test)-1, predicted = rho2)
mat.metrics[2,4] <-sum(diag(tab <- table(as.numeric(response.test)-1, rho2)))/sum(tab)

# rho3
mat.metrics[3,1] <- auc(response.test, rho3)
mat.metrics[3,2] <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho3, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
mat.metrics[3,3] <- ModelMetrics::fScore(actual = as.numeric(response.test)-1, predicted = as.numeric(rho3>0.5))
mat.metrics[3,4] <- sum(diag(tab <- table(as.numeric(response.test)-1, as.numeric(rho3>0.5))))/sum(tab)

# rho4
mat.metrics[4,1] <- auc(response.test, rho4)
mat.metrics[4,2] <-dWit(t = as.numeric(levels(y.test))[y.test], rho = rho4, s = 0.5, estimator.type = "asymptotic-tv-search", tv.seq = seq(from = 0, to = 1, by = 0.001))$tvhat
mat.metrics[4,3] <-ModelMetrics::fScore(actual = as.numeric(response.test)-1, predicted = as.numeric(rho4>0.5))
mat.metrics[4,4] <-sum(diag(tab <- table(as.numeric(response.test)-1, as.numeric(rho4>0.5))))/sum(tab)

colnames(mat.metrics) <- c("AUC","TV","Fscore","Accuracy")
rownames(mat.metrics) <- c("RF", "H0","LDA", "RFw")

# resulting measures
print(mat.metrics)
