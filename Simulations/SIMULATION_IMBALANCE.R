require(dWit)
require(distrEx)
require(ranger)
require(MASS)

# Imbalance Example





# Do proper example with replications:
B<-50
svec<-seq(0.1, 0.5, by=0.1)
N<-1000
set.seed(0)

tvhat1<-c()
tvhatB1<-c()
tvhat2 <-c()
tvhatB2 <-c()
tvhat3 <- c()
tvhatB3 <- c()

silly<-matrix(NaN, nrow=length(svec), ncol=2)
ldamat<-matrix(NaN, nrow=length(svec), ncol=2)
RFmat<-matrix(NaN, nrow=length(svec), ncol=2)

j<-0

for (s in svec){
  j<-j+1

for (b in 1:B){

  n <- rbinom(n=1, size=N, prob=1-s)
  m <- N-n
  set.seed(0)

  Y<- rnorm(n, mean=1)
  X <- rnorm(m)
  Ytest<- rnorm(n, mean=1)
  Xtest <- rnorm(m)


  # True TV
  lambda=TotalVarDist(Norm(mean = 0), Norm(mean = 1))



  t = c(rep(1, n), rep(0, m))


  # 1) Just using rho=1
  rho1= rep(1, length(t))

  tvhat1[b] <- dWit(t = t, rho = rho1, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
  tvhatB1[b] <- dWit(t = t, rho = rho1, s = 0.5, estimator.type = "binomial-test")$tvhat



  # 2) Using LDA
  l <- lda(t~x, data = data.frame(t=t, x = c(Y,X)))
  rho2 <- unlist(predict(l, data = data.frame(x = c(Ytest,Xtest)))$posterior[,"1"])

  tvhat2[b] <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
  tvhatB2[b] <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "binomial-test")$tvhat


  # 3) Using RF

  rf <- ranger(t~x, data = data.frame(t=t, x = c(Y,X)),classification = TRUE, probability = TRUE)
  rho3 <- predict(rf, data = data.frame(x = c(Ytest,Xtest)))$predictions[,1]

  tvhat3[b] <- dWit(t = t, rho = rho3, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
  tvhatB3[b] <- dWit(t = t, rho = rho3, s = 0.5, estimator.type = "binomial-test")$tvhat


}



silly[j,1]<-round(mean(tvhat1),2)
silly[j,2]<-round(mean(tvhatB1),2)

ldamat[j,1] <- round(mean(tvhat2),2)
ldamat[j,2] <- round(mean(tvhatB2),2)


RFmat[j,1] <- round(mean(tvhat3),2)
RFmat[j,2] <- round(mean(tvhatB3),2)

}

# True TV
lambda=TotalVarDist(Norm(mean = 0), Norm(mean = 1))


# Toying around:

n <- 1000
m <- 0.2*n


Y<- rnorm(n, mean=1)
X <- rnorm(m)






t = c(rep(1, n), rep(0, m))


# 1) Just using rho=1
rho1= rep(1, length(t))

tvhat1 <- dWit(t = t, rho = rho1, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhatB1 <- dWit(t = t, rho = rho1, s = 0.5, estimator.type = "binomial-test")$tvhat

# 2) Using wrongly specified Bayes
## Why are we not overestimating?
bayesRatio <- function(x) {0.5* dnorm(x, mean=1) / (0.5*dnorm(x, mean=1) + 0.5 *dnorm(x))} # (wrongly) assuming s=1/2, so that the whole binomial story actually works.
rho21 <- bayesRatio(Y)
rho22 <- bayesRatio(X)
rho2<-c(rho21, rho22)

tvhat2 <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhat2 <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "asymptotic-tv-search")$tvhat
tvhatB2 <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "binomial-test")$tvhat


# 3) Using RF
Ytest<- rnorm(n, mean=1)
Xtest <- rnorm(m)
rf <- ranger(t~x, data = data.frame(t=t, x = c(Y,X)),classification = TRUE, probability = TRUE)
rho3 <- predict(rf, data = data.frame(x = c(Ytest,Xtest)))$predictions[,1]

tvhat3 <- dWit(t = t, rho = rho3, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhatB3 <- dWit(t = t, rho = rho3, s = 0.5, estimator.type = "binomial-test")$tvhat




# 4) Using LDA

Ytest<- rnorm(n, mean=1)
Xtest <- rnorm(m)
l <- lda(t~x, data = data.frame(t=t, x = c(Y,X)))
rho4 <- unlist(predict(l, data = data.frame(x = c(Ytest,Xtest)))$posterior[,"1"])

tvhat4 <- dWit(t = t, rho = rho4, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhatB4 <- dWit(t = t, rho = rho4, s = 0.5, estimator.type = "binomial-test")$tvhat








