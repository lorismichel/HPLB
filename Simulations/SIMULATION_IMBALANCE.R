require(dWit)
require(distrEx)
require(ranger)

# Imbalance Example

n <- 1000
m <- 0.2*n
set.seed(0)

Y<- rnorm(n, mean=1)
X <- rnorm(m)


# True TV
lambda=TotalVarDist(Norm(mean = 0), Norm(mean = 1))



t = c(rep(1, n), rep(0, m))


# 1) Just using rho=1
rho1= rep(1, length(t))

tvhat1 <- dWit(t = t, rho = rho1, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhatB1 <- dWit(t = t, rho = rho1, s = 0.5, estimator.type = "binomial-test")$tvhat

# 2) Using wrongly specified Bayes
bayesRatio <- function(x, lambda) { dnorm(x, mean=1) / (dnorm(x, mean=1) + dnorm(x))} # (wrongly) assuming s=1/2
rho21 <- bayesRatio(Y)
rho22 <- bayesRatio(X)
rho2<-c(rho21, rho22)

tvhat2 <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhatB2 <- dWit(t = t, rho = rho2, s = 0.5, estimator.type = "binomial-test")$tvhat


# 3) Using RF
Ytest<- rnorm(n, mean=1)
Xtest <- rnorm(m)
rf <- ranger(t~x, data = data.frame(t=t, x = c(Y,X)),classification = TRUE, probability = TRUE)
rho3 <- predict(rf, data = data.frame(x = c(Ytest,Xtest)))$predictions[,1]

tvhat3 <- dWit(t = t, rho = rho3, s = 0.5, estimator.type = "hypergeometric-test")$tvhat
tvhatB3 <- dWit(t = t, rho = rho3, s = 0.5, estimator.type = "binomial-test")$tvhat

