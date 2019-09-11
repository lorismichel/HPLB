#
require(distrEx)
M <- 10
lambda <- 0.15

TotalVarDist(Norm(mean = 0), distr::UnivarMixingDistribution(Norm(mean=0), Norm(mean=M), mixCoeff = c(1-lambda, lambda)))

m <- 400
n <- 200
m.test <- 400
n.test <- 100
t <-c(rep(0,n), rep(1, m))
x <- c(rnorm(n,mean=0),ifelse(runif(m)<=lambda, rnorm(m,mean=M), rnorm(m, mean=0)))
t.test <-c(rep(0,n.test), rep(1, m.test))
x.test <- c(rnorm(n.test,mean=0),ifelse(runif(m.test)<=lambda, rnorm(m.test,mean=M), rnorm(m.test, mean=0)))
train <- data.frame(x=x,t=t)
test <- data.frame(x=x.test,t=t.test)
rf <- ranger::ranger(factor(t)~., train, probability = TRUE, class.weights = c(2/3,1/3))
rho <- predict(rf, test)$predictions[,2]
dWit::dWit(t = t.test, rho = rho, estimator.type = "binomial", s = 0.5, verbose.plot = TRUE)
dWit::dWit(t = t.test, rho = rho, estimator.type = "asymptotic-tv-search", s = 0.5,verbose.plot = TRUE)
dWit::dWit(t = t.test, rho = rho, estimator.type = "empirical-tv-search", s = 0.5)


# Outlier example
require(distrEx)
M <- -4
lambda <- 0.1

TotalVarDist(Norm(mean = 0), distr::UnivarMixingDistribution(Norm(mean=0), Norm(mean=M), mixCoeff = c(1-lambda, lambda)))

m <- 400
n <- 400
m.test <- 400
n.test <- 400
t <-c(rep(0,n), rep(1, m))
#x <- c(rnorm(n,mean=0),ifelse(runif(m)<=lambda, rnorm(m,mean=M), rnorm(m, mean=0)))
x<-c(rnorm(n,mean=0),ifelse(runif(m)<=lambda, runif(m,M,M+1), rnorm(m,mean=0)   ))
t.test <-c(rep(0,n.test), rep(1, m.test))
#x.test <- c(rnorm(n.test,mean=0),ifelse(runif(m.test)<=lambda, rnorm(m.test,mean=M), rnorm(m.test, mean=0)))
outliers<-runif(m)<=lambda
x.test <- c(rnorm(n,mean=0),ifelse(outliers, runif(m,M,M+1), rnorm(m,mean=0)  ))
nrofoutliers <- sum(outliers==1)
train <- data.frame(x=x,t=t)
test <- data.frame(x=x.test,t=t.test)
rf <- ranger::ranger(factor(t)~., train, probability = TRUE)
rho <- predict(rf, test)$predictions[,2]
dWit::dWit(t = t.test, rho = rho, estimator.type = "binomial", s = 0.5, verbose.plot = TRUE)
dWit::dWit(t = t.test, rho = rho, estimator.type = "asymptotic-tv-search", s = 0.5,verbose.plot = TRUE)
#dWit::dWit(t = t.test, rho = rho, estimator.type = "asymptotic-wit-search",direction="right", s = 0.5,verbose.plot = TRUE)
dWit::dWit(t = t.test, rho = rho, estimator.type = "empirical-tv-search", s = 0.5)
