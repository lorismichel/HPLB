# test file for dWit function

# univariate regression
set.seed(0)
n <- 2000
t <- runif(n)
x <- ifelse(t <= 0.5,rnorm(n),rnorm(n,2))
rho = dnorm(x,2) / (dnorm(x) + dnorm(x,2))
par(mfrow=c(2,1))
plot(t,x,pch=19)
plot(t,rho,pch=19)

require(distrEx)
distrEx::TotalVarDist(Norm(mean = 0,sd = 1),Norm(mean = 2, sd = 1))

# single split
dWit(t = t, rho = rho, s = 0.5, alpha = 0.05, estimator.type = "basic")

# multiple splits
s <- c(0.3, 0.5, 0.7)
rho <- sapply(s, function(s) (1-s)*dnorm(x,2) / (s*dnorm(x) + (1-s)*dnorm(x,2)))
d <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "basic")
par(mfrow=c(2,1))
plot(s, d$tvhat.Fs, pch=19, type="b")
lines(s, d$tvhat.Gs, pch=19, type="b",col="blue")
plot(s, d$lambdahat.Fs, pch=19, type="b")
lines(s, d$lambdahat.Gs, pch=19, type="b",col="blue")
d

# tv-search:
s <- c(0.3, 0.5, 0.7)
rho = dnorm(x,2) / (dnorm(x) + dnorm(x,2))
dsearchTV <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "tv-search")
dsearchTV


#debug(dWit)
s <- c(0.3, 0.5, 0.7)
rho <- sapply(s, function(s) (1-s)*dnorm(x,2) / (s*dnorm(x) + (1-s)*dnorm(x,2)))
dsearchWit <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "wit-search")
dsearchWit


# univariate classification (binomial wins)
set.seed(0)
n <- 50000
t <- sample(c(0,1), n, replace = TRUE)
x <- ifelse(t <= 0, rnorm(n), rnorm(n, mean = 1))
rho = dnorm(x, mean = 1) / (dnorm(x) + dnorm(x, mean = 1))
par(mfrow=c(2,1))
#plot(t, x, pch=19)
#plot(t, rho, pch=19)

require(distrEx)
distrEx::TotalVarDist(Norm(mean = 0,sd = 1),Norm(mean = 1, sd = 1))

s <- 0.5
bin_test <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "binomial-test", threshold = 0.5)
hyper_test <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "hypergeometric-test", z = sum(rho <= mean(rho)), tv.seq = seq(0,1,0.001))
as_test <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "asymptotic-tv-search")
bin_test$tvhat
hyper_test$tvhat
as_test$tvhat
bin_test$tvhat / hyper_test$tvhat

# univariate classification (asymptotic wins)
set.seed(0)
n <- 5000
lambda <- 0.03
t <- sample(c(0,1), n, replace = TRUE)
x <- ifelse(t <= 0, ifelse(runif(n) <= lambda, rnorm(n, mean = 10), rnorm(n, mean = 0)), rnorm(n, mean = 0))
rho = dnorm(x, mean = 0) / ((2-lambda)*dnorm(x, mean = 0) + lambda*dnorm(x, mean = 10))
par(mfrow=c(2,1))
#plot(t, x, pch=19)
#plot(t, rho, pch=19)


require(distrEx)
distrEx::TotalVarDist(Norm(mean = 0, sd = 1),
                      UnivarMixingDistribution(Norm(mean = 0, sd = 1), Norm(mean = 10, sd = 1), mixCoeff = c(1-lambda, lambda)))


s <- 0.5
bin_test <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "binomial-test")
hyper_test <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "hypergeometric-test", z = sum(t==0))
as_test <- dWit(t = t, rho = rho, s = s, alpha = 0.05, estimator.type = "asymptotic-tv-search", verbose.plot = TRUE)
bin_test$tvhat
hyper_test$tvhat
as_test$tvhat
bin_test$tvhat / hyper_test$tvhat





