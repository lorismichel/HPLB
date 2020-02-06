# simulation continuous change point detection

set.seed(7)
# example 1: mean shift

# PLOT_SIMULATION_MIXTURE_1.png
par(mfrow=c(3,2))

n <- 10000

t.train <- runif(n, 0, 1)
x.train <- ifelse(t.train<=1/3, rnorm(n, 0), ifelse(t.train<=2/3, rnorm(n, 1), rnorm(n, 2)))
t.test <- runif(n, 0, 1)
x.test <- ifelse(t.test<=1/3, rnorm(n, 0), ifelse(t.test<=2/3, rnorm(n, 1), rnorm(n, 2)))

plot(t.train,x.train,pch=19,cex=0.5,xlab='Time', ylab='x', font.lab=1, font.main=1)
rf <- ranger::ranger(t~., data = data.frame(t=t.train, x=x.train))

# compute the mixture TV
require(dWit)

tvhat <- dWit(t = t.test, rho = predict(rf, data.frame(x=x.test))$predictions,
              s = seq(0.05, 0.95, 0.05), estimator.type = "asymptotic-tv-search")
plot(seq(0.05, 0.95, 0.05), tvhat$tvhat,type='b',pch=19, font.lab=1, font.main=1,ylab=expression(hat(lambda)[H]),xlab="Time")

# example 2: sd shift

n <- 10000

t.train <- runif(n, 0, 1)
x.train <- ifelse(t.train<=1/3, rnorm(n, sd=1), ifelse(t.train<=2/3, rnorm(n, sd=2), rnorm(n, sd=3)))
t.test <- runif(n, 0, 1)
x.test <- ifelse(t.test<=1/3, rnorm(n, sd=1), ifelse(t.test<=2/3, rnorm(n, sd=2), rnorm(n, sd=3)))

plot(t.train,x.train,pch=19,cex=0.5,xlab='Time', ylab='x', font.lab=1, font.main=1)
rf <- ranger::ranger(t~., data = data.frame(t=t.train, x=x.train))

# compute the mixture TV
require(dWit)

tvhat <- dWit(t = t.test, rho = predict(rf, data.frame(x=x.test))$predictions,
              s = seq(0.05, 0.95, 0.05), estimator.type = "asymptotic-tv-search")
plot(seq(0.05, 0.95, 0.05), tvhat$tvhat,type='b',pch=19,font.lab=1, font.main=1,ylab=expression(hat(lambda)[H]),xlab="Time")


# example 3: linear mean shift

n <- 10000

t.train <- runif(n, 0, 1)
x.train <- rnorm(n, 2*t.train, sd=1)
t.test <- runif(n, 0, 1)
x.test <- rnorm(n, 2*t.test, sd=1)

plot(t.train,x.train,pch=19,cex=0.5,xlab='Time', ylab='x',font.lab=1, font.main=1)
rf <- ranger::ranger(t~., data = data.frame(t=t.train, x=x.train))

# compute the mixture TV
require(dWit)

tvhat <- dWit(t = t.test, rho = predict(rf, data.frame(x=x.test))$predictions,
              s = seq(0.05, 0.95, 0.05), estimator.type = "asymptotic-tv-search")
plot(seq(0.05, 0.95, 0.05), tvhat$tvhat,type='b',pch=19,font.lab=1, font.main=1,ylab=expression(hat(lambda)[H]),xlab="Time")


# PLOT_SIMULATION_MIXTURE_2.png
par(mfrow=c(3,2))
# example 4: multivariate change
n <- 10000
d <- 2
x.train <- matrix(nrow=n, ncol=d)
x.test <- matrix(nrow=n, ncol=d)
t.train <- runif(n, 0, 1)
t.test <- runif(n, 0, 1)

for (i in 1:n) {
  x.train[i, ] <-  if (t.train[i]<=0.5) MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = (matrix(0.5, ncol=d,nrow=d) + diag(0.5, d))) else MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = (matrix(-0.5, ncol=d,nrow=d) + diag(1.5, d)))
  x.test[i, ] <- if (t.test[i]<=0.5) MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = (matrix(0.5, ncol=d,nrow=d) + diag(0.5, d))) else MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = (matrix(-0.5, ncol=d,nrow=d) + diag(1.5, d)))
}

# fits
rf_x1 <- ranger::ranger(t~., data = data.frame(t=t.train, x=x.train[,1]))
rf_x2 <- ranger::ranger(t~., data = data.frame(t=t.train, x=x.train[,2]))
rf_joint <- ranger::ranger(t~., data = data.frame(t=t.train, x=x.train))

# compute the mixture TV
require(dWit)

# dwit
tvhat_joint <- dWit(t = t.test, rho = predict(rf_joint, data.frame(x=x.test))$predictions,
              s = seq(0.05, 0.95, 0.05), estimator.type = "asymptotic-tv-search")
tvhat_x1 <- dWit(t = t.test, rho = predict(rf_x1, data.frame(x=x.test[,1]))$predictions,
              s = seq(0.05, 0.95, 0.05), estimator.type = "asymptotic-tv-search")
tvhat_x2 <- dWit(t = t.test, rho = predict(rf_x2, data.frame(x=x.test[,2]))$predictions,
              s = seq(0.05, 0.95, 0.05), estimator.type = "asymptotic-tv-search")

plot(t.train,x.train[,1],pch=19,cex=0.5,xlab='Time', ylab=expression(x[1]),font.lab=1, font.main=1)
plot(seq(0.05, 0.95, 0.05), tvhat_x1$tvhat,type='b',pch=19,font.lab=1, font.main=1,ylab=expression(hat(lambda)[H]),xlab="Time")
plot(t.train,x.train[,2],pch=19,cex=0.5,xlab='Time', ylab=expression(x[2]),font.lab=1, font.main=1)
plot(seq(0.05, 0.95, 0.05), tvhat_x2$tvhat,type='b',pch=19,font.lab=1, font.main=1,ylab=expression(hat(lambda)[H]),xlab="Time")
plot(x.train[,1],x.train[,2],pch=19,cex=0.5,xlab=expression(x[1]), ylab=expression(x[2]),col=c("darkblue","grey")[(t.train>0.5)+1],font.lab=1, font.main=1)
plot(seq(0.05, 0.95, 0.05), tvhat_joint$tvhat,type='b',pch=19,font.lab=1, font.main=1,ylab=expression(hat(lambda)[H]),xlab="Time")


