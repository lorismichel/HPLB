# Introductory example for section 1.2
introExample <- function(n = 1000, p = 0.01,   d = 10,
                         m.core = c(0, 0, rep(0, d)),
                         m.outlier = c(3, 3, rep(0, d))) {
  ids.outliers <- sample(c(0, 1), size = n, prob = c(1-p,p), replace = TRUE)
  x1 <- matrix(rnorm(n*(2+d)), ncol=2+d)
  x2 <- matrix(rnorm(n*(2+d)), ncol=2+d)
  x1 <- sweep(x1, MARGIN = 2,  m.core, FUN = "+")
  x2.outlier <- sweep(x2, MARGIN = 2,  m.outlier, FUN = "+")
  x2[which(ids.outliers==1), ] <- x2.outlier[which(ids.outliers==1), ]

  bayesRate <- function(x) ((1-p) * dnorm(x[1])*dnorm(x[2]) + p * dnorm(x[1], m.outlier[1]) * dnorm(x[2], m.outlier[2]))/(((1-p) * dnorm(x[1])*dnorm(x[2]) + p * dnorm(x[1], m.outlier[1]) * dnorm(x[2], m.outlier[2])) + (dnorm(x[1], m.core[1])*dnorm(x[2], m.core[2])))

  rho <- as.numeric(apply(rbind(x1,x2), 1, function(xx) bayesRate(xx)))
  return(list(x1=x1, x2=x2, t = c(rep(0, n), rep(1, n)), rho = rho, bayesRate = bayesRate))
}

computeStats <- function(t, rho, s) {
  test <- as.numeric(rho>s)
  return(c(sum(test[t==1])/sum(t==1), sum(1-test[t==0])/sum(t==0)))
}

## plots for paper with random forest fit

# repro
set.seed(1)

# training data
intro_training <- introExample(n = 10000, p = 0.01)
require(ranger)

# forest fit
rf <- ranger(formula = factor(y)~. , data.frame(y=intro_training$t, x1=rbind(intro_training$x1, intro_training$x2)), num.trees = 500, probability = TRUE)

# testing data
intro_testing <- introExample(n = 10000, p = 0.01)

# predictions
preds <- predict(rf, data.frame(y=intro_testing$t, x1=rbind(intro_testing$x1, intro_testing$x2)))$predictions

# permutation test
permPval <- function(tab,nsim=5000) {

  cs <- apply(tab,2,sum)
  rs <- apply(tab,1,sum)

  stat <- tab[1,1]/rs[1] + tab[2,2]/rs[2] - 1


  stat.vec <- rep(NA, nsim)
  for (i in 1:nsim) {
    s <- sample(c(rep(0,cs[1]),rep(1,cs[2])))
    ss <- table(c(rep(FALSE, rs[1]), rep(TRUE, rs[2])), s)
    stat.vec[i] <- ss[1,1]/rs[1] + ss[2,2]/rs[2] - 1
  }
  return((sum(stat.vec >= stat)+1)/(nsim+1))
}

par(mfrow=c(1,1))
par(mar=rep(4,4))

png("./Plots/PLOT_EXAMPLE_DATA.png", width = 1000, height = 1000)
plot(intro_testing$x1, pch=19, col="darkgrey",cex=0.5,font.lab=1,xlab=expression(X[1]),ylab=expression(X[1]),xlim=c(-5,5),ylim=c(-5,5))
points(intro_testing$x2, pch=19, col="darkblue",cex=0.5)
dev.off()

# permutations p-values
permPval(table(preds[,2]>0.5, intro_testing$t), nsim = 1000)
permPval(table(preds[,2]>0.7, intro_testing$t), nsim = 1000)


# Actual estimated TV values

dWit(t = intro_testing$t, rho = preds[,2],
     s = 0.5,
     estimator.type = "binomial-test")


dWit(t = intro_testing$t, rho = preds[,2],
     s = 0.5,
     estimator.type = "asymptotic-tv-search")










# case where 0.7 is winning
# set.seed(1)
# perm1 <- sapply(1:1000, function(i) {ex <- introExample(n = 10000, p = 0.01); computeStats(t = ex$t, rho = ex$rho, s = 0.5)})
# set.seed(1)
# perm2 <- sapply(1:1000, function(i) {ex <- introExample(n = 10000, p = 0.01); computeStats(t = ex$t, rho = ex$rho, s = 0.7)})
#
#
# set.seed(1)
# perm3 <- sapply(1:1000, function(i) {ex <- introExample(n = 1000, p = 0.1); computeStats(t = ex$t, rho = ex$rho, s = 0.5)})
# set.seed(1)
# perm4 <- sapply(1:1000, function(i) {ex <- introExample(n = 1000, p = 0.1); computeStats(t = ex$t, rho = ex$rho, s = 0.7)})
#
#
# # range
# range(perm1[1,]+perm1[2,]-1)
# range(perm2[1,]+perm2[2,]-1)
# range(perm3[1,]+perm3[2,]-1)
# range(perm4[1,]+perm4[2,]-1)
#
# # plots
# par(mfrow=c(2,2))
#
# hist(perm1[1,]+perm1[2,]-1, xlim=c(0.0001,0.1),main="statistic at 0.5",
#      xlab="Statistic", ylab="Density",font.main=1,font.lab=1, breaks = 30)
# abline(v = quantile(perm1[1,]+perm1[2,]-1,0.05),col="red")
# hist(perm3[1,]+perm3[2,]-1, xlim=c(0.0001,0.2),main="statistic at 0.5",
#      xlab="Statistic", ylab="Density",font.main=1,font.lab=1, breaks = 30)
# abline(v = quantile(perm3[1,]+perm3[2,]-1,0.05),col="red")
# hist(perm2[1,]+perm2[2,]-1, xlim=c(0.0001,0.1),main="statistic at 0.7",
#      xlab="Statistic", ylab="Density",font.main=1,font.lab=1, breaks = 30)
# abline(v = quantile(perm2[1,]+perm2[2,]-1,0.05),col="red")
# hist(perm4[1,]+perm4[2,]-1, xlim=c(0.0001,0.2),main="statistic at 0.7",
#      xlab="Statistic", ylab="Density",font.main=1,font.lab=1, breaks = 30)
# abline(v = quantile(perm4[1,]+perm4[2,]-1,0.05),col="red")
#
# par(pty="s")
# par(mfrow=c(2,1))
# par(mar = c(1,1,0,1))
# set.seed(1)
# e1 <- introExample(n = 10000, p = 0.01)
# plot(e1$x1[,1], e1$x1[,2], asp = 2, xlim=c(-7,7),ylim=c(-7,7),pch=19,font.main=1,font.lab=1,xlab=expression(x[1]),ylab=expression(x[2]),cex=0.5,col="black")
# points(e1$x2[,1], e1$x2[,2], pch=19,col="darkblue",cex=0.5)
#
# par(mar = c(1,1,0,1))
# set.seed(1)
# e1 <- introExample(n = 10000, p = 0.1)
# plot(e1$x1[,1], e1$x1[,2], xlim=c(-7,7),ylim=c(-7,7),pch=19,font.main=1,font.lab=1,xlab=expression(x[1]),ylab=expression(x[2]),cex=0.5,col="black")
# points(e1$x2[,1], e1$x2[,2], pch=19,col="darkblue",cex=0.5)
