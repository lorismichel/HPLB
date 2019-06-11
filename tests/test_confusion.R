# example confusion with univariate data
n <- 1000
x <- matrix(c(rnorm(n), rnorm(n, mean = 2), rnorm(n, sd = 3), rnorm(n,mean = -4)),ncol=1)
labels <- factor(rep(1:4, each=n))
list.dist <- list(Norm(0,1),Norm(2,1),Norm(0,3),Norm(-4,1))
true.tv.mat <- matrix(0, ncol=4,nrow=4)
for (i in 1:4) {
  for (j in 1:4) {
    if (i != j) {
      true.tv.mat[i,j] <- TotalVarDist(list.dist[[i]], list.dist[[j]])
    }
  }
}

set.seed(0)
true.tv.mat
pairWiseWitConfusion(x = x, labels = labels)
