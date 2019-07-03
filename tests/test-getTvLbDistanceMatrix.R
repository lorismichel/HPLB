# test of the confusion
set.seed(0)
n1 <- n2 <- n3 <- n4 <-  1000
n <- n1+n2+n3+n4
x1 <- rnorm(n1,-2)
x2 <- rnorm(n2,0)
x3 <- rnorm(n3,2)
x4 <- rnorm(n4,2)
x <- c(x1,x2,x3,x4)
y <- rep(0:3, c(n1,n2,n3,n4))
d <- data.frame(x = x, y = y)


# true tvs
require(distrEx)

mat.vec <- c(-2,0,2,2)
tv.dist <- matrix(NA,4,4)
for (i in 1:4) {
  for (j in 1:4) {
    tv.dist[i,j] <- TotalVarDist(Norm(mat.vec[i],1),Norm(mat.vec[j],1))
  }
}

# constrcution of ordering array
ar <- array(dim = c(4, 4, n))


# fixed function
for (i in 1:4) {
  for (j in 1:4) {
    ar[i,j,] <- x
  }
}


# multiclass predictions
ind.train <- sample(1:nrow(d), size = nrow(d)/2, replace = FALSE)
rf <- ranger::ranger(factor(y)~., probability = TRUE, data = d[ind.train, ])
preds <- predict(rf, data = d[-ind.train, ])$predictions

for (i in 1:4) {
  for (j in 1:4) {
    ar[i,j,] <- preds[,j] - preds[,i]
  }
}


# estimate dist
dist.mat <- getTvLbDistanceMatrix(labels = d[-ind.train, ]$y, ordering.array = ar)
dist.mat
tv.dist

# compare true mat with estim
tv.dist - dist.mat
par(mfrow=c(2,1))
image(tv.dist)
image(dist.mat)

