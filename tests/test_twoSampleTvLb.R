# testing the two-sample TV lower-bound
set.seed(0)
x1 <- rnorm(1000)
x2 <- rnorm(1000, mean = -1)
d <- data.frame(cbind(x=c(x1,x2),y=rep(0:1, each = 1000)))
d <- d[sample(1:nrow(d)),]
d$y <- factor(d$y)
# half of the data for training the classifier
ind.train <- sample(1:nrow(d), size = nrow(d)/2)
# single tree
rp <- rpart::rpart(y~., data = d[ind.train,])

# building the predictions and labels
preds <- predict(rp, newdata = d[-ind.train,])[,"1"]
labels <- d[-ind.train,]$y

# computing the lower-bound
twoSampleTvLb(labels = labels,
              preds = preds)

# verify that value
require(distrEx)
TotalVarDist(e1 = Norm(mean = 0,1), e2 = Norm(mean = -1,1))

