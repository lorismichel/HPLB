# libs
require(dWit)
require(distrEx)

# dataset
d<- getDataSet(dataset.name = "UnivariateNormalShift+Contamination", lambda = 0.01, mean = 1, sd = 1, m.train = 2000, n.train = 2000, m.test = 2000, n.test = 2000)

# fit a forest
rf <- ranger::ranger(factor(t)~., d$train, probability = TRUE)
rho <- predict(rf, d$test)$predictions[,"1"]

# get tv
dWit(t = d$test$t, rho = rho, s = 0.5)
dWit(t = d$test$t, rho = rho, s = 0.5, estimator.type = "binomial")

# apply our strategy

# 1)  fit a first forest for the reordering
rf.reorder <- ranger::ranger(factor(t)~., d$train, probability =  TRUE)
rho.reorder <- predict(rf, d$test)$predictions[,"1"]


par(mfrow=c(4,2))
# no
x.train <- d$train$x
t.train <- permLabels(d$train$t, preds = rf.reorder$predictions[,"1"], prob = 0, u = 0.9, d = 0)
x.test <- d$test$x
t.test <- permLabels(d$test$t, preds = rho.reorder, prob = 0, u = 0.9, d = 0)
plot(density(x.train[t.train==1]),col="red")
lines(density(x.train[t.train==0]))
rf <- ranger::ranger(factor(t)~., data.frame(t=t.train,x=x.train), probability = TRUE)
rho <- predict(rf, data.frame(x=x.test))$predictions[,"1"]


dWit(t = t.test, rho = rho, s = 0.5, verbose.plot = TRUE)
dWit(t = t.test, rho = rho, s = 0.5, estimator.type = "binomial")

# mit
x.train <- d$train$x
t.train <- permLabels(d$train$t, preds = rf.reorder$predictions[,"1"], prob = 1/2, u = 0.9, d = 0)
x.test <- d$test$x
t.test <- permLabels(d$test$t, preds = rho.reorder, prob = 1/2, u = 0.9, d = 0)
plot(density(x.train[t.train==1]),col="red")
lines(density(x.train[t.train==0]))
rf <- ranger::ranger(factor(t)~., data.frame(t=t.train,x=x.train), probability = TRUE)
rho <- predict(rf, data.frame(x=x.test))$predictions[,"1"]


dWit(t = t.test, rho = rho, s = 0.5, verbose.plot = TRUE)
dWit(t = t.test, rho = rho, s = 0.5, estimator.type = "binomial")


# hard
x.train <- d$train$x
t.train <- permLabels(d$train$t, preds = rf.reorder$predictions[,"1"], prob = 1/2, u = 0.95, d = 0)
x.test <- d$test$x
t.test <- permLabels(d$test$t, preds = rho.reorder, prob = 1/2, u = 0.95, d = 0)
plot(density(x.train[t.train==1]),col="red")
lines(density(x.train[t.train==0]))
rf <- ranger::ranger(factor(t)~., data.frame(t=t.train,x=x.train), probability = TRUE)
rho <- predict(rf, data.frame(x=x.test))$predictions[,"1"]


dWit(t = t.test, rho = rho, s = 0.5, verbose.plot = TRUE)
dWit(t = t.test, rho = rho, s = 0.5, estimator.type = "binomial")


# very hard
x.train <- d$train$x
t.train <- permLabels(d$train$t, preds = rf.reorder$predictions[,"1"], prob = 1/2, u = 0.99, d = 0)
x.test <- d$test$x
t.test <- permLabels(d$test$t, preds = rho.reorder, prob = 1/2, u = 0.99, d = 0)
plot(density(x.train[t.train==1]),col="red")
lines(density(x.train[t.train==0]))
rf <- ranger::ranger(factor(t)~., data.frame(t=t.train,x=x.train), probability = TRUE)
rho <- predict(rf, data.frame(x=x.test))$predictions[,"1"]


dWit(t = t.test, rho = rho, s = 0.5, verbose.plot = TRUE)
dWit(t = t.test, rho = rho, s = 0.5, estimator.type = "binomial")


