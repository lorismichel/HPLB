
# source
source("./climate-preprocessing.R")

# libs
require(dWit)

# prepro
d <- climatePrePro()

# analysis
dat <- data.frame(class=d$tenyears.ind,
                  x = cbind(t(d$air),t(d$mslp),t(d$prate),t(d$shum)))

# splits train-test
ind.train <- unlist(lapply(lapply(unique(dat$class), function(i) which(dat$class==i)), function(l) l[l > quantile(l,0.25) & l <= quantile(l,0.75)]))
ind.test <- unlist(lapply(lapply(unique(dat$class), function(i) which(dat$class==i)), function(l) l[l <= quantile(l,0.25) | l > quantile(l,0.75)]))

dat.train <- dat[ind.train,]
dat.test <- dat[ind.test,]

# fit a forest
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


# build the rho function
for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)
tv.mat
