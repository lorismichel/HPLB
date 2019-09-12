# climatic data (real analysis)

# libs
library(sf)
library(ncdf4)
library(raster)
library(rasterVis)
library(RColorBrewer)



# PATH of data
ncpath <- "~/Downloads/reanalysis_jeff_loris/"


# loading the data (4 signals)
air_raster <- brick(paste(ncpath, "NCEP2.air.2m.day", ".nc", sep=""), varname="air")
mslp_raster <- brick(paste(ncpath, "NCEP2.mslp.2m.day", ".nc", sep=""), varname="mslp")
prate_raster <- brick(paste(ncpath, "NCEP2.prate.2m.day", ".nc", sep=""), varname="prate")
shum_raster <- brick(paste(ncpath, "NCEP2.shum.2m.day", ".nc", sep=""), varname="shum")



########### preprocessing

# transform to date
toDate <- function(x) {
  return(sub(sub(substr(x, start = 2, stop = 11), pattern = "\\.", replacement = "-"), pattern = "\\.", replacement = "-"))
}

# times
t_air <- as.Date(sapply(colnames(air_raster[1,1]), toDate))
t_mslp <- as.Date(sapply(colnames(mslp_raster[1,1]), toDate))
t_prate <- as.Date(sapply(colnames(prate_raster[1,1]), toDate))
t_shum <- as.Date(sapply(colnames(shum_raster[1,1]), toDate))

# subset to smaller time step
air <- getValues(air_raster)[,1:14641]
mslp <- getValues(mslp_raster)[,1:14641]
prate <- getValues(prate_raster)[,1:14641]
shum <- getValues(shum_raster)[,1:14641]


fake.series <- matrix(nrow=nrow(air),ncol=ncol(air))
for (i in 1:nrow(fake.series)) {
  fake.series[i,] <- rnorm(n=14641, mean = (10^-4)*1:14641)
}

########### analysis
# here we intend to start with a certain block size and decrease it more and more
# this can be done in regression or classification task

# build indicator of 10 years in a row
time <- t_air
tenyears.ind <- cut(time, breaks = as.Date(c("1979-01-01",
                                             "1989-01-01",
                                             "1999-01-01",
                                             "2009-01-01",
                                             "2019-02-01")), labels=FALSE)





# filtering of time


## function that makes the block size analysis
decendingBlocksReg <- function(signals = c("air"), max.p = 9) {

  if (signals == "air") {
    x <- t(air)
  } else if (signals == "mslp") {
    x <- t(mslp)
  } else if (signals == "prate") {
    x <- t(prate)
  } else if (signals == "shum") {
    x <- t(shum)
  } else {
    x <- cbind(t(air),t(mslp),t(prate),t(shum))
  }

  mean.local.tv.search <- c()
  mean.local.tv.bin <- c()
  for (p in 1:max.p) {
    ind.time <- cut(as.numeric(time), breaks = 2^p, labels = FALSE)
    dat <- data.frame(time = as.numeric(time), class=ind.time, x = x)
    ind.train <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[c(1:(length(l)/2))]))
    ind.test <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[-c(1:(length(l)/2))]))
    dat.train <- dat[ind.train,]
    dat.test <- dat[ind.test,]

    # regression forest
    regRF <- ranger::ranger(time~.-class, data = dat.train, seed = 23)


    ordering.array <- array(dim=c(length(unique(ind.time)),length(unique(ind.time)),nrow(dat.test)))
    labels <- dat.test$class-1

    # compute ordering
    preds <- predict(regRF, data = dat.test)$predictions
    ordering.array <- preds


    # get the tv mat
    source("./R/dWit.R")
    source("./R/invertBinMeanTest.R")
    source("./R/getTvLbDistanceMatrix.R")
    tv.mat.search <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array, estimator.type = "tv-search")
    tv.mat.bin <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array, estimator.type = "binomial")
    values.search <- c()
    values.bin <- c()
    for (i in 2:nrow(tv.mat.search)) {
      values.search <- c(values.search, tv.mat.search[i,i-1])
      values.bin <- c(values.bin, tv.mat.bin[i,i-1])
    }
    mean.local.tv.search <- c(mean.local.tv.search, mean(values.search))
    mean.local.tv.bin <- c(mean.local.tv.bin, mean(values.bin))
    print(p)
  }

  return(list(mean.local.tv.search = mean.local.tv.search, mean.local.tv.bin = mean.local.tv.bin))
}


## function that makes the block size analysis
decendingBlocksBinClass <- function(signals = c("air"),
                                    max.p = 9, filter = 1:length(time),
                                    odd.even.split = FALSE) {

  if (signals == "air") {
    x <- t(air[,filter])
  } else if (signals == "mslp") {
    x <- t(mslp[,filter])
  } else if (signals == "prate") {
    x <- t(prate[,filter])
  } else if (signals == "shum") {
    x <- t(shum[,filter])
  } else {
    x <- cbind(t(air[,filter]),t(mslp[,filter]),t(prate[,filter]),t(shum[,filter]))
  }

  mean.local.tv.search <- list()
  mean.local.tv.bin <- list()
  for (p in 1:max.p) {
    ind.time <- cut(as.numeric(time[filter]), breaks = 2^p, labels = FALSE)
    dat <- data.frame(time = as.numeric(time[filter]), class=ind.time, x = x)

    if (!odd.even.split) {
      ind.train <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[c(1:(length(l)/2))]))
      ind.test <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[-c(1:(length(l)/2))]))
    } else {
      ind.train <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[c(1:length(l))%%2==0]))
      ind.test <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[c(1:length(l))%%2==1]))
    }

    dat.train <- dat[ind.train,]
    dat.test <- dat[ind.test,]


    tv.bin <- c()
    tv.search <- c()
    for (i in 1:c(length(unique(ind.time))-1)) {

      class.ind.train <- as.numeric(dat.train$class==(i+1))[dat.train$class%in%c(i,i+1)]
      class.ind.test <- as.numeric(dat.test$class==(i+1))[dat.test$class%in%c(i,i+1)]
      bcRF <- ranger::ranger(class.ind.train~.-class-time, data = dat.train[dat.train$class%in%c(i,i+1),], seed = 23, probability = TRUE)
      preds <- predict(bcRF, dat.test[dat.test$class%in%c(i,i+1),])$predictions
      tv.bin <- c(tv.bin , dWit(t = class.ind.test, rho = preds[,2], estimator.type = "binomial")$tvhat)
      tv.search <- c(tv.search ,dWit(t = class.ind.test, rho = preds[,2], estimator.type = "tv-search")$tvhat)
      print(i)
    }
    print(paste("p done."))
    mean.local.tv.bin[[p]] <- tv.bin
    mean.local.tv.search[[p]] <- tv.search

  }
  return(list(mean.local.tv.search = mean.local.tv.search,
              mean.local.tv.bin = mean.local.tv.bin))
}


## function that makes the block size analysis
decendingBlocksBinClassMix <- function(signals = c("air"),
                                    max.p = 9, min.p = 1, filter = 1:length(time),
                                    prob = 0, locations = c(1:nrow(air)),
                                    use.lda = FALSE) {
  if (signals == "fake") {
    x <- t(fake.series[locations,filter])
  } else if (signals == "air") {
    x <- t(air[locations,filter])
  } else if (signals == "mslp") {
    x <- t(mslp[locations,filter])
  } else if (signals == "prate") {
    x <- t(prate[locations,filter])
  } else if (signals == "shum") {
    x <- t(shum[locations,filter])
  } else {
    x <- cbind(t(air[locations,filter]),t(mslp[locations,filter]),t(prate[locations,filter]),t(shum[locations,filter]))
  }

  mean.local.tv.search <- list()
  mean.local.tv.bin <- list()
  mean.local.tv.empirical <- list()
  for (p in min.p:max.p) {
    ind.time <- cut(as.numeric(time[filter]), breaks = 2^p, labels = FALSE)
    dat <- data.frame(time = as.numeric(time[filter]), class=ind.time, x = x)

    ind.train <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[l > quantile(l,0.25) & l <= quantile(l,0.75)]))
    ind.test <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[l <= quantile(l,0.25) | l > quantile(l,0.75)]))

    dat.train <- dat[ind.train,]
    dat.test <- dat[ind.test,]


    tv.bin <- c()
    tv.search <- c()
    tv.empirical <- c()

    if (p >= 5) {
      vals <- unique(table(ind.time[ind.test]))

      e11 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[1], n = vals[1])

      if (length(vals)>1) {
        e12 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[1], n = vals[2])
        e21 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[2], n = vals[1])
        e22 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[2], n = vals[2])
      }
    }


    for (i in 1:c(length(unique(ind.time))-1)) {

      class.ind.train <- as.numeric(dat.train$class==(i+1))[dat.train$class%in%c(i,i+1)]
      class.ind.test <- as.numeric(dat.test$class==(i+1))[dat.test$class%in%c(i,i+1)]

      class.ind.train <- ifelse(sample(c(0,1),
                                       prob = c(1-prob, prob),
                                       replace = TRUE,
                                       size = length(class.ind.train))==1, 1-class.ind.train, class.ind.train)
      class.ind.test <- ifelse(sample(c(0,1),
                                       prob = c(1-prob, prob),
                                       replace = TRUE,
                                       size = length(class.ind.test))==1, 1-class.ind.test, class.ind.test)



      if (use.lda) {
        require(MASS)
        bcLDA <- lda(class.ind.train~.-class-time, data = dat.train[dat.train$class%in%c(i,i+1),])
        preds <- predict(bcLDA, dat.test[dat.test$class%in%c(i,i+1),])$posterior[,"1"]
      } else {
        bcRF <- ranger::ranger(class.ind.train~.-class-time, data = dat.train[dat.train$class%in%c(i,i+1),], seed = 23, probability = TRUE)
        preds <- predict(bcRF, dat.test[dat.test$class%in%c(i,i+1),])$predictions
      }

      tv.bin <- c(tv.bin , dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "binomial")$tvhat)
      tv.search <- c(tv.search ,dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "asymptotic-tv-search")$tvhat)

      if (p >= 5) {
        b <- if (sum(class.ind.test) == vals[1]) {
          if (sum(1-class.ind.test) == vals[1]) {
            e11
          } else {
            e21
          }
        } else {
          if (sum(1-class.ind.test) == vals[1]) {
            e12
          } else {
            e22
          }
        }

        tv.empirical <- c(tv.empirical, dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "custom-tv-search",
                                             custom.bounding.seq = b)$tvhat)
      } else {
        tv.empirical <- c(tv.empirical, NA)
      }
      print(i)
    }
    print(paste("p done."))
    mean.local.tv.bin[[p]] <- tv.bin
    mean.local.tv.search[[p]] <- tv.search
    mean.local.tv.empirical[[p]] <- tv.empirical

  }
  return(list(mean.local.tv.search = mean.local.tv.search,
              mean.local.tv.bin = mean.local.tv.bin,
              mean.local.tv.empirical = mean.local.tv.empirical))
}




## function that makes the block size analysis
runClimateAnalysis <- function(signals = c("air"),
                                       max.p = 9, min.p = 1, filter = 1:length(time),
                                       prob = 0, locations = c(1:nrow(air)),
                                       use.lda = FALSE,
                                       spatial.centering = TRUE) {
  if (signals == "fake") {
    x <- t(fake.series[locations,filter])
  } else if (signals == "air") {
    x <- t(air[locations,filter])
  } else if (signals == "mslp") {
    x <- t(mslp[locations,filter])
  } else if (signals == "prate") {
    x <- t(prate[locations,filter])
  } else if (signals == "shum") {
    x <- t(shum[locations,filter])
  } else {
    x <- cbind(t(air[locations,filter]),t(mslp[locations,filter]),t(prate[locations,filter]),t(shum[locations,filter]))
  }

  if (spatial.centering) {
    x <- apply(x, 2, function(xx) xx-mean(xx,na.rm=T))
  }

  mean.local.tv.search <- list()
  mean.local.tv.bin <- list()
  mean.local.tv.empirical <- list()
  for (p in min.p:max.p) {
    ind.time <- cut(as.numeric(time[filter]), breaks = 2^p, labels = FALSE)
    dat <- data.frame(time = as.numeric(time[filter]), class=ind.time, x = x)

    ind.train <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[l > quantile(l,0.25) & l <= quantile(l,0.75)]))
    ind.test <- unlist(lapply(lapply(unique(ind.time), function(i) which(ind.time==i)), function(l) l[l <= quantile(l,0.25) | l > quantile(l,0.75)]))

    dat.train <- dat[ind.train,]
    dat.test <- dat[ind.test,]


    tv.bin <- c()
    tv.search <- c()
    tv.empirical <- c()

    if (p >= 5) {
      vals <- unique(table(ind.time[ind.test]))

      e11 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[1], n = vals[1])

      if (length(vals)>1) {
        e12 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[1], n = vals[2])
        e21 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[2], n = vals[1])
        e22 <- empiricalBF(tv.seq = seq(from = 0, to = 1, by = 0.01), nrep = 1000, m = vals[2], n = vals[2])
      }
    }


    for (i in 1:c(length(unique(ind.time))-1)) {

      class.ind.train <- as.numeric(dat.train$class==(i+1))[dat.train$class%in%c(i,i+1)]
      class.ind.test <- as.numeric(dat.test$class==(i+1))[dat.test$class%in%c(i,i+1)]

      class.ind.train <- ifelse(sample(c(0,1),
                                       prob = c(1-prob, prob),
                                       replace = TRUE,
                                       size = length(class.ind.train))==1, 1-class.ind.train, class.ind.train)
      class.ind.test <- ifelse(sample(c(0,1),
                                      prob = c(1-prob, prob),
                                      replace = TRUE,
                                      size = length(class.ind.test))==1, 1-class.ind.test, class.ind.test)



      if (use.lda) {
        require(MASS)
        bcLDA <- lda(class.ind.train~.-class-time, data = dat.train[dat.train$class%in%c(i,i+1),])
        preds <- predict(bcLDA, dat.test[dat.test$class%in%c(i,i+1),])$posterior[,"1"]
      } else {
        bcRF <- ranger::ranger(class.ind.train~.-class-time, data = dat.train[dat.train$class%in%c(i,i+1),], seed = 23, probability = TRUE)
        preds <- predict(bcRF, dat.test[dat.test$class%in%c(i,i+1),])$predictions
      }

      tv.bin <- c(tv.bin , dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "binomial")$tvhat)
      tv.search <- c(tv.search ,dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "asymptotic-tv-search")$tvhat)

      if (p >= 5) {
        b <- if (sum(class.ind.test) == vals[1]) {
          if (sum(1-class.ind.test) == vals[1]) {
            e11
          } else {
            e21
          }
        } else {
          if (sum(1-class.ind.test) == vals[1]) {
            e12
          } else {
            e22
          }
        }

        tv.empirical <- c(tv.empirical, dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "custom-tv-search",
                                             custom.bounding.seq = b)$tvhat)
      } else {
        tv.empirical <- c(tv.empirical, NA)
      }
      print(i)
    }
    print(paste("p done."))
    mean.local.tv.bin[[p]] <- tv.bin
    mean.local.tv.search[[p]] <- tv.search
    mean.local.tv.empirical[[p]] <- tv.empirical

  }
  return(list(mean.local.tv.search = mean.local.tv.search,
              mean.local.tv.bin = mean.local.tv.bin,
              mean.local.tv.empirical = mean.local.tv.empirical))
}


## runs

# with binary class
max.p <- 9
res2 <- decendingBlocksBinClass(signals = "shum", max.p = max.p, odd.even.split = TRUE)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

# checkpoints
abline(v=log(10),col="blue",lty=2)
abline(v=log(5),col="blue",lty=2)
abline(v=log(1),col="blue",lty=2)
abline(v=log(1/2),col="blue",lty=2)
abline(v=log(1/4),col="blue",lty=2)
abline(v=log(1/12),col="blue",lty=2)
length(time)/2^{1:max.p}

# with binary class just Dec, Jan, Feb, Match
max.p <- 5
res3 <- decendingBlocksBinClass(signals = c("prate"), max.p = max.p, filter = months(time)%in%c("February"))
plot(log(length(time)/2^{1:max.p}/365),lapply(res3$mean.local.tv.search,quantile,0.95) , ylim=c(0,1), type="l", col="black",pch=19,lty=2,xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res3$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res3$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res3$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res3$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res3$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


# analysis with mixing

par(mfrow=c(3,1))
for (prob in c(0,0.1,0.3)) {
max.p <- 5
filter <- months(time)%in%c("February")
res3 <- decendingBlocksBinClassMix(signals = c("mslp"),
                                   max.p = max.p,
                                   filter = filter,
                                   prob = prob,
                                   odd.even.split = TRUE)
plot(log(length(time[filter])/2^{1:max.p}/365),lapply(res3$mean.local.tv.search,quantile,0.95) , ylim=c(0,1), type="l", col="black",pch=19,lty=2,xlab="block size in log(years)",ylab="TV")
lines(log(length(time[filter])/2^{1:max.p}/365),lapply(res3$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time[filter])/2^{1:max.p}/365),lapply(res3$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time[filter])/2^{1:max.p}/365),lapply(res3$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time[filter])/2^{1:max.p}/365),lapply(res3$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time[filter])/2^{1:max.p}/365),lapply(res3$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
}
require(distrEx)
TotalVarDist(Norm(mean=(10^-4)*5), Norm(mean=(10^-4)*))
TotalVarDist(
UnivarMixingDistribution(lapply(1:30, FUN = function(i) Norm(mean=10^{-4}*i))),
UnivarMixingDistribution(lapply(31:60, FUN = function(i) Norm(mean=10^{-4}*i))))
block.size <- 30
TotalVarDist(do.call(what = "UnivarMixingDistribution", lapply(1:block.size, FUN = function(i) Norm(mean=10^{-4}*i))),
             do.call(what = "UnivarMixingDistribution", lapply((block.size+1):(2*block.size), FUN = function(i) Norm(mean=10^{-4}*i))))

# rough viz analysis
plot(shum[1,months(time)%in%c("February")],pch=19,type="l")
plot(shum[1,months(time)%in%c("February")][1:1000],pch=19,type="l")
plot(shum[1,months(time)%in%c("February")][1:500],pch=19,type="l")
plot(shum[1,months(time)%in%c("February")][1:100],pch=19,type="l")
plot(shum[1,months(time)%in%c("February")][1:30],pch=19,type="l")


plot(mslp[1,],pch=19,type="l")
plot(mslp[1,1:(365*10)],pch=19,type="l")
plot(mslp[1,1:(365*5)],pch=19,type="l")
plot(mslp[1,1:(365*3)],pch=19,type="l")
plot(mslp[1,1:(365*1)],pch=19,type="l")
plot(mslp[1,1:(365*0.5)],pch=19,type="l")




### ana
par(mfrow=c(1,1))
max.p <- 9
res2 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = TRUE, prob = 0)


plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,0.1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

block.size.vec <- round(length(time)/2^{1:max.p})
tv.vec <- c()
for (block.size in block.size.vec[-c(1:4)]) {
tv.vec <- c(tv.vec, TotalVarDist(do.call(what = "UnivarMixingDistribution", lapply(1:block.size, FUN = function(i) Norm(mean=10^{-4}*i))),
             do.call(what = "UnivarMixingDistribution", lapply((block.size+1):(2*block.size), FUN = function(i) Norm(mean=10^{-4}*i)))))
}
lines(log(length(time)/2^{1:max.p}/365), c(rep(-1,4), tv.vec),type="b",col="blue",pch=19)


# analysis of shum:
max.p <- 9
res_shum_odd_0 <- decendingBlocksBinClassMix(signals = "shum", max.p = max.p, odd.even.split = TRUE, prob = 0)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

max.p <- 9
res_shum_odd_01 <- decendingBlocksBinClassMix(signals = "shum", max.p = max.p, odd.even.split = TRUE, prob = 0.1)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


max.p <- 9
res_shum_odd_03 <- decendingBlocksBinClassMix(signals = "shum", max.p = max.p, odd.even.split = TRUE, prob = 0.3)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

save(res_shum_odd_0, res_shum_odd_01, res_shum_odd_03, file = "~/phdWork/dWit/Data/climate_res_shum_odd.Rdata")


max.p <- 9
res_shum_block_0 <- decendingBlocksBinClassMix(signals = "shum", max.p = max.p, odd.even.split = FALSE, prob = 0)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

max.p <- 9
res_shum_block_01 <- decendingBlocksBinClassMix(signals = "shum", max.p = max.p, odd.even.split = FALSE, prob = 0.1)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


max.p <- 9
res_shum_block_03 <- decendingBlocksBinClassMix(signals = "shum", max.p = max.p, odd.even.split = FALSE, prob = 0.3)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

save(res_shum_block_0, res_shum_block_01, res_shum_block_03, file = "~/phdWork/dWit/Data/climate_res_shum_block.Rdata")



# analysis of air:
max.p <- 9
res_air_odd_0 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = TRUE, prob = 0)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

max.p <- 9
res_air_odd_01 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = TRUE, prob = 0.1)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


max.p <- 9
res_air_odd_03 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = TRUE, prob = 0.3)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

save(res_air_odd_0, res_air_odd_01, res_air_odd_03, file = "~/phdWork/dWit/Data/climate_res_air_odd.Rdata")


max.p <- 9
res_air_block_0 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = FALSE, prob = 0)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

max.p <- 9
res_air_block_01 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = FALSE, prob = 0.1)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


max.p <- 9
res_air_block_03 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = FALSE, prob = 0.3)
plot(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res2$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

save(res_air_block_0, res_air_block_01, res_air_block_03, file = "~/phdWork/dWit/Data/climate_res_air_block.Rdata")


# by months
max.p <- 5
res_air_odd_Feb_0 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = TRUE, prob = 0, filter = months(time)=="February")

max.p <- 5
res_air_odd_Feb_01 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = TRUE, prob = 0.1, filter = months(time)=="February")

max.p <- 5
res_air_odd_Feb_03 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = TRUE, prob = 0.3, filter = months(time)=="February")

save(res_air_odd_Feb_0, res_air_odd_Feb_01, res_air_odd_Feb_03, file = "~/phdWork/dWit/Data/climate_res_air_odd_Feb.Rdata")


max.p <- 5
res_air_block_Feb_0 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = FALSE, prob = 0, filter = months(time)=="February")

max.p <- 5
res_air_block_Feb_01 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = FALSE, prob = 0.1, filter = months(time)=="February")

max.p <- 5
res_air_block_Feb_03 <- decendingBlocksBinClassMix(signals = "air", max.p = max.p, odd.even.split = FALSE, prob = 0.3, filter = months(time)=="February")

save(res_air_block_Feb_0, res_air_block_Feb_01, res_air_block_Feb_03, file = "~/phdWork/dWit/Data/climate_res_air_block_Feb.Rdata")



# analysis of prate:
max.p <- 9
res_prate_odd_0 <- decendingBlocksBinClassMix(signals = "prate", max.p = max.p, odd.even.split = TRUE, prob = 0)

max.p <- 9
res_prate_odd_01 <- decendingBlocksBinClassMix(signals = "prate", max.p = max.p, odd.even.split = TRUE, prob = 0.1)

max.p <- 9
res_prate_odd_03 <- decendingBlocksBinClassMix(signals = "prate", max.p = max.p, odd.even.split = TRUE, prob = 0.3)

save(res_prate_odd_0, res_prate_odd_01, res_prate_odd_03, file = "~/phdWork/dWit/Data/climate_res_prate_odd.Rdata")


max.p <- 9
res_prate_block_0 <- decendingBlocksBinClassMix(signals = "prate", max.p = max.p, odd.even.split = FALSE, prob = 0)

max.p <- 9
res_prate_block_01 <- decendingBlocksBinClassMix(signals = "prate", max.p = max.p, odd.even.split = FALSE, prob = 0.1)

max.p <- 9
res_prate_block_03 <- decendingBlocksBinClassMix(signals = "prate", max.p = max.p, odd.even.split = FALSE, prob = 0.3)

save(res_prate_block_0, res_prate_block_01, res_prate_block_03, file = "~/phdWork/dWit/Data/climate_res_prate_block.Rdata")


# analysis of mslp:
max.p <- 9
res_mslp_odd_0 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 7, prob = 0)


max.p <- 9
res_mslp_odd_0_loc_1 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 7, prob = 0,
                                                   locations = 1)


max.p <- 9
res_mslp_odd_01 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = TRUE, prob = 0.1)

max.p <- 9
res_mslp_odd_03 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = TRUE, prob = 0.3)

save(res_mslp_odd_0, res_mslp_odd_01, res_mslp_odd_03, file = "~/phdWork/dWit/Data/climate_res_mslp_odd.Rdata")


max.p <- 9
res_mslp_block_0 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = FALSE, prob = 0)

max.p <- 9
res_mslp_block_01 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = FALSE, prob = 0.1)

max.p <- 9
res_mslp_block_03 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = FALSE, prob = 0.3)

save(res_mslp_block_0, res_mslp_block_01, res_mslp_block_03, file = "~/phdWork/dWit/Data/climate_res_mslp_block.Rdata")


# by months
max.p <- 9
res_mslp_odd_Feb_0 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = TRUE, prob = 0, filter = months(time)=="February")

max.p <- 9
res_mslp_odd_Feb_01 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = TRUE, prob = 0.1, filter = months(time)=="February")

max.p <- 9
res_mslp_odd_Feb_03 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = TRUE, prob = 0.3, filter = months(time)=="February")

save(res_mslp_odd_Feb_0, res_mslp_odd_Feb_01, res_mslp_odd_Feb_03, file = "~/phdWork/dWit/Data/climate_res_mslp_odd_Feb.Rdata")


max.p <- 9
res_mslp_block_Feb_0 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = FALSE, prob = 0, filter = months(time)=="February")

max.p <- 9
res_mslp_block_Feb_01 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = FALSE, prob = 0.1, filter = months(time)=="February")

max.p <- 9
res_mslp_block_Feb_03 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, odd.even.split = FALSE, prob = 0.3, filter = months(time)=="February")

save(res_mslp_block_Feb_0, res_mslp_block_Feb_01, res_mslp_block_Feb_03, file = "~/phdWork/dWit/Data/climate_res_mslp_block_Feb.Rdata")



# analysis of fake:
max.p <- 9
res_fake_odd_0 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = TRUE, prob = 0)

max.p <- 9
res_fake_odd_01 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = TRUE, prob = 0.1)

max.p <- 9
res_fake_odd_03 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = TRUE, prob = 0.3)

save(res_fake_odd_0, res_fake_odd_01, res_fake_odd_03, file = "~/phdWork/dWit/Data/climate_res_fake_odd.Rdata")


max.p <- 9
res_fake_block_0 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = FALSE, prob = 0)

max.p <- 9
res_fake_block_01 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = FALSE, prob = 0.1)

max.p <- 9
res_fake_block_03 <- decendingBlocksBinClassMix(signals = "fake", max.p = max.p, odd.even.split = FALSE, prob = 0.3)

save(res_fake_block_0, res_fake_block_01, res_fake_block_03, file = "~/phdWork/dWit/Data/climate_res_fake_block.Rdata")



# analysis of all:
max.p <- 9
res_all_odd_0 <- decendingBlocksBinClassMix(signals = "all", min.p = 5, max.p = max.p, odd.even.split = TRUE, prob = 0)

max.p <- 9
res_all_odd_01 <- decendingBlocksBinClassMix(signals = "all", max.p = max.p, odd.even.split = TRUE, prob = 0.1)

max.p <- 9
res_all_odd_03 <- decendingBlocksBinClassMix(signals = "all", max.p = max.p, odd.even.split = TRUE, prob = 0.3)

save(res_all_odd_0, res_all_odd_01, res_all_odd_03, file = "~/phdWork/dWit/Data/climate_res_all_odd.Rdata")


max.p <- 9
res_all_block_0 <- decendingBlocksBinClassMix(signals = "all", max.p = max.p, odd.even.split = FALSE, prob = 0)

max.p <- 9
res_all_block_01 <- decendingBlocksBinClassMix(signals = "all", max.p = max.p, odd.even.split = FALSE, prob = 0.1)

max.p <- 9
res_all_block_03 <- decendingBlocksBinClassMix(signals = "all", max.p = max.p, odd.even.split = FALSE, prob = 0.3)

save(res_all_block_0, res_all_block_01, res_all_block_03, file = "~/phdWork/dWit/Data/climate_res_all_block.Rdata")

# plots
max.p <- 9
load("./Data/climate_res_mslp_block.Rdata")
load("./Data/climate_res_mslp_odd.Rdata")
## mslp
{

  par(mfrow=c(3,1))
plot(mslp[1,],type="l",main="mslp, 40years")
plot(mslp[1,1:(5*365)],type="l",main="mslp, 5years")
plot(mslp[1,1:(1*365)],type="l",main="mslp, 1years")

#par(mfrow=c(3,1))
#plot(main="mslp, block, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
par(mfrow=c(1,1))
plot(main="mslp, middle, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.empirical,mean), type="b", col="blue",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.empirical,quantile,0.95) , type="l", col="blue",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_0$mean.local.tv.empirical,quantile,0.05) , type="l", col="blue",pch=19,lty=2)

#plot(main="mslp, block, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

#plot(main="mslp, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
#lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="mslp, block, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="mslp, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_mslp_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

# air
par(mfrow=c(3,2))
plot(main="air, block, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="air, odd-even, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="air, block, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="air, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="air, block, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="air, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
}##

# prate
max.p <- 9
load("./Data/climate_res_prate_block.Rdata")
load("./Data/climate_res_prate_odd.Rdata")
{
  par(mfrow=c(3,1))
  plot(prate[1,],type="l",main="prate, 40years")
  plot(prate[1,1:(5*365)],type="l",main="prate, 5years")
  plot(prate[1,1:(1*365)],type="l",main="prate, 1years")


  par(mfrow=c(1,1))
  #plot(air[1,],type="l",main="air, 40years")
  plot(prate[1,],type="l",main="air, 5years")
  res <- decendingBlocksBinClassMix(signals = "prate", max.p = 1, odd.even.split = FALSE, prob = 0)
  res$mean.local.tv.search
  res$mean.local.tv.bin

par(mfrow=c(3,2))
plot(main="prate, block, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="prate, odd-even, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="prate, block, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="prate, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="prate, block, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="prate, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_prate_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
}##

# shum
max.p <- 9
load("./Data/climate_res_shum_block.Rdata")
load("./Data/climate_res_shum_odd.Rdata")
{
par(mfrow=c(3,1))
plot(shum[1,],type="l",main="shum, 40years")
plot(shum[1,1:(5*365)],type="l",main="shum, 5years")
plot(shum[1,1:(1*365)],type="l",main="shum, 1years")


par(mfrow=c(3,2))
plot(main="shum, block, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="shum, odd-even, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="shum, block, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="shum, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="shum, block, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="shum, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_shum_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
}##

# air
max.p <- 9
load("./Data/climate_res_air_block.Rdata")
load("./Data/climate_res_air_odd.Rdata")
{
par(mfrow=c(3,1))
plot(air[1,],type="l",main="air, 40years")
plot(air[1,1:(5*365)],type="l",main="air, 5years")
plot(air[1,1:(1*365)],type="l",main="air, 1years")

par(mfrow=c(1,1))
#plot(air[1,],type="l",main="air, 40years")
plot(air[1,c(57:(57+55))],type="l",main="air, 5years")
res <- decendingBlocksBinClassMix(signals = "air", max.p = 1, odd.even.split = FALSE, prob = 0, filter = 1:nrow(air)%in%c(57:(57+55)))
res$mean.local.tv.search
res$mean.local.tv.bin
#plot(air[1,1:(1*365)],type="l",main="air, 1years")

par(mfrow=c(3,2))
plot(main="air, block, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="air, odd-even, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="air, block, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="air, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="air, block, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="air, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
}##

# fake
max.p <- 9
load("./Data/climate_res_fake_block.Rdata")
load("./Data/climate_res_fake_odd.Rdata")
{
par(mfrow=c(3,1))
plot(fake.series[1,],type="l",main="air, 40years")
plot(fake.series[1,1:(5*365)],type="l",main="air, 5years")
plot(fake.series[1,1:(1*365)],type="l",main="air, 1years")

par(mfrow=c(3,2))
plot(main="fake, block, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="fake, odd-even, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="fake, block, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="fake, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


plot(main="fake, block, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

plot(main="fake, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
lines(log(length(time)/2^{1:max.p}/365),lapply(res_fake_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)
}##

# air (Feb only)
max.p <- 5
load("./Data/climate_res_air_block_Feb.Rdata")
load("./Data/climate_res_air_odd_Feb.Rdata")
{par(mfrow=c(3,2))
  plot(main="air, block, prob=0, Feb",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

  plot(main="air, odd-even, prob=0, Feb",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


  plot(main="air, block, prob=0.1, Feb",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

  plot(main="air, odd-even, prob=0.1, Feb",log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_odd_Feb_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


  plot(main="air, block, prob=0.3, Feb", log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

  plot(main="air, odd-even, prob=0.3, Feb", log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_air_block_Feb_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)}


# all (odd-even only)
max.p <- 9
load("./Data/climate_res_all_odd.Rdata")
{par(mfrow=c(3,1))

  plot(main="all, odd-even, prob=0",log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_0$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_0$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_0$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_0$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_0$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_0$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)


  plot(main="all, odd-even, prob=0.1",log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_01$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_01$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_01$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_01$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_01$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_01$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)

  plot(main="all, odd-even, prob=0.3", log(length(time)/2^{1:max.p}/365),lapply(res_all_odd_03$mean.local.tv.search,quantile,0.95) , type="l", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_block_03$mean.local.tv.search,mean), type="b", col="black",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_block_03$mean.local.tv.search,quantile,0.05) , type="l", col="black",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_block_03$mean.local.tv.bin,mean), type="b", col="red",pch=19)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_block_03$mean.local.tv.bin,quantile,0.95) , type="l", col="red",pch=19,lty=2)
  lines(log(length(time)/2^{1:max.p}/365),lapply(res_all_block_03$mean.local.tv.bin,quantile,0.05) , type="l", col="red",pch=19,lty=2)}

# length in days
length(time)/2^{1:max.p}







# Meeting with Nicolai and Jeff

# - corrected split test train
# - run the analysis for all locations together
# - some alone
# analysis of mslp:

# all locations, rf
max.p <- 9
res_mslp_all_rf_5 <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = FALSE)

# all locations, lda
max.p <- 9
res_mslp_all_lda <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = TRUE)

# one location, rf
max.p <- 9
res_mslp_1_rf <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = FALSE, locations = 1)

# one location, lda
max.p <- 9
res_mslp_1_lda <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = TRUE, locations = 1)


# two location, rf
max.p <- 9
res_mslp_2_rf <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = FALSE, locations = 1:2)

# two location, lda
max.p <- 9
res_mslp_2_lda <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = TRUE, locations = 1:2)

# ten location, rf
max.p <- 9
res_mslp_10_rf <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = FALSE, locations = 1:10)

# ten location, lda
max.p <- 9
res_mslp_10_lda <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = TRUE, locations = 1:10)

# 500 location, rf
max.p <- 9
res_mslp_500_rf <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = FALSE, locations = 1:500)

# 500 location, lda
max.p <- 9
res_mslp_500_lda <- decendingBlocksBinClassMix(signals = "mslp", max.p = max.p, min.p = 5, prob = 0, use.lda = TRUE, locations = 1:500)


save(res_mslp_all_rf, res_mslp_all_lda, res_mslp_500_rf, res_mslp_500_lda, res_mslp_2_rf, res_mslp_2_lda, res_mslp_10_rf, res_mslp_10_lda, file = "./Data/res_mslp.Rdata")




# vizu

# two locations
par(mfrow=c(4,2))
min.p <- 5
plot(main="2 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_2_rf$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_2_rf$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_2_rf$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)

min.p <- 5
plot(main="2 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_2_lda$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_2_lda$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_2_lda$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)

# ten locations
min.p <- 5
plot(main="10 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_10_rf$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_10_rf$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_10_rf$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)

min.p <- 5
plot(main="10 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_10_lda$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_10_lda$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_10_lda$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)

# 500 locations
min.p <- 5
plot(main="500 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_500_rf$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_500_rf$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_500_rf$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)

min.p <- 5
plot(main="500 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_500_lda$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_500_lda$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_500_lda$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)


# all locations
min.p <- 5
plot(main="all locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_all_rf$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_all_rf$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_all_rf$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)

min.p <- 5
plot(main="all locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_all_lda$mean.local.tv.search,mean)[5:9], type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_all_lda$mean.local.tv.bin,mean)[5:9], type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(res_mslp_all_lda$mean.local.tv.empirical,mean)[5:9], type="b", col="blue", pch=19)




# small simulation to see whats happening


# move in the core of the distribution
n <- 500
x <- c(rnorm(n),rnorm(n,2))
t <- c(rep(0,n),rep(1,n))
rhofun <- function(x) dnorm(x,2) / (dnorm(x,2) + dnorm(x))
rho <- rhofun(x)
plot(density(x[1:n]),xlim=c(-10,10))
lines(density(x[-c(1:n)]), col="red")
par(mfrow=c(1,1))
d <- dWit(t = t, rho = rho, s = 0.5, estimator.type = "asymptotic-tv-search", verbose.plot = TRUE)


# move more localized
n <- 500
x <- c(rnorm(n),ifelse(runif(n) > 0.95, rnorm(x,4), rnorm(x)))
t <- c(rep(0,n),rep(1,n))
rhofun <- function(x) (0.1*dnorm(x,4)+0.95*dnorm(x)) / (0.05*dnorm(x,4) + (1.95)*dnorm(x))
rho <- rhofun(x)
plot(density(x[1:n]),xlim=c(-10,10))
lines(density(x[-c(1:n)]), col="red")
plot(rhofun, -10, 10)
par(mfrow=c(1,1))
d <- dWit(t = t, rho = rho, s = 0.5, estimator.type = "asymptotic-tv-search", verbose.plot = TRUE)
d
d <- dWit(t = t, rho = rho, s = 0.5, estimator.type = "binomial", verbose.plot = TRUE)
d


# take largest block compare with lienar and rf, for each day subtract the mean spatially (depends on all location or not?) probably non-zero
# might shift to tails?

# Objective for tuesday
@Loris
# - fucking finish the climatic data. same analysis but remove spatial mean
# - do smth with the cifar10 that is interesting to show. ...
# - (I consolidate the rest if the simulations). ...
# - sims for the gaps of lowerbound. Simulations similar to power but for (consistent)
# - sims supportive for corollary. (Adapt on the way)

@Jeff, @Loris
# - Move on paper for binomial writing (currently done by Jeff)
# - Proof of consistency (we can do it quickly tomorrow)



# final runs
min.p <- 1
max.p <- 9
final_mslp_rf <- runClimateAnalysis(signals = "mslp", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE)

min.p <- 1
max.p <- 9
final_air_rf <- runClimateAnalysis(signals = "air", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE)

min.p <- 1
max.p <- 9
final_shum_rf <- runClimateAnalysis(signals = "shum", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE)


min.p <- 1
max.p <- 9
final_prate_rf <- runClimateAnalysis(signals = "prate", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE)


# final runs
min.p <- 1
max.p <- 9
final_mslp_lda <- runClimateAnalysis(signals = "mslp", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE)

min.p <- 1
max.p <- 9
final_air_lda <- runClimateAnalysis(signals = "air", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE)

min.p <- 1
max.p <- 9
final_shum_lda <- runClimateAnalysis(signals = "shum", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE)


min.p <- 1
max.p <- 9
final_prate_lda <- runClimateAnalysis(signals = "prate", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE)



save(final_air_lda, final_mslp_lda, final_prate_rf, final_shum_rf, final_air_rf, final_mslp_rf, file = "./Data/climate_final.Rdata")


# plots of what we have
plot(main="all locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)



# final runs
min.p <- 1
max.p <- 9
final_mslp_rf_10 <- runClimateAnalysis(signals = "mslp", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:10)

min.p <- 1
max.p <- 9
final_air_rf_10 <- runClimateAnalysis(signals = "air", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:10)

min.p <- 1
max.p <- 9
final_shum_rf_10 <- runClimateAnalysis(signals = "shum", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:10)


min.p <- 1
max.p <- 9
final_prate_rf_10 <- runClimateAnalysis(signals = "prate", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:10)


# final runs
min.p <- 1
max.p <- 9
final_mslp_lda_10 <- runClimateAnalysis(signals = "mslp", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:10)

min.p <- 1
max.p <- 9
final_air_lda_10 <- runClimateAnalysis(signals = "air", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:10)

min.p <- 1
max.p <- 9
final_shum_lda_10 <- runClimateAnalysis(signals = "shum", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:10)


min.p <- 1
max.p <- 9
final_prate_lda_10 <- runClimateAnalysis(signals = "prate", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:10)


save(final_air_lda_10, final_mslp_lda_10, final_prate_rf_10, final_shum_rf_10, final_air_rf_10, final_mslp_rf_10, file = "./Data/climate_final_10.Rdata")



# final runs
min.p <- 1
max.p <- 9
final_mslp_rf_500 <- runClimateAnalysis(signals = "mslp", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:500)

min.p <- 1
max.p <- 9
final_air_rf_500 <- runClimateAnalysis(signals = "air", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:500)

min.p <- 1
max.p <- 9
final_shum_rf_500 <- runClimateAnalysis(signals = "shum", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:500)


min.p <- 1
max.p <- 9
final_prate_rf_500 <- runClimateAnalysis(signals = "prate", max.p = max.p, min.p = min.p, prob = 0, use.lda = FALSE, spatial.centering = TRUE, locations = 1:500)


# final runs
min.p <- 1
max.p <- 9
final_mslp_lda_500 <- runClimateAnalysis(signals = "mslp", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:500)

min.p <- 1
max.p <- 9
final_air_lda_500 <- runClimateAnalysis(signals = "air", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:500)

min.p <- 1
max.p <- 9
final_shum_lda_500 <- runClimateAnalysis(signals = "shum", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:500)


min.p <- 1
max.p <- 9
final_prate_lda_500 <- runClimateAnalysis(signals = "prate", max.p = max.p, min.p = min.p, prob = 0, use.lda = TRUE, spatial.centering = TRUE, locations = 1:500)

save(final_air_lda_500, final_mslp_lda_500, final_prate_rf_500, final_shum_rf_500, final_air_rf_500, final_mslp_rf_500, file = "./Data/climate_final_500.Rdata")


n <- 50
p <- 10
X <- matrix(rnorm(n * p), n, p)
Y <- X[, 1] * rnorm(n)
c.forest <- custom_forest(X, Y)




# pressure
par(mfrow=c(3,2))
plot(main="pressure, all locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="pressure, all locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)


plot(main="pressure, 500 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf_500$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf_500$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf_500$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="pressure, 500 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda_500$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda_500$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda_500$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)


plot(main="pressure, 10 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf_10$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf_10$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_rf_10$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="pressure, 10 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda_10$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda_10$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_mslp_lda_10$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)


# temperature
par(mfrow=c(3,2))
plot(main="temperature, all locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="temperature, all locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)


plot(main="temperature, 500 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf_500$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf_500$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf_500$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="temperature, 500 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda_500$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda_500$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda_500$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)


plot(main="temperature, 10 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf_10$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf_10$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_rf_10$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="temperature, 10 locations, LDA",log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda_10$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda_10$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_air_lda_10$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)


par(mfrow=c(3,1))
plot(main="rainfall, all locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="rainfall, 500 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf_500$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf_500$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf_500$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

plot(main="rainfall, 10 locations, RF",log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf_10$mean.local.tv.search,mean), type="b", col="black",pch=19,lty=2,ylim=c(0,1),xlab="block size in log(years)",ylab="TV")
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf_10$mean.local.tv.bin,mean), type="b", col="red",pch=19)
lines(log(length(time)/2^{min.p:max.p}/365),lapply(final_prate_rf_10$mean.local.tv.empirical,mean), type="b", col="blue", pch=19)

b <- c(90,0)
lon.pts <- c(0,90)
lat.pts <- c(90,180)

lon.center <- c(8.2275)
lat.center <- c(46.8182)

genCoordXY <- function(deltaX = 1, deltaY = 1, nb.steps = 0) {
  return(expand.grid(lon.center+seq(-nb.steps*deltaX, nb.steps*deltaX, deltaX), lat.center+seq(-nb.steps*deltaY, nb.steps*deltaY, deltaY)))
}

for (nb in 25:30) {
  xy.coords <- genCoordXY(nb.steps = nb)
  print(dim(e <- extract(air_raster, xy.coords,method="bilinear")))
}
xy.center <- cbind(lon.center,lat.center)
e <- extract(air_raster, xy.coords,method="bilinear")
dim(e)
