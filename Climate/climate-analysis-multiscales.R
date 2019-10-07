# Multi-scales analysis

# source
source("./climate-preprocessing.R")

runClimateAnalysis <- function(signals = c("air"),
                               max.p = 9, min.p = 1, filter = 1:length(time),
                               prob = 0, locations = c(1:nrow(air)),
                               use.lda = FALSE,
                               spatial.centering = TRUE,
                               x = NULL) {
  if (is.null(x)) {
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

      tv.bin <- c(tv.bin , dWit(t = class.ind.test, rho = if (use.lda) preds else preds[,2], estimator.type = "binomial-test")$tvhat)
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

# prepro
d <- climatePrePro()

# generate params
lon.center <- c(8.2275)
lat.center <- c(46.8182)

genCoordXY <- function(deltaX = 1, deltaY = 1, nb.steps = 0) {
  return(expand.grid(lon.center+seq(-nb.steps*deltaX, nb.steps*deltaX, deltaX), lat.center+seq(-nb.steps*deltaY, nb.steps*deltaY, deltaY)))
}


# lopping
runC.list <- list()
nb.loc.list <- list()
for (nb in 1:10) {

  # get the locations
  xy.coords <- genCoordXY(nb.steps = nb)

  # get the series
  x <- extract(d$prate_raster, xy.coords,method="bilinear")
  x <- x[,1:14641]
  ind <- apply(x,1,function(xx) any(is.na(xx)))
  x <- x[which(!ind),]
  nb.loc.list[[i]] <- nrow(x)
  runC.list[[nb]] <- runClimateAnalysis(x = t(x), max.p = 1, min.p = 1, use.lda = FALSE)
  print(nb)
}

search <- unlist(lapply(runC.list, function(l) l$mean.local.tv.search))
bin <- unlist(lapply(runC.list, function(l) l$mean.local.tv.bin))

# plotting
par(mfrow=c(1,1))
plot(search,type="b",ylim=c(0,1))
lines(bin,col="red")

