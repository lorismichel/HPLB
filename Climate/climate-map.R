# climate map

# options
RANK.TRANSFORMATION <- FALSE
SPLIT.TRAIN.TEST <- 1

# source
source("./Climate/climate-preprocessing.R")

# repro
set.seed(1)

# libs
require(dWit)
require(fields)
require(ranger)

# prepro
d <- climatePrePro()


# create the data
data <- list()
data$air <- d$air
data$prate <- d$prate
data$shum <- d$shum
data$mslp <- d$mslp


# transformations
if (RANK.TRANSFORMATION) {
  data$air   <- apply(data$air, 2, rank)
  shum  <- apply(data$shum, 2, rank)
  prate <- apply(data$prate, 2, rank)
  mslp  <- apply(data$mslp, 2, rank)
}


# looping
p <- nrow(data$air)
n <- ncol(data$air)-1

resmat <- matrix(nrow=p,ncol=2)
tvseq <- sort(unique(c(seq(0,1,by=0.01),seq(0,0.1,by=0.001),seq(0,0.2,by=0.005),c(0,0.0001,0.001,0.002,0.005,0.01,0.015,0.02))))

#sample(1:p,p)
for (k in 283:p){
#for (k in 1:20) {
  X <- apply(cbind( data$air[k,], data$mslp[k,], data$prate[k,], data$shum[k,]),2,diff)

  # splits train-test
  if (SPLIT.TRAIN.TEST == 0) {
    ind.train <- which(1:nrow(dat)%%2 == 0)
    ind.test <- which(1:nrow(dat)%%2 == 1)
  } else if (SPLIT.TRAIN.TEST == 1) {
    n <- length(air)
    ind.train <- (round(n/4):round(n*3/4))
    ind.test <- (1:n)[-ind.train]
  }

  Y <- as.factor((1:n)>round(n/2))

  rf <- ranger(formula = y~., data = data.frame(x=X[ind.train,], y=Y[ind.train]), probability = TRUE)

  pred <- predict(rf, data = data.frame(x=X[ind.test,]))$predictions[,"TRUE"]

  tvhat_bin <- dWit( t=as.numeric(Y[ind.test])-1, rho = pred, tv.seq  = tvseq, estimator.type = "binomial-test")$tvhat
  tvhat_search <- dWit( t=as.numeric(Y[ind.test])-1, rho = pred, tv.seq  = tvseq, estimator.type = "asymptotic-tv-search")$tvhat

  resmat[k,] <- c(tvhat_bin, tvhat_search)

  if (k%%10==0) {
    print(k)
  }
}


# get coordinates
matrix.coords <- matrix(nrow=nrow(d$shum),ncol=2)

for (i in 1:nrow(d$shum)) {
  matrix.coords[i,] <- xyFromCell(object = d$shum_raster, i)
}

# saving the data
save(resmat, file = paste0("Data/DATA_CLIMATE_MAP_", RANK.TRANSFORMATION, "_", SPLIT.TRAIN.TEST, ".Rdata"))



# # producing plots
# info <- load("Data/DATA_CLIMATE_MAP.Rdata")
#
# # pdf(file="results_differenced.pdf",width=12,height=6)
# # par(mfrow=c(1,3))
# # plot(resmat[1:k,1],resmat[1:k,2],xlab="TV-bound (conf table)",ylab="TV-bound")
# # abline(c(0,1))
# # plot(matrix.coords[,1]-180, matrix.coords[,2], cex=rexp(nrow(matrix.coords)),xlab="longitute",ylab="latitude",font.main=1,font.lab=1)
# # plot(matrix.coords[,1]-180, matrix.coords[,2], cex=rexp(nrow(matrix.coords)),xlab="longitute",ylab="latitude",font.main=1,font.lab=1)
# #
# # i#mage.plot( matrix( resmat[,1],nrow=72),axes=FALSE,main="TV-bound (conf table)")
# # #image.plot( matrix( resmat[,2],nrow=72),axes=FALSE,main="TV-bound")
# # dev.off()
# #
# #
# #Now Layer the cities on top
# library("ggmap")
# library(maptools)
# library(maps)
# #
# #
# x <-  matrix.coords[,1]-180
# y <-  matrix.coords[,2]
#
# # binomial
# #Using GGPLOT, plot the Base World Map
# mp_bin <- NULL
# mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
# mp_bin <- ggplot() +   mapWorld
#
# #Now Layer the cities on top
# tv <- resmat[,1]
# mp_bin <- mp_bin + geom_point(aes(x=x, y=y, color = tv), size=3) + scale_colour_gradient(low = "white", high="black")
# mp_bin
#
# # tvsearch
# mp_search <- NULL
# mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
# mp_search <- ggplot() +   mapWorld
#
# #Now Layer the cities on top
# mp_search <- mp_search + geom_point(aes(x=x, y=y, color = resmat[,2]), size=1) + scale_colour_gradient(low = "red", high="blue")
# mp_search
#
#
# plot(resmat[,1], resmat[,2], col="black",pch=19,xlab="binomial-test",ylab="asymptotic-tv-search",font.lab=1,font.main=1, cex = 0.7)
# abline(0,1,col="blue",lty=2)
#
