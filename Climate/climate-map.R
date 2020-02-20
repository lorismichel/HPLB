# Climatic map analysis

# options
PATH.CLIMATE.DATA <- "~/Downloads/reanalysis_jeff_loris/"


PRODUCE.PLOTS       <- FALSE
RUN.ANALYSIS        <- FALSE
RUN.COMBINED.METHOD <- TRUE

PREPRO <- 2
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
d <- climatePrePro(PATH.CLIMATE.DATA)


# create the data
data       <- list()
data$air   <- d$air
data$prate <- d$prate
data$shum  <- d$shum
data$mslp  <- d$mslp

# preprocessing
if (PREPRO == 1) {
  data$air   <- apply(data$air, 1, rank)
  data$shum  <- apply(data$shum, 1, rank)
  data$prate <- apply(data$prate, 1, rank)
  data$mslp  <- apply(data$mslp, 1, rank)
} else if (PREPRO == 2) {
  data$air   <- apply(data$air, 1, diff)
  data$shum  <- apply(data$shum, 1, diff)
  data$prate <- apply(data$prate, 1, diff)
  data$mslp  <- apply(data$mslp, 1, diff)
} else if (PREPRO == 3) {
  data$air   <- apply(data$air, 1, function(x) rank(diff(x)))
  data$shum  <- apply(data$shum, 1, function(x) rank(diff(x)))
  data$prate <- apply(data$prate, 1, function(x) rank(diff(x)))
  data$mslp  <- apply(data$mslp, 1, function(x) rank(diff(x)))
}


# getting the dimensions
p <- ncol(data$air)
n <- nrow(data$air)



# results
resmat <- matrix(nrow=p,
                 ncol = ifelse(RUN.COMBINED.METHOD,
                               3,
                               2))
# tv seq
tvseq <- sort(unique(c(seq(0,1,by=0.01),
                       seq(0,0.1,by=0.001),
                       seq(0,0.2,by=0.005),
                       c(0,0.0001,0.001,0.002,0.005,0.01,0.015,0.02))))

if (RUN.ANALYSIS) {

  # main for loop
  for (k in 1:p){

    X <- cbind(data$air[,k],
               data$mslp[,k],
               data$prate[,k],
               data$shum[,k])

    # splits train-test
    if (SPLIT.TRAIN.TEST == 0) {
      ind.train <- which(1:n%%2 == 0)
      ind.test <- which(1:n%%2 == 1)
    } else if (SPLIT.TRAIN.TEST == 1) {
      ind.train <- (round(n/4):round(n*3/4))
      ind.test <- (1:n)[-ind.train]
    }

    # creating the data and fitting forest
    Y <- as.factor((1:n)>round(n/2))
    rf <- ranger(formula = y~., data = data.frame(x=X[ind.train,], y=Y[ind.train]), probability = TRUE)
    pred <- predict(rf, data = data.frame(x=X[ind.test,]))$predictions[,"TRUE"]

    tvhat_bin <- dWit(t=as.numeric(Y[ind.test])-1, rho = pred, tv.seq  = tvseq,
                      estimator.type = "binomial-test")$tvhat
    tvhat_search <- dWit(t=as.numeric(Y[ind.test])-1, rho = pred, tv.seq  = tvseq,
                         estimator.type = "asymptotic-tv-search")$tvhat

    if (RUN.COMBINED.METHOD) {
      tvhat_combined <- dWit(t=as.numeric(Y[ind.test])-1, rho = pred,
                             tv.seq  = c(tvhat_bin, Filter(x = tvseq, f = function(x) x > tvhat_bin)),
                             estimator.type = "asymptotic-tv-search")$tvhat
      resmat[k,] <- c(tvhat_bin,
                      tvhat_search,
                      tvhat_combined)
    } else {
      resmat[k,] <- c(tvhat_bin,
                      tvhat_search)
    }

    # verbose
    if (k%%10 == 0) {
      print(k)
    }
  }

  # saving the data
  save(resmat,
       file = paste0("Data/DATA_CLIMATE_MAP_SPLIT_",
                     SPLIT.TRAIN.TEST,
                     "_PREPRO_",
                     PREPRO,
                     ".Rdata"))

}



if (PRODUCE.PLOTS) {

  # load the data
  info <- load(paste0("Data/DATA_CLIMATE_MAP_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",
                      PREPRO,
                      ".Rdata"))

  # produce the plots
  library("ggmap")
  library(maptools)
  library(maps)

  # coords
  x <-  d$matrix.coords[,1]
  y <-  d$matrix.coords[,2]

  png(filename = paste0("Plots/PLOT_CLIMATE_MAP_SPLIT_",
                        SPLIT.TRAIN.TEST,
                        "_PREPRO_",
                        PREPRO,
                        "_PLOT_1.png"),
      width = 1000)
  ## binomial tv
  mp_bin <- NULL
  mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
  mp_bin <- ggplot() +   mapWorld

  tv <- resmat[,1]
  mp_bin <- mp_bin + geom_point(aes(x=x, y=y, color = tv), size=2) + scale_colour_gradient(low = "white", high="black")
  mp_bin
  dev.off()

  png(filename = paste0("Plots/PLOT_CLIMATE_MAP_SPLIT_",
                        SPLIT.TRAIN.TEST, "_PREPRO_",
                        PREPRO, "_PLOT_2.png"), width = 1000)
  ## tvsearch tv
  mp_search <- NULL
  mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
  mp_search <- ggplot() +   mapWorld

  tv <- resmat[,2]
  mp_search <- mp_search + geom_point(aes(x=x, y=y, color = tv), size=2) + scale_colour_gradient(low = "white", high="black")
  mp_search
  dev.off()

  png(filename = paste0("Plots/PLOT_CLIMATE_MAP_SPLIT_", SPLIT.TRAIN.TEST, "_PREPRO_", PREPRO, "_PLOT_3"), width = 1000)
  plot(resmat[,1], resmat[,2], col="black",pch=19,xlab="binomial-test",ylab="asymptotic-tv-search",font.lab=1,font.main=1, cex = 0.7)
  abline(0,1,col="blue",lty=2)
  dev.off()


  if (RUN.COMBINED.METHOD) {
    mp_bin <- NULL
    mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
    mp_bin <- ggplot() +   mapWorld

    tv <- resmat[,3]
    mp_bin <- mp_bin + geom_point(aes(x=x, y=y, color = tv), size=2) + scale_colour_gradient(low = "white", high="black")
    mp_bin
    dev.off()

    png(filename = paste0("Plots/PLOT_CLIMATE_MAP_SPLIT_",
                          SPLIT.TRAIN.TEST,
                          "_PREPRO_", PREPRO,
                          "_PLOT_3.png"), width = 1000)
    # # tvsearch tv
    mp_search <- NULL
    mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
    mp_search <- ggplot() +   mapWorld

  }


}
