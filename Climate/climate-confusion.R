# confusion matrix 4 decades

## options
PATH.CLIMATE.DATA <- "~/Downloads/reanalysis_jeff_loris/"
PREPRO <- 0
SPLIT.TRAIN.TEST <- 0

# source
source("./Climate/climate-preprocessing.R")

# repro
set.seed(1)

# libs
require(dWit)

# preprocessing and loading of the data
d <- climatePrePro(path = PATH.CLIMATE.DATA)

# zurich coordinates and data
zurich.coord <- c(8.5391825, 47.3686498)
air   <- extract(d$air_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]
shum  <- extract(d$shum_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]
prate <- extract(d$prate_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]
mslp  <- extract(d$mslp_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]

# preprocessing
if (PREPRO == 1) {
  air   <- rank(air)
  shum  <- rank(shum)
  prate <- rank(prate)
  mslp  <- rank(mslp)
} else if (PREPRO == 2) {
  air   <- diff(air)
  shum  <- diff(shum)
  prate <- diff(prate)
  mslp  <- diff(mslp)
} else if (PREPRO == 3) {
  air   <- rank(diff(air))
  shum  <- rank(diff(shum))
  prate <- rank(diff(prate))
  mslp  <- rank(diff(mslp))
}

## core of analysis
dat <- data.frame(class = d$tenyears.ind,
                  x = cbind(t(air),t(mslp),t(prate),t(shum)))

# run analysis

# splits train-test
if (SPLIT.TRAIN.TEST == 0) {
  ind.train <- which(1:nrow(dat)%%2 == 0)
  ind.test <- which(1:nrow(dat)%%2 == 1)
} else if (SPLIT.TRAIN.TEST == 1) {
  n <- length(air)
  ind.train <- (round(n/4):round(n*3/4))
  ind.test <- (1:n)[-ind.train]
}


# define the multiclass by quantile splits
b <- as.numeric(quantile(1:length(air), c(seq(0,1,length.out = 5))))
b[1] <- 0
class <- cut(1:length(air), breaks = b, labels = FALSE)

# fit a forest
mRF_air <- ranger::ranger(class~air, data = data.frame(class = class, air=air)[ind.train,], probability = TRUE)
mRF_mslp <- ranger::ranger(class~mslp, data = data.frame(class = class, mslp=mslp)[ind.train,],  probability = TRUE)
mRF_prate <- ranger::ranger(class~prate, data = data.frame(class = class, prate=prate)[ind.train,],  probability = TRUE)
mRF_shum <- ranger::ranger(class~shum, data = data.frame(class = class, shum=shum)[ind.train,],  probability = TRUE)
mRF_joint <- ranger::ranger(class~., data = data.frame(class = class, air=air, mslp=mslp, prate=prate, shum=shum)[ind.train,], importance = 'permutation',  probability = TRUE)


# construct the test set
labels <- (class - 1)[-ind.train]
dat.test <- data.frame(class = class, t=1:length(shum), air=air, mslp=mslp, prate=prate, shum=shum)[-ind.train,]


ordering.array.air <- array(dim=c(4,4,nrow(dat.test)))
ordering.array.mslp <- array(dim=c(4,4,nrow(dat.test)))
ordering.array.prate <- array(dim=c(4,4,nrow(dat.test)))
ordering.array.shum <- array(dim=c(4,4,nrow(dat.test)))
ordering.array.joint <- array(dim=c(4,4,nrow(dat.test)))


# build the rho function
for (i in 1:4) {
  for (j in 1:4) {
    ordering.array.air[i,j,] <- predict(mRF_air, data = dat.test)$predictions[,j]-predict(mRF_air, data = dat.test)$predictions[,i]
    ordering.array.mslp[i,j,] <- predict(mRF_mslp, data = dat.test)$predictions[,j]-predict(mRF_mslp, data = dat.test)$predictions[,i]
    ordering.array.prate[i,j,] <- predict(mRF_prate, data = dat.test)$predictions[,j]-predict(mRF_prate, data = dat.test)$predictions[,i]
    ordering.array.shum[i,j,] <- predict(mRF_shum, data = dat.test)$predictions[,j]-predict(mRF_shum, data = dat.test)$predictions[,i]
    ordering.array.joint[i,j,] <- predict(mRF_joint, data = dat.test)$predictions[,j]-predict(mRF_joint, data = dat.test)$predictions[,i]
  }
  print(i)
}

# get the TV matrices
tv.mat.air   <- dWitMatrix(labels = labels, ordering.array = ordering.array.air)
tv.mat.mslp  <- dWitMatrix(labels = labels, ordering.array = ordering.array.mslp)
tv.mat.prate <- dWitMatrix(labels = labels, ordering.array = ordering.array.prate)
tv.mat.shum  <- dWitMatrix(labels = labels, ordering.array = ordering.array.shum)
tv.mat.joint <- dWitMatrix(labels = labels, ordering.array = ordering.array.joint)

save(tv.mat.air,
     tv.mat.mslp,
     tv.mat.prate,
     tv.mat.shum,
     tv.mat.joint,
     file = paste0("./Data/DATA_CLIMATE_CONFUSION_SPLIT", SPLIT.TRAIN.TEST,"_PREPRO_",PREPRO,".Rdata"))
