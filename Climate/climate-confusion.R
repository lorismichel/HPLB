# confusion matrix 4 decades

## options
PATH.CLIMATE.DATA <- "~/Downloads/reanalysis_jeff_loris/"
PREPRO <- 2
SPLIT.TRAIN.TEST <- 1

# source
source("./Climate/climate-preprocessing.R")

# repro
set.seed(1)

# libs
require(dWit)

# preprocessing and loading of the data
d <- climatePrePro(path = PATH.CLIMATE.DATA)

# location coordinates and data
loc.coord <- c(8.5391825, 47.3686498)
air   <- extract(d$air_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]
shum  <- extract(d$shum_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]
prate <- extract(d$prate_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]
mslp  <- extract(d$mslp_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]

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

# run analysis

# splits train-test
if (SPLIT.TRAIN.TEST == 0) {

  ind.train <- which(1:nrow(dat)%%2 == 0)
  ind.test <- which(1:nrow(dat)%%2 == 1)
} else if (SPLIT.TRAIN.TEST == 1) {

  n <- length(air)
  ind <- cut(1:n, c(-1, split.ids, n+1), FALSE)
  ind.train <- Filter(x = 1:n,f =  function(i) {
    s <- which(ind==ind[i])[1]
    l <- length(which(ind==ind[i]))
    i %in% ((round(s + l/4):round(s + l*3/4)))
  })
  ind.test <- (1:n)[-ind.train]
}


# define the class by quantile splits
b <- as.numeric(quantile(1:length(air), c(seq(0,1,length.out = 5))))
b[1] <- 0
class <- cut(1:length(air), breaks = b, labels = FALSE)

# create train and test sets
train <- data.frame(class=class, air = air, mslp = mslp, shum = shum, prate = prate)[ind.train,]
test  <- data.frame(class=class, air = air, mslp = mslp, shum = shum, prate = prate)[ind.test,]


# fit a forest
mRF_air <- ranger::ranger(class~air, data = train, probability = TRUE)
mRF_mslp <- ranger::ranger(class~mslp, data = train,  probability = TRUE)
mRF_prate <- ranger::ranger(class~prate, data = train,  probability = TRUE)
mRF_shum <- ranger::ranger(class~shum, data = train,  probability = TRUE)
mRF_joint <- ranger::ranger(class~., data = train, importance = 'permutation',  probability = TRUE)


# construct the test set labels
labels <- (class - 1)[-ind.train]


ordering.array.air <- array(dim=c(4,4,nrow(test)))
ordering.array.mslp <- array(dim=c(4,4,nrow(test)))
ordering.array.prate <- array(dim=c(4,4,nrow(test)))
ordering.array.shum <- array(dim=c(4,4,nrow(test)))
ordering.array.joint <- array(dim=c(4,4,nrow(test)))


# getting prediction and building the rho
preds_air <- predict(mRF_air, data = test)$predictions
preds_mslp <- predict(mRF_mslp, data = test)$predictions
preds_prate <- predict(mRF_prate, data = test)$predictions
preds_shum <- predict(mRF_shum, data = test)$predictions
preds_joint <- predict(mRF_joint, data = test)$predictions


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array.air[i,j,] <- preds_air[,j]-preds_air[,i]
    ordering.array.mslp[i,j,] <- preds_mslp[,j]-preds_mslp[,i]
    ordering.array.prate[i,j,] <- preds_prate[,j]-preds_prate[,i]
    ordering.array.shum[i,j,] <- preds_shum[,j]-preds_shum[,i]
    ordering.array.joint[i,j,] <- preds_joint[,j]-preds_joint[,i]
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
     file = paste0("./Data/DATA_CLIMATE_CONFUSION_SPLIT",
                   SPLIT.TRAIN.TEST,
                   "_PREPRO_",
                   PREPRO,
                   ".Rdata"))
