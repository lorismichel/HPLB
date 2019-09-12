library(sf)
library(ncdf4)
library(raster)
library(rasterVis)
library(RColorBrewer)

# set path and filename
ncpath <- "~/Downloads/reanalysis_jeff_loris/"
#dname <- "air"  # note: tmp means temperature (not temporary)

air_raster <- brick(paste(ncpath, "NCEP2.air.2m.day", ".nc", sep=""), varname="air")
mslp_raster <- brick(paste(ncpath, "NCEP2.mslp.2m.day", ".nc", sep=""), varname="mslp")
prate_raster <- brick(paste(ncpath, "NCEP2.prate.2m.day", ".nc", sep=""), varname="prate")
shum_raster <- brick(paste(ncpath, "NCEP2.shum.2m.day", ".nc", sep=""), varname="shum")

t_air <- colnames(air_raster[1,1])

toDate <- function(x) {
  return(sub(sub(substr(x, start = 2, stop = 11), pattern = "\\.", replacement = "-"), pattern = "\\.", replacement = "-"))
}


t_air <- as.Date(sapply(colnames(air_raster[1,1]), toDate))
t_mslp <- as.Date(sapply(colnames(mslp_raster[1,1]), toDate))
t_prate <- as.Date(sapply(colnames(prate_raster[1,1]), toDate))
t_shum <- as.Date(sapply(colnames(shum_raster[1,1]), toDate))

# look at length
length(intersect(t_prate, t_air))
length(intersect(t_prate, t_mslp))
length(intersect(t_prate, t_shum))
length(t_air)
length(t_mslp)
length(t_shum)

# look at range
range(t_prate)
range(t_air)
range(t_mslp)
range(t_shum)




any(is.na(air_raster))

air <- as.array(air_raster)
dim(air_raster)
dim(air)
mslp <- as.array(mslp_raster)
dim(mslp_raster)
dim(mslp)
prate <- as.array(prate_raster)
dim(prate_raster)
dim(prate)
shum <- as.array(shum_raster)
dim(shum_raster)
dim(shum)


# subset to smaller date
air <- getValues(air_raster)[,1:14641]
mslp <- getValues(mslp_raster)[,1:14641]
prate <- getValues(prate_raster)[,1:14641]
shum <- getValues(shum_raster)[,1:14641]

# check nas
any(is.na(air))
any(is.na(mslp))
any(is.na(prate))
any(is.na(shum))


#tmp_array <- getValues(tmp_raster)
#dim(tmp_array)

#length(tmp_raster[1,1][1])
#tmp_raster[1,1][3]
#length(tmp_raster)



# we can look at time
#colnames(tmp_raster[1,1])

# create years time
time <- t_air
tenyears.ind <- cut(time, breaks = as.Date(c("1979-01-01","1989-01-01","1999-01-01","2009-01-01","2019-02-01")), labels=FALSE)



# just air
dat <- data.frame(class=tenyears.ind, x = cbind(t(air)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)

# joint analysis
dat <- data.frame(class=tenyears.ind, x = cbind(t(air),t(mslp),t(prate),t(shum)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)
tv.mat


# just mslp
dat <- data.frame(class=tenyears.ind, x = cbind(t(mslp)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)

# joint analysis
dat <- data.frame(class=tenyears.ind, x = cbind(t(air),t(mslp),t(prate),t(shum)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)
tv.mat



# just prate
dat <- data.frame(class=tenyears.ind, x = cbind(t(mslp)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)

# joint analysis
dat <- data.frame(class=tenyears.ind, x = cbind(t(air),t(mslp),t(prate),t(shum)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)
tv.mat


# just shum
dat <- data.frame(class=tenyears.ind, x = cbind(t(mslp)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)

# joint analysis
dat <- data.frame(class=tenyears.ind, x = cbind(t(air),t(mslp),t(prate),t(shum)))

dat.train <- dat[1:nrow(dat)%%2 == 0,]
dat.test <- dat[1:nrow(dat)%%2 == 1,]
mRF <- ranger::ranger(class~., data = dat.train, probability = TRUE)

ordering.array <- array(dim=c(4,4,nrow(dat.test)))
labels <- dat.test$class-1


for (i in 1:4) {
  for (j in 1:4) {
    ordering.array[i,j,] <- predict(mRF, data = dat.test)$predictions[,j]-predict(mRF, data = dat.test)$predictions[,i]
  }
}

# get the tv mat
tv.mat <- getTvLbDistanceMatrix(labels = labels, ordering.array = ordering.array)
tv.mat
