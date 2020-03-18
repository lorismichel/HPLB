# climate continuous change point detection

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

# loc coordinates and data
loc.coord <- c(-45,-8) #c(8.5391825, 47.3686498)
air   <- extract(d$air_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]
shum  <- extract(d$shum_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]
prate <- extract(d$prate_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]
mslp  <- extract(d$mslp_raster, matrix(loc.coord,ncol=2),method="bilinear")[1,1:14641]


# time
time <- as.Date(d$time)

# split ids (plots)
split.ids   <- as.numeric(quantile(1:length(air), seq(0,1,length.out = 5)[-c(1,5)]))
split.dates <- time[split.ids]

# PLOT_CLIMATE_MIXTURE_1.png
png(filename = paste0("./Plots/PLOT_CLIMATE_MIXTURE_RAW_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",
                      PREPRO,
                      "_PLOT_1.png"),
    width = 1000, height = 700)
par(mar=c(5.4,4.5,4.5,2.3))
par(mfrow=c(2,2))
plot(as.Date(d$time), air,type="l",xlab='Time',ylab='Temperature',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue')

plot(as.Date(d$time), mslp,type="l",xlab='Time',ylab='Pressure',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue')

plot(as.Date(d$time), prate,type="l",xlab='Time',ylab='Precipitation',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue',font.lab = 1,font.main=1)

plot(as.Date(d$time), shum,type="l",xlab='Time',ylab='Humidity',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue')
dev.off()

# preprocessing
if (PREPRO == 1) {
  air   <- rank(air)
  shum  <- rank(shum)
  prate <- rank(prate)
  mslp  <- rank(mslp)
} else if (PREPRO == 2) {
  time <- time[-1]
  air   <- diff(air)
  shum  <- diff(shum)
  prate <- diff(prate)
  mslp  <- diff(mslp)
} else if (PREPRO == 3) {
  time <- time[-1]
  air   <- rank(diff(air))
  shum  <- rank(diff(shum))
  prate <- rank(diff(prate))
  mslp  <- rank(diff(mslp))
}


# Plot preprocessed data

if (PREPRO == 2) {
# PLOT_CLIMATE_MIXTURE_1.png
png(filename = paste0("./Plots/PLOT_CLIMATE_MIXTURE_TRANS_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",
                      PREPRO,
                      "_PLOT_1.png"),
    width = 1000, height = 700)
par(mar=c(5.4,4.5,4.5,2.3))
par(mfrow=c(2,2))
par
plot(time, air,type="l",xlab='Time',ylab='Diff. temperature',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue')

plot(time, mslp,type="l",xlab='Time',ylab='Diff. pressure',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue')

plot(time, prate,type="l",xlab='Time',ylab='Diff. precipitation',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue',font.lab = 1,font.main=1)

plot(time, shum,type="l",xlab='Time',ylab='Diff. humidity',font.lab = 1,font.main=1, cex.axis = 2, cex.lab = 2)
abline(v = split.dates,lty=2,col='blue')
dev.off()
}

# run analysis

# split ids (analysis)
split.ids   <- as.numeric(quantile(1:length(air), seq(0,1,length.out = 9)[-c(1,9)]))
split.dates <- time[split.ids]

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

# create train and test sets
train <- data.frame(t=1:length(air), air = air, mslp = mslp, shum = shum, prate = prate)[ind.train,]
test  <- data.frame(t=1:length(air), air = air, mslp = mslp, shum = shum, prate = prate)[ind.test,]

# fit the different forests
mRF_air   <- ranger::ranger(t~air, data = train)
mRF_mslp  <- ranger::ranger(t~mslp, data = train)
mRF_prate <- ranger::ranger(t~prate, data = train)
mRF_shum  <- ranger::ranger(t~shum, data = train)
mRF_joint <- ranger::ranger(t~., data = train,
                            importance = 'permutation')


# marginals
# PLOT_CLIMATE_MIXTURE_2.png
png(filename = paste0("./Plots/PLOT_CLIMATE_MIXTURE_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",
                      PREPRO,
                      "_PLOT_2.png"),
    width = 1000)

par(mfrow=c(2,2))
plot(split.dates,
     dWit(t = c(1:length(air))[ind.test], rho = predict(mRF_air, test)$predictions,
     s = split.ids,
     estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[abc]),type="b",pch=19,
     xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Diff. temperature",cex.lab=2.5,cex.main=2.5,cex.axis=2)
plot(split.dates,
     dWit(t = c(1:length(mslp))[ind.test], rho = predict(mRF_mslp, test)$predictions,
     s = split.ids,
     estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[abc]),type="b",pch=19,
     xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Diff. pressure",cex.lab=2.5,cex.main=2.5,cex.axis=2)
plot(split.dates,
     dWit(t = c(1:length(prate))[ind.test], rho = predict(mRF_prate, test)$predictions,
     s = split.ids,
     estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[abc]),type="b",pch=19,
     xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Diff. precipitation",cex.lab=2.5,cex.main=2.5,cex.axis=2)
plot(split.dates,
     dWit(t = c(1:length(shum))[ind.test], rho = predict(mRF_shum, test)$predictions,
     s = split.ids,
     estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[abc]),type="b",pch=19,
     xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Diff. humidity",cex.lab=2.5,cex.main=2.5,cex.axis=2)
dev.off()

# joint
# PLOT_CLIMATE_MIXTURE_3.png
png(filename = paste0("./Plots/PLOT_CLIMATE_MIXTURE_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",
                      PREPRO,"_PLOT_3.png"),
    width = 1000)

par(mfrow=c(1,1))
plot(split.dates,
     dWit(t = c(1:length(air))[ind.test], rho = predict(mRF_joint, test)$predictions,
     s = split.ids, estimator.type = "asymptotic-tv-search")$tvhat,
     type="b",pch=19,ylim=c(0,0.2),font.lab = 1,font.main=1,xlab="Time",
     ylab=expression(hat(lambda)[H]),main="All differenced variables",cex.lab=2.5,cex.main=2.5,cex.axis=2)
dev.off()

# PLOT_CLIMATE_MIXTURE_4.png
png(filename = paste0("./Plots/PLOT_CLIMATE_MIXTURE_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",PREPRO,
                      "_PLOT_4.png"),
    width = 1000)

par(mfrow=c(2,2))
plot(density(air[(1:length(shum)>7320.50)]),lwd=2,col="grey",xlab='Temperature',
     ylab='Density', main="Temperature (marginal)",font.lab = 1,font.main=1,cex.lab=2.5,cex.main=2.5,cex.axis=2,ylim=c(0,0.9))
lines(density(air[(1:length(shum)<=7320.50)]),col="darkblue")

plot(density(mslp[(1:length(shum)>7320.50)]),lwd=2,col="grey",xlab='Pressure',
     ylab='Density', main="Pressure (marginal)",font.lab = 1,font.main=1,cex.lab=2.5,cex.main=2.5,cex.axis=2)
lines(density(mslp[(1:length(shum)<=7320.50)]),col="darkblue")

plot(density(prate[(1:length(shum)>7320.50)]),lwd=2,col="grey",xlab='',
     ylab='Precipitation', main="Precipitation (marginal)",font.lab = 1,font.main=1,cex.lab=2.5,cex.main=2.5,cex.axis=2)
lines(density(prate[(1:length(shum)<=7320.50)]),col="darkblue")

plot(density(shum[(1:length(shum)>7320.50)]),lwd=2,col="grey",xlab='Humidity',
     ylab='Density', main="Humidity (marginal)",font.lab = 1,font.main=1,cex.lab=2.5,cex.main=2.5,cex.axis=2, ylim=c(0,800))
lines(density(shum[(1:length(shum)<=7320.50)]),col="darkblue")
dev.off()

# investigation plot
par(mfrow=c(4,3))
for (s in split.ids) {
  plot(density(shum[(1:length(shum)>s)]),col="grey",xlab='Humidity',
       ylab='Density', main="Humidity (marginal)",font.lab = 1,font.main=1,cex.lab=2.5,cex.main=2.5,cex.axis=2)
  lines(density(shum[(1:length(shum)<=s)]),col="darkblue")
}

# PLOT_CLIMATE_MIXTURE_5.png
png(filename = paste0("./Plots/PLOT_CLIMATE_MIXTURE_SPLIT_",
                      SPLIT.TRAIN.TEST,
                      "_PREPRO_",
                      PREPRO,
                      "_PLOT_5.png"),
    width = 1000)

par(mfrow=c(1,1))
pairs(cbind(temperature=air, pressure=mslp, precipitation=prate, humidity=shum),
      col=c('darkblue', 'grey')[(1:length(shum)>7320.50)+1],pch=19,cex=0.5,font.lab = 1,font.main=1,cex.lab=2.5,cex.main=2.5,cex.axis=2)
dev.off()

