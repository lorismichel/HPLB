# climate continuous change point detection

# options
RANK.TRANSFORMATION <- FALSE
SPLIT.TRAIN.TEST <- 1

# source
source("./Climate/climate-preprocessing.R")

# repro
set.seed(1)

# libs
require(dWit)

# prepro
d <- climatePrePro()

# zurich coordinates
zurich.coord <- c(8.5391825, 47.3686498)
air   <- extract(d$air_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]
shum  <- extract(d$shum_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]
prate <- extract(d$prate_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]
mslp  <- extract(d$mslp_raster, matrix(zurich.coord,ncol=2),method="bilinear")[1,1:14641]


if (RANK.TRANSFORMATION) {
  air   <- rank(air)
  shum  <- rank(shum)
  prate <- rank(prate)
  mslp  <- rank(mslp)
}

# analysis
dat <- data.frame(class=d$tenyears.ind,
                  x = cbind(t(air),t(mslp),t(prate),t(shum)))

# air time series

# PLOT_CLIMATE_MIXTURE_1.png
par(mfrow=c(2,2))
plot(as.Date(d$time), air,type="l",xlab='Time',ylab='Temperature', main="Zürich temperature",font.lab = 1,font.main=1)
abline(v = as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))],lty=2,col='blue')

plot(as.Date(d$time), mslp,type="l",xlab='Time',ylab='Pressure', main="Zürich pressure",font.lab = 1,font.main=1)
abline(v = as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))],lty=2,col='blue')

plot(as.Date(d$time), prate,type="l",xlab='Time',ylab='Precipitation', main="Zürich precipitation",font.lab = 1,font.main=1)
abline(v = as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))],lty=2,col='blue',font.lab = 1,font.main=1)

plot(as.Date(d$time), shum,type="l",xlab='Time',ylab='Humidity', main="Zürich humidity",font.lab = 1,font.main=1)
abline(v = as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))],lty=2,col='blue')


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


# fit a forest
mRF_air <- ranger::ranger(t~., data = data.frame(t=1:length(air), x=air)[ind.train,])
mRF_mslp <- ranger::ranger(t~., data = data.frame(t=1:length(mslp), x=mslp)[ind.train,])
mRF_prate <- ranger::ranger(t~., data = data.frame(t=1:length(prate), x=prate)[ind.train,])
mRF_shum <- ranger::ranger(t~., data = data.frame(t=1:length(shum), x=shum)[ind.train,])

mRF_joint <- ranger::ranger(t~., data = data.frame(t=1:length(shum), air=air, mslp=mslp, prate=prate, shum=shum)[ind.train,], importance = 'permutation')

# run dWit
require(dWit)

# marginals
# PLOT_CLIMATE_MIXTURE_2.png
par(mfrow=c(2,2))
plot(as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))], dWit(t = c(1:length(air))[ind.test], rho = predict(mRF_air, data.frame(x=air[ind.test]))$predictions,
     s = as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10))), estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[H]),type="b",pch=19,xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Zürich temperature")


plot(as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))],dWit(t = c(1:length(mslp))[ind.test], rho = predict(mRF_mslp, data.frame(x=mslp[ind.test]))$predictions,
     s = as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10))), estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[H]),type="b",pch=19,xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Zürich pressure")


plot(as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))], dWit(t = c(1:length(prate))[ind.test], rho = predict(mRF_prate, data.frame(x=prate[ind.test]))$predictions,
     s = as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10))), estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[H]),type="b",pch=19,xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Zürich precipitation")


plot(as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))], dWit(t = c(1:length(shum))[ind.test], rho = predict(mRF_shum, data.frame(x=shum[ind.test]))$predictions,
     s = as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10))), estimator.type = "asymptotic-tv-search")$tvhat,ylab=expression(hat(lambda)[H]),type="b",pch=19,xlab="Time",ylim=c(0,0.1),font.lab = 1,font.main=1,main="Zürich humidity")


# PLOT_CLIMATE_MIXTURE_3.png
par(mfrow=c(1,1))
plot(as.Date(d$time)[as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))], dWit(t = c(1:length(air))[ind.test], rho = predict(mRF_joint, data.frame(air=air, mslp=mslp, prate=prate, shum=shum)[ind.test,])$predictions,
     s = as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10))), estimator.type = "asymptotic-tv-search")$tvhat,type="b",pch=19,ylim=c(0,0.1),font.lab = 1,font.main=1,xlab="Time",ylab=expression(hat(lambda)[H]),main="Zürich (all signals)")

# PLOT_CLIMATE_MIXTURE_4.png
par(mfrow=c(2,2))
plot(density(air[(1:length(shum)>7971.667)]),col="grey",xlab='Temperature',ylab='Density', main="Zürich temperature (marginal)",font.lab = 1,font.main=1)
lines(density(air[(1:length(shum)<=7971.667)]),col="darkblue")

plot(density(mslp[(1:length(shum)>7971.667)]),col="grey",xlab='Pressure',ylab='Density', main="Zürich pressure (marginal)",font.lab = 1,font.main=1)
lines(density(mslp[(1:length(shum)<=7971.667)]),col="darkblue")

plot(density(prate[(1:length(shum)>7971.667)]),col="grey",xlab='',ylab='Precipitation', main="Zürich precipitation (marginal)",font.lab = 1,font.main=1)
lines(density(prate[(1:length(shum)<=7971.667)]),col="darkblue")

plot(density(shum[(1:length(shum)>7971.667)]),col="grey",xlab='Humidity',ylab='Density', main="Zürich humidity (marginal)",font.lab = 1,font.main=1)
lines(density(shum[(1:length(shum)<=7971.667)]),col="darkblue")


# investigation plot
par(mfrow=c(4,3))
for (s in as.numeric(quantile(1:length(air), seq(0.1,0.9,length.out = 10)))) {
  plot(density(shum[(1:length(shum)>s)]),col="grey",xlab='Humidity',ylab='Density', main="Zürich humidity (marginal)",font.lab = 1,font.main=1)
  lines(density(shum[(1:length(shum)<=s)]),col="darkblue")
}

# PLOT_CLIMATE_MIXTURE_5.png
par(mfrow=c(1,1))
pairs(cbind(temperature=air,pressure=mslp,precipitation=prate,humidity=shum), col=c('darkblue', 'grey')[(1:length(shum)>7971.667)+1],pch=19,cex=0.5,font.lab = 1,font.main=1)

