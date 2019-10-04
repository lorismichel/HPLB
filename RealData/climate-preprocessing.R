# climatic data (real analysis)

# libs
library(sf)
library(ncdf4)
library(raster)
library(rasterVis)
library(RColorBrewer)


climatePrePro <- function(path = "~/Downloads/reanalysis_jeff_loris/") {
# loading the data (4 signals)
air_raster <- brick(paste(path, "NCEP2.air.2m.day", ".nc", sep=""), varname="air")
mslp_raster <- brick(paste(path, "NCEP2.mslp.2m.day", ".nc", sep=""), varname="mslp")
prate_raster <- brick(paste(path, "NCEP2.prate.2m.day", ".nc", sep=""), varname="prate")
shum_raster <- brick(paste(path, "NCEP2.shum.2m.day", ".nc", sep=""), varname="shum")



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
return(list(time = time, tenyears.ind = tenyears.ind, air = air, mslp = mslp, prate = prate, shum = shum,
            air_raster = air_raster, mslp_raster = mslp_raster, prate_raster = prate_raster, shum_raster = shum_raster))
}

