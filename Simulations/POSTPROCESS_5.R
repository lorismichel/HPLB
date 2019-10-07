# POSTPROCESS 5
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_5.Rdata"))
#load("C:/Users/jeffr/Dropbox/tvForest/SimulationsData/DATA_SIMULATION_5.Rdata")

# libs
require(data.table)


# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_5.png"))
par(mfrow=c(3,3))

# aggregate the res for same p
res[,tvhat_search := mean(tvhat_search,na.rm=T),by=c("p","dataset")]
res[,tvhat_binomial := mean(tvhat_binomial,na.rm=T),by=c("p","dataset")]

# Boston
plot(res[dataset=="Boston",]$p, res[dataset=="Boston",]$tvhat_search, ylim=c(0,1),pch=19, xlab="p", ylab="TV",main="Boston",type="l", cex = 0.7, font.main=3)
points(res[dataset=="Boston",]$p,pmax(0, res[dataset=="Boston",]$tvhat_binomial), pch=19,col="red",type="l")

# titanic
plot(res[dataset=="titanic",]$p, res[dataset=="titanic",]$tvhat_search, pch=19, ylim=c(0,1), xlab="p", ylab="TV",main="Titanic", type="l", cex = 0.7, font.main=3)
points(res[dataset=="titanic",]$p,pmax(0, res[dataset=="titanic",]$tvhat_binomial), pch=19,col="red", type="l")

plot(res[dataset=="BreastCancer",]$p, res[dataset=="BreastCancer",]$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Breast Cancer", type="l", cex = 0.7, font.main=3)
points(res[dataset=="BreastCancer",]$p,pmax(0, res[dataset=="BreastCancer",]$tvhat_binomial), pch=19,col="red", type="l")

plot(res[dataset=="Ionosphere",]$p, res[dataset=="Ionosphere",]$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Ionosphere", type="l", cex = 0.7, font.main=3)
points(res[dataset=="Ionosphere",]$p,pmax(0, res[dataset=="Ionosphere",]$tvhat_binomial), pch=19,col="red", type="l")

plot(res[dataset=="abalone",]$p, res[dataset=="abalone",]$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Abalone", type="l", cex = 0.7, font.main=3)
points(res[dataset=="abalone",]$p,pmax(0, res[dataset=="abalone",]$tvhat_binomial), pch=19,col="red", type="l")

#plot(power.table.adult$p, power.table.adult$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Adult", type="l", cex = 0.7, font.main=3)
#points(power.table.adult$p,pmax(0, power.table.adult$tvhat_binomial), pch=19,col="red", type="l")

plot(res[dataset=="banknotes",]$p, res[dataset=="banknotes",]$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Bank Notes", type="l", cex = 0.7, font.main=3)
points(res[dataset=="banknotes",]$p,pmax(0, res[dataset=="banknotes",]$tvhat_binomial), pch=19,col="red", type="l")

plot(res[dataset=="Default",]$p, res[dataset=="Default",]$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Default (ISLR)", type="l", cex = 0.7, font.main=3)
points(res[dataset=="Default",]$p,pmax(0, res[dataset=="Default",]$tvhat_binomial), pch=19,col="red", type="l")

plot(res[dataset=="credit",]$p, res[dataset=="credit",]$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="German Credit Default", type="l", cex = 0.7, font.main=3)
points(res[dataset=="credit",]$p,pmax(0, res[dataset=="credit",]$tvhat_binomial), pch=19,col="red", type="l")
dev.off()
