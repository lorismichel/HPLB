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


res[,tvhat_search := mean(tvhat_search,na.rm=T),by=c("p","dataset")]
res[,tvhat_binomial := mean(tvhat_binomial,na.rm=T),by=c("p","dataset")]

dataset.names = unique(res$dataset)

for (n in dataset.names) {

 plot(res[dataset==n,]$p, res[dataset==n,]$tvhat_search, ylim=c(0,1),pch=19, xlab="p", ylab="TV",main=n,type="l", cex = 0.7, font.main=3)
 points(res[dataset==n,]$p,pmax(0, res[dataset==n,]$tvhat_binomial), pch=19,col="red",type="l")

}
# aggregate the res for same p
#res[,tvhat_search := mean(tvhat_search,na.rm=T),by=c("p","dataset")]
#res[,tvhat_binomial := mean(tvhat_binomial,na.rm=T),by=c("p","dataset")]

# Boston
dev.off()
