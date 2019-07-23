# POSTPROCESS 4
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_4.Rdata"))

# libs
require(data.table)


# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_4.png"))
par(mfrow=c(1,1))
plot(power.table$tvhat_search,pch=19, xlab="p", ylab="tv")
points(pmax(0, power.table$tvhat_binomial),pch=19,col="red")
dev.off()
