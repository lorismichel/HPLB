# POSTPROCESS 5
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_5.Rdata"))

# libs
require(data.table)


# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_5.png"))
par(mfrow=c(1,2))
# Boston
plot(power.table.boston$p, power.table$tvhat_search, pch=19, xlab="p", ylab="tv",main="Boston")
points(power.table.boston$p,pmax(0, power.table$tvhat_binomial), pch=19,col="red")
# titanic
plot(power.table.titanic$p, power.table$tvhat_search, pch=19, xlab="p", ylab="tv",main="titanic")
points(power.table.titanic$p,pmax(0, power.table$tvhat_binomial), pch=19,col="red")
dev.off()
