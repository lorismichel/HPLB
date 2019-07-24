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
plot(power.table$tvhat_search_rf, pch=19, xlab="p", ylab="tv", ylim=c(-0,0.4))
points(pmax(0, power.table$tvhat_search_mmd), pch=19,col="red")
points(pmax(0, power.table$tvhat_binomial_rf), pch=19,col="blue")
dev.off()
