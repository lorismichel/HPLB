# POSTPROCESS 2
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_2.Rdata"))

# libs
require(data.table)


		 
# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_2.png"))
plot(power.table[n == 10000, ]$p, power.table[n==10000,]$tvhat_search,pch=19, xlab="p", ylab="tv")
points(power.table[n == 10000, ]$p, power.table[n==10000,]$tvhat_binomial,pch=19,col="red")
dev.off()
