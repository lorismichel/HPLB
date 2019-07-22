# POSTPROCESS 1
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_1.Rdata"))

# libs
require(data.table)


# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_1.png"))
plot(power.table$logn,power.table$loglambda,
     cex = power.table$power_search*2, 
     pch=19,xlab="log(n)",
     ylab="-gamma*log(n)")
points(power.table$logn,power.table$loglambda, 
       cex = power.table$power_binomial*2, pch=19,col="red")
abline(a = 0, b = -1,col="blue")
dev.off()
