# POSTPROCESS 1
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_1.Rdata"))

# libs
require(data.table)


colfunc <- colorRampPalette(c("black", "white"))
cols = rev(colfunc(length(seq(-0.1,1.1,0.1))))
print(cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1)))
		 
# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_1.png"))
par(mfrow=c(1,2))
plot(power.table$logn,power.table$loglambda,
     col = cols[cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1))],
     pch=19,xlab="log(n)",
     ylab="-gamma*log(n)")
plot(power.table$logn,power.table$loglambda,
     col = cols[cut(power.table$power_binomial, labels = FALSE, breaks = seq(-0.1,1.1,0.1))],
     pch=19,xlab="log(n)",
     ylab="-gamma*log(n)")
dev.off()
