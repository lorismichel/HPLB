# POSTPROCESS 2
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_2.Rdata"))

# libs
require(data.table)


colfunc <- colorRampPalette(c("black", "white"))
cols = rev(colfunc(length(seq(-0.1,1.1,0.1))))
print(cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1)))
		 
# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_2.png"))
plot(power.table[n == 10000, ]$p, power.table[n==10000,]$tvhat,pch=19, xlab="p", ylab="tv")
dev.off()
