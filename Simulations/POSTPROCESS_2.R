# POSTPROCESS 2
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_2.Rdata"))

# libs
require(data.table)


print(summary(power.table[n==10000,]$tvhat_binomial))		 
# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_2.png"))
par(mfrow=c(2,2))
plot(power.table[m == 0.5, ]$p, power.table[m==0.5,]$tvhat_search,pch=19, xlab="p", ylab="tv")
points(power.table[m == 0.5, ]$p, pmax(0, power.table[m==0.5,]$tvhat_binomial),pch=19,col="red")

plot(power.table[m == 2, ]$p, power.table[m==2,]$tvhat_search,pch=19, xlab="p", ylab="tv")
points(power.table[m == 2, ]$p, pmax(0, power.table[m==2,]$tvhat_binomial),pch=19,col="red")

plot(power.table[m == 4, ]$p, power.table[m==4,]$tvhat_search,pch=19, xlab="p", ylab="tv")
points(power.table[m == 4, ]$p, pmax(0, power.table[m==4,]$tvhat_binomial),pch=19,col="red")

plot(power.table[m == 10, ]$p, power.table[m==10,]$tvhat_search,pch=19, xlab="p", ylab="tv")
points(power.table[m == 10, ]$p, pmax(0, power.table[m==10,]$tvhat_binomial),pch=19,col="red")
dev.off()
