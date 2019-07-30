# POSTPROCESS 2
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_3.Rdata"))

# libs
require(data.table)


# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_2.png"))
par(mfrow=c(2,2))
plot(power.table[m==0.5,]$tvhat_search,pch=19, xlab="p", ylab="TV", main="Null hypothesis")
points(pmax(0, power.table[m==0.5,]$tvhat_binomial),pch=19,col="red")

plot(power.table[m==2,]$tvhat_search,pch=19, xlab="p", ylab="TV",  main="Null hypothesis")
points(pmax(0, power.table[m==2,]$tvhat_binomial),pch=19,col="red")

plot(power.table[m==4,]$tvhat_search,pch=19, xlab="p", ylab="TV",  main="Null hypothesis")
points(pmax(0, power.table[m==4,]$tvhat_binomial),pch=19,col="red")

plot(power.table[m==10,]$tvhat_search,pch=19, xlab="p", ylab="TV",  main="Null hypothesis")
points(pmax(0, power.table[m==10,]$tvhat_binomial),pch=19,col="red")
dev.off()
