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
par(mfrow=c(2,2))
# Boston
plot(power.table.boston$p, power.table.boston$tvhat_search, pch=19, xlab="p", ylab="tv",main="Boston dataset",type="b")
points(power.table.boston$p,pmax(0, power.table.boston$tvhat_binomial), pch=19,col="red",type="b")
# titanic
plot(power.table.titanic$p, power.table.titanic$tvhat_search, pch=19, xlab="p", ylab="tv",main="Titanic dataset", type="b")
points(power.table.titanic$p,pmax(0, power.table.titanic$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.breast$p, power.table.breast$tvhat_search, pch=19, xlab="p", ylab="tv",main="Breast cancer dataset", type="b")
points(power.table.breast$p,pmax(0, power.table.breast$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.iono$p, power.table.iono$tvhat_search, pch=19, xlab="p", ylab="tv",main="Ionosphere dataset", type="b")
points(power.table.iono$p,pmax(0, power.table.iono$tvhat_binomial), pch=19,col="red", type="b")
dev.off()
