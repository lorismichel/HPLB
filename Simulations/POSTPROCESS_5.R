# POSTPROCESS 5
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "DATA_SIMULATION_5.Rdata"))
#load("C:/Users/jeffr/Dropbox/tvForest/SimulationsData/DATA_SIMULATION_5.Rdata")

# libs
require(data.table)


# loglogplot
png(filename = paste0(PATH.PLOTS,"PLOT_SIMULATION_5.png"))
par(mfrow=c(3,3))

# Boston
plot(power.table.boston$p, power.table.boston$tvhat_search, ylim=c(0,1),pch=19, xlab="p", ylab="TV",main="Boston",type="b", cex = 0.7, font.main=3)
points(power.table.boston$p,pmax(0, power.table.boston$tvhat_binomial), pch=19,col="red",type="b")

# titanic
plot(power.table.titanic$p, power.table.titanic$tvhat_search, pch=19, ylim=c(0,1), xlab="p", ylab="TV",main="Titanic", type="b", cex = 0.7, font.main=3)
points(power.table.titanic$p,pmax(0, power.table.titanic$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.breast$p, power.table.breast$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Breast Cancer", type="b", cex = 0.7, font.main=3)
points(power.table.breast$p,pmax(0, power.table.breast$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.iono$p, power.table.iono$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Ionosphere", type="b", cex = 0.7, font.main=3)
points(power.table.iono$p,pmax(0, power.table.iono$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.abalone$p, power.table.abalone$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Abalone", type="b", cex = 0.7, font.main=3)
points(power.table.abalone$p,pmax(0, power.table.abalone$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.adult$p, power.table.adult$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Adult", type="b", cex = 0.7, font.main=3)
points(power.table.adult$p,pmax(0, power.table.adult$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.banknotes$p, power.table.banknotes$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Bank Notes", type="b", cex = 0.7, font.main=3)
points(power.table.banknotes$p,pmax(0, power.table.banknotes$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.default$p, power.table.default$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="Default (ISLR)", type="b", cex = 0.7, font.main=3)
points(power.table.default$p,pmax(0, power.table.default$tvhat_binomial), pch=19,col="red", type="b")

plot(power.table.credit$p, power.table.credit$tvhat_search, ylim=c(0,1), pch=19, xlab="p", ylab="TV",main="German Credit Default", type="b", cex = 0.7, font.main=3)
points(power.table.credit$p,pmax(0, power.table.credit$tvhat_binomial), pch=19,col="red", type="b")
dev.off()
