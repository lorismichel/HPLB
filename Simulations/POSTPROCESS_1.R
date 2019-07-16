PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "powerStudy.Rdata"))

require(data.table)

png(filename = paste0(PATH.PLOTS,"sim1.png"))
plot(power.table$logn,power.table$loglambda,cex = power.table$power_search*2, pch=19,xlab="log(n)",ylab="-gamma*log(n)")
points(power.table$logn,power.table$loglambda,cex = power.table$power_binomial*2, pch=19,col="red")
abline(a = 0, b = -1,col="blue")
dev.off()
