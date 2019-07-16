PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"
load(paste0(PATH.DATA, "powerStudy.Rdata"))

png(filename = paste0(PATH.PLOTS,"sim1.png"))
plot(power.table$logn,power.table$loglambda,cex = power.table$power*2, pch=19,xlab="log(n)",ylab="-gamma*log(n)")
abline(a = 0, b = -1,col="blue")
dev.off()
