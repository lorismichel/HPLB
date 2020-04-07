# POSTPROCESS 1
# type: producing plots for the loglog scatter plot of sample size against tv

# PATHS
PATH.DATA <- "../Data/"
PATH.PLOTS <- "../Plots/"

fnames.in <- c("DATA_SIMULATION_1_SC1.Rdata", "DATA_SIMULATION_1_SC2.Rdata", "DATA_SIMULATION_1_SC3.Rdata")
fnames.out <- c("PLOT_SIMULATION_1_SC1.png", "PLOT_SIMULATION_1_SC2.png", "PLOT_SIMULATION_1_SC3.png")


i<-2

load(paste0(PATH.DATA, fnames.in[i]))

power.table$power_search<-ifelse(power.table$power_search <= 0.05, 0, power.table$power_search)
power.table$power_binomial<-ifelse(power.table$power_binomial <= 0.05, 0, power.table$power_binomial)


# libs
require(data.table)


colfunc <- colorRampPalette(c("black", "white"))
cols = rev(colfunc(length(seq(-0.1,1.1,0.1))))
print(cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1)))

# loglogplot
png(filename = paste0(PATH.PLOTS,fnames.out[i]))
par(mfrow=c(2,1))
plot(power.table$logn,power.table$loglambda,
     col = cols[cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1))],
     pch=19,xlab=expression(log(N)),
     ylab=expression(-gamma %.% log(N)), main=expression(hat(lambda)[abc]),font.main=2)
abline(a = 1.3, b = -1,col="red")
#abline(a = 0, b = -1/2, col="blue")
plot(power.table$logn,power.table$loglambda,
     col = cols[cut(power.table$power_binomial, labels = FALSE, breaks = seq(-0.1,1.1,0.1))],
     pch=19,xlab=expression(log(N)),
     ylab=expression(-gamma %.% log(N)), main=expression(hat(lambda)[bc]),font.main=2)
#abline(a = 1.3, b = -1,col="red")
abline(a = 0, b = -1/2, col="blue")
#plot(seq(-0.1,1.1,0.1), rep(1, length(seq(-0.1,1.1,0.1))), col=cols,pch=19,cex=3)
dev.off()





i<-3

load(paste0(PATH.DATA, fnames.in[i]))

power.table$power_search<-ifelse(power.table$power_search <= 0.05, 0, power.table$power_search)
power.table$power_binomial<-ifelse(power.table$power_binomial <= 0.05, 0, power.table$power_binomial)


# libs
require(data.table)


colfunc <- colorRampPalette(c("black", "white"))
cols = rev(colfunc(length(seq(-0.1,1.1,0.1))))
print(cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1)))

# loglogplot
png(filename = paste0(PATH.PLOTS,fnames.out[i]))
par(mfrow=c(2,1))
plot(power.table$logn,power.table$loglambda,
     col = cols[cut(power.table$power_search, labels = FALSE, breaks = seq(-0.1,1.1,0.1))],
     pch=19,xlab=expression(log(N)),
     ylab=expression(-gamma %.% log(N)), main=expression(hat(lambda)[abc]),font.main=2)
abline(a = 1.3, b = -1,col="red")
#abline(a = 0, b = -1/2, col="blue")
plot(power.table$logn,power.table$loglambda,
     col = cols[cut(power.table$power_binomial, labels = FALSE, breaks = seq(-0.1,1.1,0.1))],
     pch=19,xlab=expression(log(N)),
     ylab=expression(-gamma %.% log(N)), main=expression(hat(lambda)[bc]),font.main=2)
abline(a = 1.3, b = -1,col="red")
#abline(a = 0, b = -1/2, col="blue")
#plot(seq(-0.1,1.1,0.1), rep(1, length(seq(-0.1,1.1,0.1))), col=cols,pch=19,cex=3)
dev.off()



