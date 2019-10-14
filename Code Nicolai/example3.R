library(mlbench)

library(randomForest)
source("lambdahat_code.R")

data(BostonHousing)
Y = as.factor( BostonHousing$crim > median(BostonHousing$crim))
X = BostonHousing[,-1]
for (k in 1:ncol(X))X[,k] = as.numeric(X[,k])

n <- nrow(X)

train <- 1:round(n/2)
test <- (1:n)[-train]
rf <- randomForest(X[train,],Y[train])
print(rf)

pred <- predict(rf, X[test,],type="prob")
conf = table( Y[test], predict(rf))

tvhyper = gethypertv(conf)
            
tvhat = dWit( t=as.numeric(Y[test]), rho=as.numeric(pred[,2]-pred[,1]),s =1.5, tv.seq   = sort(unique(c(seq(0,1,by=0.01),c(0,0.0001,0.001,0.002,0.005,0.01,0.015,0.02)))))$tvhat




gethypertv <- function(tab,alpha=0.05){

    nsim = 200
    stat = sum(diag(tab))
    statsim = numeric(nsim)
    indr = c( rep(0, sum(tab[1,])), rep(1, sum(tab[2,])))
    indc = c( rep(0, sum(tab[,1])), rep(1, sum(tab[,2])))
    for (sim in 1:nsim){
        sampr = sample( indr, length(indr))
        sampc = sample(indc, length(indc))
        tabsim = table(sampr, sampc)
        statsim[sim] = sum(diag(tabsim))
    }
    ## still needs binomial adjustment to go from witness bound to TV bound
    return( max(0, (stat - quantile( statsim, 1-alpha))  /sum(tab)))
    
}
