# SIMULATION_6
# type: checking use of estimator over several splits.

## saving paths
PATH.SAVE = "../Data/"

## requirements
require(dWit)
require(parallel)
require(data.table)
require(stats)
require(ranger)
require(MASS)
require(mlbench)
require(distrEx)





# parameters for the simulation setting
m1 <- 0
m2 <- 6
m3 <- 6
m4 <- 0

startpoints<-c(0, 1/3,1/2, 2/3)
sigma1 <- 1
sigma2 <- 1
sigma3 <- 1
sigma4 <- 1

# mixing weights
# From the left:
# Before entering the interval: 0
# If we are in the interval of the class: (s- starting point)/s (where starting point can be 0)!
# After leaving the interval: (length of interval)/s


mix_left<- function(s, startpoints){

  nrclasses<-length(startpoints)
  startpointss<-data.frame(val=c(startpoints,s), label=c(1:length(startpoints),-100)) # c(paste(c(1:length(startpoints))),"s")
  startpointss<-startpointss[order(startpointss$val),]

  rowofs<-which(startpointss$label==-100)
  group<-startpointss[rowofs-1,"label"]
  sp <-startpointss[rowofs-1,"val"]

  weights<-rep(0,nrclasses)

  #After leaving the interval: (length of interval)/s
  if (group == 1){
    weights[1]=1
  }else{
    for (i in 1:(group-1)){weights[i]<-(startpoints[i+1] -startpoints[i])/s}
  }

  #If we are in the interval of the class: (s- starting point)/s (where starting point can be 0)!
  weights[group]=(s- sp)/s

  #Before entering the interval: 0 (as initialized)


  return(weights)
}

mix_right<-function(s, startpoints){
  weights<-rev(mix_left(1-s,head(cumsum(c(0,rev(diff(c(startpoints,1))))),-1)))

  #tmp<-data.frame(weights=weights, level=1:length(weights))
  #weights<-tmp[order(-tmp$level),"weights"]

  return(weights)
}

#mix_left <- function(s) { if (s<=1/3) c(1,0,0,0) else if (s<=1/3+1/6) c((1/3)/s,(s-(1/3))/s, 0,0) else if (s<=2/3) c((1/3)/s, (s-(1/3))/s, (s-(1/3 + 1/6))/s ,0) else c((1/3)/s, (1/3)/s, 1-((2/3)/s))}
#mix_right <- function(s) { if (s<=1/3) c(((1/3)-s)/(1-s),(1/3)/(1-s),(1/3)/(1-s)) else if (s<=2/3) c(0,(2/3-s)/(1-s),(1/3)/(1-s)) else c(0,0,1)}


# densities of F_s and G_s
f_s <- function(s,startpoints) { return(function(x)   sapply(x, function(xx) sum(mix_left(s,startpoints)*c(dnorm(xx, m1,  sd=sigma1),dnorm(xx, m2,sd=sigma2),dnorm(xx, m3,sd=sigma3), dnorm(xx, m4,sd=sigma4)))))}
g_s <- function(s,startpoints) { return(function(x)   sapply(x, function(xx) sum(mix_right(s,startpoints)*c(dnorm(xx, m1,sd=sigma1),dnorm(xx, m2,sd=sigma2),dnorm(xx, m3,sd=sigma3), dnorm(xx, m4,sd=sigma4)))))}


# Distributions as implemented in distrEX
F_s <- function(s, startpoints) { if (s <= startpoints[2]) {Norm(m1,sigma1)} else {UnivarMixingDistribution(Norm(m1,sigma1), Norm(m2,sigma2),Norm(m3,sigma3),Norm(m4,sigma4),
                                                                                                            mixCoeff= mix_left(s,startpoints))}}
G_s <- function(s, startpoints) { if (s >= startpoints[length(startpoints)]) {Norm(m4,sigma4)} else {UnivarMixingDistribution(Norm(m1,sigma1), Norm(m2,sigma2),Norm(m3,sigma3),Norm(m4,sigma4),
                                                                                                                              mixCoeff= mix_right(s,startpoints))}}
# sampling the data with the witnesses as well as the split s
set.seed(0)
#par(mfrow=c(1,1))
n=2000
t <- runif(n = n, min = 0, max = 1)

#for (j in 1:lenght(startpoints)){
#    eval(as.name(paste0("m",j)))
#}

#x <- ifelse(t<=startpoints[2], rnorm(sum(t<=startpoints[2]), mean = m1, sd = sigma1), ifelse(startpoints[2]<t & t<=startpoints[3], rnorm(sum(startpoints[2]<t & t<=startpoints[3]), mean = m2, sd = sigma2),ifelse(startpoints[3]<t & t<=startpoints[4], rnorm(sum(startpoints[3]<t & t<=startpoints[4]), mean = m3, sd = sigma3),rnorm(sum(startpoints[4]<t & t<=1), mean = m4, sd = sigma4))))
x<-c()
x[t<=startpoints[2]] <- rnorm(sum(t<=startpoints[2]), mean = m1, sd = sigma1)
x[startpoints[2]<t & t<=startpoints[3]] <- rnorm(sum(startpoints[2]<t & t<=startpoints[3]), mean = m2, sd = sigma2)
x[startpoints[3]<t & t<=startpoints[4]] <- rnorm(sum(startpoints[3]<t & t<=startpoints[4]), mean = m3, sd = sigma3)
x[startpoints[4]<t & t<=1] <- rnorm(sum(startpoints[4]<t & t<=1), mean = m4, sd = sigma4)




# combined data
plot(t,x,pch=19, font.main=3,font.lab=3,main="Univariate Gaussian Change in Mean", xlab="time", ylab="x")

# sample splitting between test and train
id <- sample(1:length(t), size = length(t)/2, replace = F)

# combined data
d <- data.frame(x,t)
dtrain <-d[id,]
dtest <- d[-id,]



# total variation and densities
set.seed(0)
s_vec <- runif(100)
tv_s <- c()
for (s in s_vec) {
  tv_s <- c(tv_s, distrEx::TotalVarDist(F_s(s,startpoints),G_s(s,startpoints)))
}
plot(sort(s_vec), tv_s[order(s_vec)],type="l", main = "total variation distance")

for (j in startpoints){
  abline(v=j, col="red", lty=2)

}




## evaluation quickly at multiple splits the lower-bound estimate for TV
set.seed(0)
s_vec <- seq(0.1,0.9,by=0.05) #runif(20, min=0.1, max=0.9)
tv_s <- c()
tvhat_s <- c()
withat_s <- c()
for (s in s_vec) {

  tv_s <- c(tv_s, distrEx::TotalVarDist(F_s(s, startpoints),G_s(s,startpoints)))

}


type="regression"
#type="classification"


if (type=="regression"){
  rf <- ranger::ranger(t~., data = dtrain)
  # getting the ranks of the predictions on the test set
  preds.test <- predict(rf, data = dtest)$predictions

  s <- s_vec
  dsearchTV <- dWit(t = dtest$t, rho = preds.test, s = s, alpha = 0.05, estimator.type = "tv-search")
  dsearchTV$tvhat

}




if (type=="classification"){

  tvhat<-c()
  l=0
  preds.test <- c()
  for (s in s_vec){

    l<-l+1
    dtrains<-dtrain
    dtrains$t<-NULL
    dtrains$ytrain <- as.factor(as.numeric(dtrain$t > s))

  rf <- ranger::ranger(ytrain~., data = dtrains, probability = T)
  # getting the ranks of the predictions on the test set
  preds.test <- cbind(preds.test, predict(rf, data = dtest)$predictions[,"1"])

  #ytest <- as.numeric(dtest$t <= s)

  }

  dsearchTV <- dWit(t = dtest$t, rho = preds.test, s = s_vec, alpha = 0.05, estimator.type = "tv-search")


}



N=1
ind1 <- (s_vec <= 0.5*N)
ind2 <- ((rep(1,length(s)) - ind1)==1)
tvhat_s<-rep(0,length(s_vec))
tvhat_s[ind1]<-l$TVhat_asymptoticFs[ind1]
tvhat_s[ind2]<-l$TVhat_asymptoticGs[ind2]

withat_s[ind1]<-l$lambdahat_asymptoticFs[ind1]
withat_s[ind2]<-l$lambdahat_asymptoticGs[ind2]


# ploting the lower bound estimates as well as the true ones
plot(sort(s_vec), tv_s[order(s_vec)],type="b", main = "Total Variation Distance",ylim=c(0,1),pch=19, xlab="s", ylab="")
lines(x = sort(s_vec), y = tvhat_s[order(s_vec)], col="blue",type="b",pch=19)

# Plotting the behavior of the cleaned Vtildes!

#Vtildes<- sapply(s_vec, function(s) {wits<-witness_pointing(x = x, f = f_s(s, startpoints), g = g_s(s, startpoints)); return(cleaning(x, t,wits$wF, wits$wG , s)$Vtildezmax)   })
#plot(sort(s_vec),Vtildes[order(s_vec)], type="l")


# witnesses lower bound
wit_s <- sapply(s_vec, function(s) {ifelse(s <= 0.5, sum(witness_pointing(seed = NULL, x = x, f = f_s(s, startpoints), g = g_s(s, startpoints))$wF[-id][dtest$t <= s]), sum(witness_pointing(x = x, f = f_s(s, startpoints), g = g_s(s, startpoints))$wG[-id][dtest$t > s]))})
plot(sort(s_vec), wit_s[order(s_vec)],type="b", main = "Witness Count",pch=19, xlab="s", ylab="")
lines(x = sort(s_vec), y = withat_s[order(s_vec)], col="blue",type="b",pch=19)





