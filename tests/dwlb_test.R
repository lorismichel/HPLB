# simulation example
require(distrEx)
require(ranger)
require(data.table)

set.seed(0)
n <- 2000
t <- runif(n)
x <- ifelse(t <= 0.4, rnorm(n), rnorm(n, mean = 2))
plot(t,x,pch=19)

# get a prediction
rf <- ranger::ranger(formula = t~x, data = data.frame(t = t[1:(n/2)], x = x[1:(n/2)]))
preds <- predict(rf, data = data.frame(t = t[-c(1:(n/2))], x = x[-c(1:(n/2))]))$predictions

# estimate of lower-bound
s <- seq(0.1, 0.9, 0.1)


# generate the dwit
omega_left <- sapply(s, function(ss) {ifelse(ss <= 0.4, 1, 0.4/ss)})
omega_right <- sapply(s, function(ss) {ifelse(ss <= 0.4, (0.4-ss)/(1-ss), 0)})
wit_info <- lapply(1:length(s), function(i) {
              dw_sampling(x = x[-c(1:(n/2))], f = function(x) omega_left[i]*dnorm(x) + (1-omega_left[i])*dnorm(x, mean = 2),
              g = function(x) omega_right[i]*dnorm(x) + (1-omega_right[i])*dnorm(x, mean = 2))
              })
ss <- 0.4
i <- match(ss, s)
plot(t,x,col=ifelse(wit_info[[i]]$wF==1 & t <= 0.4, "red", ifelse(wit_info[[i]]$wG==1 & t > 0.4, "blue", "black")),pch=19)

# total variation distance between the two data
tv_info <- sapply(1:length(s), function(i) {
  TotalVarDist(e1 = if(omega_left[i] == 1) Norm(0, 1) else UnivarMixingDistribution(Norm(0,1),Norm(2,1), mixCoeff = c(omega_left[i],1-omega_left[i])),
               e2 = if(omega_right[i] == 0) Norm(2, 1) else UnivarMixingDistribution(Norm(0,1),Norm(2,1), mixCoeff = c(omega_right[i],1-omega_right[i])))
})

# fit dwlb
res <- dwlb_min(times = t[-c(1:(n/2))], preds = preds, s = s, verbose.plot = FALSE)

# witness lower bounds
plot(s,unlist(lapply(wit_info, function(l) sum(l$wF))),type="b",col="red",ylim=c(0, n/2))
lines(s, res$lambdahat_asymptoticFs,type="b")

plot(s,unlist(lapply(wit_info, function(l) sum(l$wG))),type="b",col="red",ylim=c(0,n/2))
lines(s, res$lambdahat_asymptoticGs,type="b")

# tv lower bounds
plot(s, tv_info, type="b",col="red",ylim=c(0,1))
lines(s, res$TVhat_asymptoticFs,type="b")
