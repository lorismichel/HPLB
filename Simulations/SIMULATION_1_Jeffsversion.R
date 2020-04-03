

## requirements
require(dWit)
require(parallel)
require(data.table)

## data generation & bayes ratio
generateUnifMixturesData <- function(n, boundaries = matrix(c(-10,-9, -1, 0, 0, 1),ncol=2, byrow = TRUE), weights = c(0.1, 0.4, 0.5)) {

  force(boundaries)
  force(weights)

  # generate data
  ind.mix <- sample(1:nrow(boundaries), size = n, prob = weights, replace = TRUE)
  x <- runif(n = n, min = boundaries[ind.mix, 1], max = boundaries[ind.mix, 2])

  return(x)
}

generateUnifMixturesDensities <- function(boundaries = matrix(c(-10,-9, -1, 0, 0, 1),ncol=2, byrow = TRUE), weights = c(0.1, 0.4, 0.5)) {

  force(boundaries)
  force(weights)

  # generate density function
  d <- function(y) sum(weights * sapply(1:length(weights), function(i) dunif(x = y, min = boundaries[i,1], max = boundaries[i,2])))

  return(d)
}


nrep <- 50
grid.gamma <- rep(seq(0.3,  0.8, by = 0.1), nrep)
#grid.n <- rep(round(10^{seq(3,6,length.out = 10)}), nrep)
grid.n <- rep(round(10^{seq(3,5.5,length.out = 6)}), nrep)
grid <- expand.grid(grid.n, grid.gamma)

RUN_SC_1=T


if (RUN_SC_1) {
##### scenario 1: bad for both methods

# running simulations
set.seed(123, "L'Ecuyer")
res_tv_search_1 <- mapply(FUN = function(n, gamma) {

  print(paste(toString(n),toString(gamma) ))

  p1 <- (1/2)*n^{-gamma}
  weights <- c((1/2)-p1, (1/2)+p1)

  x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
  x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
  d.right <- generateUnifMixturesDensities(boundaries = matrix(c(1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))

  tv.seq= seq(from = 0, to = 1, by = 1)
  rho<-sapply(c(x1,x2), function(x) bayesRatio(x))

  dWit(t = rep(0:1, each = n), rho =rho, seed=NULL ,tv.seq=tv.seq , estimator.type = "asymptotic-tv-search")$tvhat
}, grid[,1], grid[,2])

set.seed(123, "L'Ecuyer")
res_binomial_1 <- mapply(FUN = function(n, gamma) {

  print(paste(toString(n),toString(gamma) ))

  p1 <- (1/2)*n^{-gamma}
  weights <- c((1/2)-p1, (1/2)+p1)

  x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
  x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
  d.right <- generateUnifMixturesDensities(boundaries = matrix(c(1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

  bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))

  rho<-sapply(c(x1,x2), function(x) bayesRatio(x))

  dWit(t = rep(0:1, each = n), rho = rho, seed=NULL, estimator.type = "binomial-test")$tvhat
}, grid[,1], grid[,2])


# gathering data
power.data <- data.table(logn = log(grid[,1], base = 10),
                         loglambda = -log(grid[,1], base = 10)*grid[,2],
                         tv_search = res_tv_search_1, tv_binomial = res_binomial_1)
power.table <- power.data[,.(power_search = mean(tv_search>0), power_binomial = mean(tv_binomial>0)),by=c("logn","loglambda")]

save.image(file="DATA_SIMULATION_1_SC1.Rdata")


}



if (RUN_SC_2) {
  ##### scenario 2: asymptotic-tv-search should win

  # running simulations
  set.seed(123, "L'Ecuyer")
  res_tv_search_2 <- mapply(FUN = function(n, gamma) {

    print(paste(toString(n),toString(gamma) ))

    p1 <- (1/2)*n^{-gamma}
    p2 <- (1/2) + (1/2)*n^{-gamma}
    weights <- c(p1, (1-p1)*p2, (1-p1)*(1-p2))

    x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-4, -3, -2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
    x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(3, 4, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

    d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-4, -3, -2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
    d.right <- generateUnifMixturesDensities(boundaries = matrix(c(3, 4, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

    bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))

    tv.seq=seq(from = 0, to = 1, by = 1)
    rho=sapply(c(x1,x2), function(x) bayesRatio(x))

    dWit(t = rep(0:1, each = n), rho = rho, seed=NULL, tv.seq=tv.seq, estimator.type = "asymptotic-tv-search")$tvhat
  }, grid[,1], grid[,2])

  set.seed(123, "L'Ecuyer")
  res_binomial_2 <- mapply(FUN = function(n, gamma) {

    print(paste(toString(n),toString(gamma) ))

    p1 <- (1/2)*n^{-gamma}
    p2 <- (1/2) + (1/2)*n^{-gamma}
    weights <- c(p1, (1-p1)*p2, (1-p1)*(1-p2))

    x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-4, -3, -2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
    x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(3, 4, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)


    d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-4, -3, -2, -1, 1, 2), ncol=2, byrow = TRUE), weights = weights)
    d.right <- generateUnifMixturesDensities(boundaries = matrix(c(3, 4, 1, 2, -2, -1), ncol=2, byrow = TRUE), weights = weights)

    bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))


    rho<-sapply(c(x1,x2), function(x) bayesRatio(x))

    dWit(t = rep(0:1, each = n), rho =rho, seed=NULL , estimator.type = "binomial-test")$tvhat
  }, grid[,1], grid[,2])


  # gathering data
  power.data <- data.table(logn = log(grid[,1], base = 10),
                           loglambda = -log(grid[,1], base = 10)*grid[,2],
                           tv_search = res_tv_search_2, tv_binomial = res_binomial_2)
  power.table <- power.data[,.(power_search = mean(tv_search>0), power_binomial = mean(tv_binomial>0)),by=c("logn","loglambda")]


  save.image(file="DATA_SIMULATION_1_SC2.Rdata")

}




if (RUN_SC_3) {
  # scenario 3: comtamination case

  # running simulations
  set.seed(123, "L'Ecuyer")
  res_tv_search_3 <- mapply(FUN = function(n, gamma) {

    print(paste(toString(n),toString(gamma) ))

    p1 <- (1/2)*n^{-gamma}
    weights <- c(p1, 1-p1)

    x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-2, -1, 1 , 2), ncol=2, byrow = TRUE), weights = weights)
    x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(1, 2), ncol=2, byrow = TRUE), weights = c(1))

    d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-2, -1, 1 , 2), ncol=2, byrow = TRUE), weights = weights)
    d.right <- generateUnifMixturesDensities(boundaries = matrix(c(1, 2), ncol=2, byrow = TRUE), weights = c(1))

    bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))

    tv.seq= seq(from = 0, to = 1, by = 1)
    rho=sapply(c(x1,x2), function(x) bayesRatio(x))

    dWit(t = rep(0:1, each = n), rho =rho, seed=NULL , tv.seq=tv.seq, estimator.type = "asymptotic-tv-search")$tvhat
  }, grid[,1], grid[,2])

  set.seed(123, "L'Ecuyer")
  res_binomial_3 <- mapply(FUN = function(n, gamma) {

    print(paste(toString(n),toString(gamma) ))

    p1 <- (1/2)*n^{-gamma}
    weights <- c(p1, 1-p1)

    x1 <- generateUnifMixturesData(n=n, boundaries = matrix(c(-2, -1, 1 , 2), ncol=2, byrow = TRUE), weights = weights)
    x2 <- generateUnifMixturesData(n=n, boundaries = matrix(c(1, 2), ncol=2, byrow = TRUE), weights = c(1))

    d.left <- generateUnifMixturesDensities(boundaries = matrix(c(-2, -1, 1 , 2), ncol=2, byrow = TRUE), weights = weights)
    d.right <- generateUnifMixturesDensities(boundaries = matrix(c(1, 2), ncol=2, byrow = TRUE), weights = c(1))

    bayesRatio <- function(x) d.right(x) / (d.left(x) + d.right(x))

    rho=sapply(c(x1,x2), function(x) bayesRatio(x))

    dWit(t = rep(0:1, each = n), rho = rho, seed=NULL, estimator.type = "binomial-test")$tvhat
  }, grid[,1], grid[,2])


  # gathering data
  power.data <- data.table(logn = log(grid[,1], base = 10),
                           loglambda = -log(grid[,1], base = 10)*grid[,2],
                           tv_search = res_tv_search_3, tv_binomial = res_binomial_3)
  power.table <- power.data[,.(power_search = mean(tv_search>0), power_binomial = mean(tv_binomial>0)),by=c("logn","loglambda")]


  save.image(file="DATA_SIMULATION_1_SC3.Rdata")

}



