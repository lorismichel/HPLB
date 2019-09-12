getDataSet <- function(dataset.name,
                       m.train = 100,
                       n.train = 100,
                       m.test = 100,
                       n.test = 100,
                       ...) {
  if (dataset.name == "UnivariateNormalShift") {
    t.train <- c(rep(0,m.train), rep(1,n.train))
    t.test <- c(rep(0,m.test), rep(1,n.test))
    x.train <- c(rnorm(m.train, mean = 0, sd = 1), rnorm(n.train, ...))
    x.test <- c(rnorm(m.test, mean = 0, sd = 1), rnorm(n.test, ...))
  } else if (dataset.name == "UnivariateNormalShift+Contamination") {
    params <- list(...)
    t.train <- c(rep(0,m.train), rep(1,n.train))
    t.test <- c(rep(0,m.test), rep(1,n.test))
    x.train <- c(rnorm(m.train, mean = 0, sd = 1), ifelse(runif(n.train)<=params$lambda,
                                                          rnorm(n.train, sd=params$sd, mean=params$mean + 4 * params$sd),
                                                          rnorm(n.train, sd=params$sd, mean=params$mean)))
    x.test <- c(rnorm(m.test, mean = 0, sd = 1), ifelse(runif(n.test)<=params$lambda,
                                                         rnorm(n.test, sd=params$sd, mean=params$mean + 4 * params$sd),
                                                         rnorm(n.test, sd=params$sd, mean=params$mean)))
  }

  return(list(train=data.frame(t=t.train, x = x.train),
              test=data.frame(t=t.test, x = x.test)))
}
