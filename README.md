# HPLB

[![Build Status](https://travis-ci.org/lorismichel/HPLB.svg?branch=master)](https://travis-ci.org/lorismichel/HPLB)
[![Build status](https://ci.appveyor.com/api/projects/status/jirtk3gmc4sdp0gl?svg=true)](https://ci.appveyor.com/project/lorismichel/hplb)

## Overview

HPLB is a package intended to provided high-probability lower bounds (HPLB) for the total variance distance (TV) based on finite samples. In particular, it implements the abc and bc estimators described in [Michel et al. 2020](https://arxiv.org/abs/2005.06006). The main idea is to compute HPLBs for TV from uni-dimensional projections that would practically be obtained from standard learning algorithms. For more information  the user can refer to the original paper. Examples of use of the library are shown below.


## Installation

The package should be (soon) available on CRAN, To install the package from github you can run
``` r
install.packages("devtools")
devtools::install_github("lorismichel/HPLB")
```

## Examples: 


We provide two examples, a shift in mean and a contamination example.

``` r
library(HPLB)
library(stats)
library(ranger)
library(distrEx)

## univariate shift in mean for normals using random forest as a lower-dimensional projceion
m <- n <- 500

x.train <- c(rnorm(n = m, mean = 0), rnorm(n = n, mean = 2))
y.train <- c(rep(0, m), rep(1, n))

x.test <- c(rnorm(n = m, mean = 0), rnorm(n = n, mean = 2))
y.test <- c(rep(0, m), rep(1, n))

# fitting a classification forest
rf <- ranger(factor(y)~.,data = data.frame(y = y.train, x = x.train))
rf.prob <- ranger(factor(y)~.,data = data.frame(y = y.train, x = x.train), probability = TRUE)

# getting the predictions on the test set (hard or soft)
preds.hard <- as.numeric(predict(rf, data.frame(x = x.test))$predictions)-1
preds.soft <- predict(rf.prob, data.frame(x = x.test))$predictions[,"1"]


# Total variation distance between N(0,1) and N(2,1)
TotalVarDist(Norm(0,1),Norm(2,1))

# getting lower-bounds on total variation distance N(0,1) and N(2,1) with different estimators

# binary classifier (bc)
HPLB(t = y.test, rho = preds.hard, estimator.type = "bc")

# adaptive binary classifier (bc)
HPLB(t = y.test, rho = preds.hard, estimator.type = "abc")


## contamination  using random forest as a lower-dimensional projceion
m <- n <- 500

x.train <- c(rnorm(n = m, mean = 0), ifelse(runif(n = n) <= 0.05, rnorm(n = n, mean = 5), rnorm(n = n, mean = 0)))
y.train <- c(rep(0, m), rep(1, n))

x.test <- c(rnorm(n = m, mean = 0), ifelse(runif(n = n) <= 0.05, rnorm(n = n, mean = 5), rnorm(n = n, mean = 0)))
y.test <- c(rep(0, m), rep(1, n))

# fitting a classification forest
rf <- ranger(factor(y)~.,data = data.frame(y = y.train, x = x.train))
rf.prob <- ranger(factor(y)~.,data = data.frame(y = y.train, x = x.train), probability = TRUE)

# getting the predictions on the test set (hard or soft)
preds.hard <- as.numeric(predict(rf, data.frame(x = x.test))$predictions)-1
preds.soft <- predict(rf.prob, data.frame(x = x.test))$predictions[,"1"]


# Total variation distance between N(0,1) and 0.05 x N(5,1) + 0.95 x N(0,1)
TotalVarDist(Norm(0,1), UnivarMixingDistribution(Norm(0,1),Norm(5,1), mixCoeff = c(0.95,0.05)))

# getting lower-bounds on total variation distance between N(0,1) and 0.05 x N(5,1) + 0.95 x N(0,1) with different estimators

# binary classifier (bc)
HPLB(t = y.test, rho = preds.hard, estimator.type = "bc")

# adaptive binary classifier (bc)
HPLB(t = y.test, rho = preds.hard, estimator.type = "abc")
```

## Issues

To report an issue, please use the [issue tracker](http://github.com/lorismichel/HPLB/issues) on github.com.
