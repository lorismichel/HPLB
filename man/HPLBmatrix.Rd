% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HPLBmatrix.R
\name{HPLBmatrix}
\alias{HPLBmatrix}
\title{Pairwise Total Variation Distance Lower Bound Matrix for the Multi-Class Setting}
\usage{
HPLBmatrix(
  labels,
  ordering.array,
  alpha = 0.05,
  computation.type = "non-optimized",
  seed = 0,
  ...
)
}
\arguments{
\item{labels}{a numeric vector value. The labels of the classes, should be encoded in [0,nclass-1].}

\item{ordering.array}{a numeric array of size (nclass, nclass, nobs) such that the value (i,j,k) represents a propensity of being of class j instead of i for observation k.}

\item{alpha}{a numeric value. The type-I error level.}

\item{computation.type}{a character value. For the moment only "non-optimized" (default) available.}

\item{seed}{an integer value. The seed for reproducility.}

\item{...}{additional parameters to be passed to the HPLB function.}
}
\value{
a numeric matrix of size (nclass, nclass) giving the matrix of pairwise total variation lower bounds.
}
\description{
Pairwise Total Variation Distance Lower Bound Matrix for the Multi-Class Setting
}
\examples{
 # iris example
 require(HPLB)
 require(ranger)

 # training a multi-class classifier on iris and getting tv lower bounds between classes
 data("iris")

 ind.train <- sample(1:nrow(iris), size = nrow(iris)/2, replace = FALSE)

 rf <- ranger(Species~., data = iris[ind.train, ], probability = TRUE)
 preds <- predict(rf, iris[-ind.train,])$predictions

 # creating the ordering array based on prediction differences
 ar <- array(dim = c(3, 3, nrow(preds)))
 for (i in 1:3) {
   for (j in 1:3) {
    ar[i,j,] <- preds[,j] - preds[,i]
   }
 }

 # encoding the class response
 y <- factor(iris$Species)
 levels(y) <- c(0,1,2)
 y <- as.numeric(y)-1

 # getting the lower bound matrix
 tvhat.iris <- HPLBmatrix(labels = y[-ind.train], ordering.array = ar)
 tvhat.iris
}
\references{
L. Michel, J. Naef and N. Meinshausen (2020). High-Probability Lower Bounds for the Total Variation Distance \cr
}
\author{
Loris Michel, Jeffrey Naef
}
