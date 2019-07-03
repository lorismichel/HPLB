#' Get Distance
#' @param labels the labels of the classes, should be encoded as 0-->nclass-1
#' @param ordering.array array of size (nclass, nclass, nobs) such that the value (i,j,k) represents a propensity of beeing of class j instead of i for observation k.
#' @param alpha type I error control
#' @param computation.type
#' @export
getTvLbDistanceMatrix <- function(labels,
                                  ordering.array,
                                  alpha = 0.05,
                                  computation.type = "non-optimized") {

  if (dim(ordering.array)[1] != dim(ordering.array)[2]) {
    stop("Invalid ordering array.")
  }
  # rest of tests



  # returned value
  dist.mat <- matrix(0, nrow = dim(ordering.array)[1], ncol = dim(ordering.array)[2])
  if (computation.type == "non-optimized") {

    nclass <- dim(ordering.array)[1]


    for (i in 0:(nclass-2)) {
      for (j in (i+1):(nclass-1)) {

        ind.pairs <- which(labels %in% c(i,j))


        dist.mat[i+1, j+1] <- dWit(t = as.numeric(labels == j)[ind.pairs],
                                    rho = ordering.array[i+1,j+1,ind.pairs],
                                    s = 0.5, alpha = alpha,
                                    estimator.type = "tv-search")$tvhat

      }
    }

  }


  return((dist.mat + t(dist.mat)))
}
