#' Get Distance
#' @param labels the labels of the classes, should be encoded as 0-->nclass-1
#' @param ordering.array array of size (nclass, nclass, nobs) such that the value (i,j,k) represents a propensity of beeing of class j instead of i for observation k.
#' @param alpha type I error control
#' @param computation.type
#' @param ...
#' @export
getTvLbDistanceMatrix <- function(labels,
                                  ordering.array,
                                  alpha = 0.05,
                                  computation.type = "non-optimized",
                                  ...) {

  if (is.null(dim(ordering.array))) {
    dist.mat <- matrix(0, nrow = length(unique(labels)), ncol = length(unique(labels)))

    nclass <- length(unique(labels))


    for (i in 0:(nclass-2)) {
      for (j in (i+1):(nclass-1)) {

        ind.pairs <- which(labels %in% c(i,j))

        dist.mat[i+1, j+1] <- dWit(t = as.numeric(labels == j)[ind.pairs],
                                   rho = ordering.array[ind.pairs],
                                   s = 0.5, threshold = mean(ordering.array[ind.pairs]),...)$tvhat
       # dist.mat[i+1, j+1] <- dWit(t = as.numeric(labels == j)[ind.pairs],
        #                           rho = ordering.array[ind.pairs],
        #                           s = mean(ordering.array[ind.pairs]), ...)$tvhat

      }
    }

    return((dist.mat + t(dist.mat)))

  }
  # rest of tests
  if (dim(ordering.array)[1] != dim(ordering.array)[2]) {
    stop("Invalid ordering array.")
  }


  # returned value
  dist.mat <- matrix(0, nrow = dim(ordering.array)[1], ncol = dim(ordering.array)[2])
  if (computation.type == "non-optimized") {

    nclass <- dim(ordering.array)[1]


    for (i in 0:(nclass-2)) {
      for (j in (i+1):(nclass-1)) {

        if (sum(labels==i)<=1 | sum(labels==j)<=1) {
          warning("at least one class has less than one representative, lower-bound set to 0 by default.")
          dist.mat[i+1, j+1] <- 0
        } else {
          ind.pairs <- which(labels %in% c(i,j))

          dist.mat[i+1, j+1] <- dWit(t = as.numeric(labels == j)[ind.pairs],
                                     rho = ordering.array[i+1,j+1,ind.pairs],
                                     s = 0.5)$tvhat
        }
      }
    }

  }


  return((dist.mat + t(dist.mat)))
}
