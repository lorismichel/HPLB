#' Pairwise Total Variation Distance Lower Bound Matrix for the Multi-Class Setting
#'
#' @param labels a numeric vector value. The labels of the classes, should be encoded in [0,nclass-1].
#' @param ordering.array a numeric array of size (nclass, nclass, nobs) such that the value (i,j,k) represents a propensity of being of class j instead of i for observation k.
#' @param alpha a numeric value. The type-I error level.
#' @param computation.type a character value. For the moment only "non-optimized" (default) available.
#' @param seed an integer value. The seed for reproducility.
#' @param ... additional parameters to be passed to the HPLB function.
#'
#' @return a numeric matrix of size (nclass, nclass) giving the matrix of pairwise total variation lower bounds.
#'
#' @author Loris Michel, Jeffrey Naef
#'
#' @references
#' L. Michel, J. Naef and N. Meinshausen (2020). High-Probability Lower Bounds for the Total Variation Distance \cr
#'
#' @examples
#'  # iris example
#'  require(HPLB)
#'  require(ranger)
#'
#'  # training a multi-class classifier on iris and getting tv lower bounds between classes
#'  data("iris")
#'
#'  ind.train <- sample(1:nrow(iris), size = nrow(iris)/2, replace = FALSE)
#'
#'  rf <- ranger(Species~., data = iris[ind.train, ], probability = TRUE)
#'  preds <- predict(rf, iris[-ind.train,])$predictions
#'
#'  # creating the ordering array based on prediction differences
#'  ar <- array(dim = c(3, 3, nrow(preds)))
#'  for (i in 1:3) {
#'    for (j in 1:3) {
#'     ar[i,j,] <- preds[,j] - preds[,i]
#'    }
#'  }
#'
#'  # encoding the class response
#'  y <- factor(iris$Species)
#'  levels(y) <- c(0,1,2)
#'  y <- as.numeric(y)-1
#'
#'  # getting the lower bound matrix
#'  tvhat.iris <- HPLBmatrix(labels = y[-ind.train], ordering.array = ar)
#'  tvhat.iris
#' @export
HPLBmatrix <- function(labels,
                       ordering.array,
                       alpha            = 0.05,
                       computation.type = "non-optimized",
                       seed             = 0,
                       ...) {

  ## parameters checks

  # seed check
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("Invalid seed, please see ?HPLB.")
  }

  # setting seed for reproducility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # labels check
  if (!is.numeric(labels) || any(!is.finite(labels)) || is.matrix(labels) || any(!(labels %in% c(0:length(unique(labels)-1))))) {
    stop("Invalid values or type for labels, please see ?HPLBmatrix.")
  }

  # ordering.array check
  if (!is.numeric(ordering.array) || any(!is.finite(ordering.array)) || !is.array(ordering.array) || (length(dim(ordering.array)) != 3)
      || (dim(ordering.array)[1] != dim(ordering.array)[2]) || (dim(ordering.array)[1] != length(unique(labels))) || (length(labels) != dim(ordering.array)[3])) {
    stop("Invalid values or type for ordering.array or compatibilities with labels, please see ?HPLBmatrix.")
  }

  # computation.type check
  if (!is.character(computation.type) || (length(computation.type) != 1) || !(computation.type %in% c("non-optimized"))) {
    stop("Invalid values for computation.type, please see ?HPLBmatrix.")
  }

  # alpha check
  if (!is.numeric(alpha) || (length(alpha) != 1) || (alpha > 1) || (alpha <= 0)) {
    stop("Invalid type I error level alpha, please see ?HPLBmatrix.")
  }



  # get the number of classes
  nclass <- length(unique(labels))

  # create the tv matrix
  if (is.null(dim(ordering.array))) {

    tvhat.mat <- matrix(0, nrow = nclass, ncol = nclass)

    # loop over the different class pairs
    for (i in 0:(nclass-2)) {

      for (j in (i+1):(nclass-1)) {

        # selecting the observations with i or j
        ind.pairs <- which(labels %in% c(i, j))

        # lower bound
        tvhat.mat[i+1, j+1] <- HPLB(t         = as.numeric(labels == j)[ind.pairs],
                                    rho       = ordering.array[ind.pairs],
                                    s         = 0.5,
                                    ...)$tvhat
      }
    }

    # return a symmetrized matrix
    return((tvhat.mat + t(tvhat.mat)))

  }

  # rest of tests
  if (dim(ordering.array)[1] != dim(ordering.array)[2]) {
    stop("Invalid ordering array.")
  }

  # returned value
  tvhat.mat <- matrix(0,
                     nrow = dim(ordering.array)[1],
                     ncol = dim(ordering.array)[2])


  if (computation.type == "non-optimized") {

    for (i in 0:(nclass-2)) {

      for (j in (i+1):(nclass-1)) {

        if (sum(labels == i) <= 1 | sum(labels == j) <= 1) {

          warning("at least one class has less than one representative, lower bound set to 0 by default, see ?HPLBmatrix.")

          tvhat.mat[i+1, j+1] <- 0
        } else {

          ind.pairs <- which(labels %in% c(i,j))

          tvhat.mat[i+1, j+1] <- HPLB(t   = as.numeric(labels == j)[ind.pairs],
                                      rho = ordering.array[i+1,j+1,ind.pairs],
                                      s   = 0.5,
                                      ...)$tvhat
        }
      }
    }

  }

  # return a symmetrized matrix
  return((tvhat.mat + t(tvhat.mat)))
}
