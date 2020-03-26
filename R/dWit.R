#' Lower Bounds for Total Variation Distance (and Distributional Witnesses) Based on Samples
#'
#' @param t a numeric vector value corresponding to an ordering of the observations. For a two-sample
#'   test 0-1 numeric values values should be provided.
#' @param rho a numeric vetor value giving an ordering. This could be
#'   a classifier, a regressor, a witness function from a MMD kernel or anything else that would witness a distributional difference.
#' @param s a numeric vector value giving split points on t.
#' @param estimator.type a character value indicating which estimator to use.
#'        For total variation lower-bounds can be either "binomial-test",
#'        "asymptotic-tv-search", "empirical-tv-search", "custom-tv-search" or "hypergeometic-test".
#'        For distributional witness estimation "asymptotic-dwit-search".
#' @param alpha a numeric value giving the overall type-I error control level.
#' @param tv.seq a sequence of values between 0 and 1 used as the grid search for the total variation distance in case of tv-search.
#' @param custom.bounding.seq a list of bounding functions respecting the order of tv.seq used in case of estimator.type "custom-tv-search".
#' @param direction a character vector value made of "left" or "right" giving which distribution witness count to estimate (t<=s or t>s?).
#' @param threshold a numeric value. This is the threshold used if misclassification error in used.
#' @param z an integer value. This is the z at which the "hypergeometric-test" is applied.
#' @param verbose.plot a boolean value for additional plots.
#' @param seed an integer value. The seed for reproducility.
#'
#' @return a list containing the relevant lower-bounds estimates.
#'
#' @import data.table
#'
#' @examples
#'
#' ## libs
#' library(dWit)
#' library(ranger)
#' library(distrEx)
#'
#' ## reproducibility
#' set.seed(0)
#'
#' # TV lower bound based on two samples (binomial-test), Gaussian mean-shift example
#'
#' n <- 1000
#' means <- rep(c(0,2), each = n / 2)
#' x <- stats::rnorm(n, mean = means)
#' t <- rep(c(0,1), each = n / 2)
#'
#' bayesRate <- function(x) return(stats::dnorm(x, mean = 2) / (stats::dnorm(x, mean = 2) + stats::dnorm(x, mean = 0)))
#'
#' tvhat <- dWit(t = t, rho = bayesRate(x), estimator.type = "binomial-test")
#'
#' # optimal mixture detection (asymptotic-tv-search), Gaussian mean-shift example
#'
#' n <- 1000
#' mean.shift <- 2
#' t.train <- runif(n, 0 ,1)
#' x.train <- ifelse(t.train>0.5, stats::rnorm(n, mean.shift), stats::rnorm(n))
#' rf <- ranger::ranger(t~x, data.frame(t=t.train,x=x.train))
#'
#' n <- 1000
#' t.test <- runif(n, 0 ,1)
#'x.test <- ifelse(t.test>0.5, stats::rnorm(n, mean.shift), stats::rnorm(n))
#' rho <- predict(rf, data.frame(t=t.test,x=x.test))$predictions

#' ## out-of-sample
#' tv.oos <- dWit(t = t.test, rho = rho, s = seq(0.1,0.9,0.1), estimator.type = "asymptotic-tv-search")
#'
#' # oob
#' tv.oob <- dWit(t = t, rho = rf$predictions, s = seq(0.1,0.9,0.1), estimator.type = "asymptotic-tv-search")
#'
#'
#'
#' ## total variation values
#' tv <- c()
#' for (s in seq(0.1,0.9,0.1)) {
#'
#'  if (s<=0.5) {
#'
#'    D.left <- Norm(0,1)
#'  } else {
#'
#'    D.left <- UnivarMixingDistribution(Dlist = list(Norm(0,1),Norm(mean.shift,1)),
#'                mixCoeff = c(ifelse(s<=0.5, 1, 0.5/s), ifelse(s<=0.5, 0, (s-0.5)/s)))
#'  }
#'  if (s<0.5) {
#'
#'    D.right <- UnivarMixingDistribution(Dlist = list(Norm(0,1),Norm(mean.shift,1)),
#'                mixCoeff = c(ifelse(s<=0.5, (0.5-s)/(1-s), 0), ifelse(s<=0.5, (0.5/(1-s)), 1)))
#'  } else {
#'
#'    D.right <- Norm(mean.shift,1)
#'  }
#' tv <- c(tv, TotalVarDist(e1 = D.left, e2 = D.right))
#' }
#'
#' ## plot
#' par(mfrow=c(2,1))
#' plot(t.test,x.test,pch=19,xlab="t",ylab="x")
#' plot(seq(0.1,0.9,0.1), tv.oos$tvhat,type="l",ylim=c(0,1),xlab="t", ylab="TV")
#' lines(seq(0.1,0.9,0.1), tv.oob$tvhat,type="l",col="blue")
#' lines(seq(0.1,0.9,0.1), tv, col="red",type="l")

#'
#' @export
dWit <- function(t,
                 rho,
                 s                   = 0.5,
                 estimator.type      = "asymptotic-tv-search",
                 alpha               = 0.05,
                 tv.seq              = seq(from = 0, to = 1, by = 1/length(t)),
                 custom.bounding.seq = NULL,
                 direction           = rep("left", length(s)),
                 threshold           = 0.5,
                 verbose.plot        = FALSE,
                 z                   = length(t)/2,
                 seed                = 0,
                 ...) {

  ## parameters checks

  # seed check
  if (!is.numeric(seed) || any(!is.finite(seed))  || length(seed) != 1) {
    stop("Invalid seed, please see ?dWit.")
  }

  # setting seed for reproducility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # estimator.type check
  if (!is.character(estimator.type) || !(estimator.type %in% c("asymptotic-tv-search",
                                                               "custom-tv-search",
                                                               "empirical-tv-search",
                                                               "binomial-test",
                                                               "asymptotic-binomial-test",
                                                               "hypergeometric-test",
                                                               "asymptotic-dwit-search"
  ))) {
    stop("Invalid estimator type, please see ?dWit.")
  }

  # t check
  if (!is.numeric(t) || any(!is.finite(t)) || is.matrix(t)) {
    stop("Invalid values or type for t, please see ?dWit.")
  }

  # s check
  if (!is.numeric(s) || any(!is.finite(s)) || any(s > max(t)) || any(s < min(t)) || is.matrix(s)) {
    stop("Invalid values or type for s, please see ?dWit.")
  }

  # z check
  if ((estimator.type == "hypergeometric-test") && (!is.numeric(z) || any(!is.finite(z)) || (length(s) != length(z)) || any(z < 0) || any(z > length(t)))) {
    stop("Invalid values or type for z, please check ?dWit.")
  }

  # rho check
  if (!is.numeric(rho) || any(!is.finite(rho))) {
    stop("Invalid values for rho, please see ?dWit.")
  }

  # two cases depending on the dimension of rho
  if (is.matrix(rho)) {

    if (nrow(rho) != length(t)) {
      stop("Incompatible dimensions between t and rho, please see ?dWit.")
    }
    if (ncol(rho) != length(s)) {
      stop("Incompatible dimensions between s and rho, please see ?dWit.")
    }
    ordering.type = "multiple"

  } else {
    if (length(rho) != length(t)) {
      stop("Incompatible dimensions between t and rho, please see ?dWit.")
    }
    ordering.type = "single"
  }

  # ordering.type check
  if (estimator.type == "binomial-test" && ordering.type == "multiple") {
    stop("Invalid ordering type for binomial-test estimator, please see ?dWit.")
  }

  # alpha check
  if (!is.numeric(alpha) || (length(alpha) != 1) || (alpha > 1) || (alpha <= 0)) {
    stop("Invalid type I error level alpha, please see ?dWit.")
  }

  # tv.seq check
  if (!is.numeric(tv.seq) || any(!is.finite(tv.seq)) || any(tv.seq > 1) || any(tv.seq < 0) || any(duplicated(tv.seq))) {
    stop("Invalid grid for tv, please see ?dWit.")
  }

  # direction check
  if ((estimator.type == "asymptotic-dwit-search") && (!is.character(direction) || !(direction %in% c("left", "right")) || (length(s) != length(direction)))) {
    stop("Invalid direction, please see ?dWit.")
  }

  # check that with user-defined bounding sequence or empirical-tv-seach s should be uni-dim
  if (estimator.type %in% c("empirical-tv-search", "custom-tv-search")) {
    if (length(s) != 1) {
      stop("Incompatible value of s and estimator.type, please see ?dWit.")
    }
    if (estimator.type == "custom-tv-search") {
      if (is.null(custom.bounding.seq)) {
        stop("Invalid custom.bounding.seq, please see ?dWit.")
      } else {
        if (length(custom.bounding.seq) != length(tv.seq)) {
          stop("Incompatible dimensions between custom.bounding.seq and tv.seq, please see ?dWit.")
        }
      }
    }

    if (estimator.type %in% c("empirical-tv-search")) {
      custom.bounding.seq <- empiricalBF(tv.seq = tv.seq,
                                         alpha = alpha,
                                         m = sum(t <= s),
                                         n = sum(t > s), ...)
    }

    # dimensions of s, k and tresholds check
    if (length(s) != length(k)) {
      stop("Incompatible dimensions between s, and k, please see ?dWit.")
    }
  }

  ## core of the code

  # sort in increasing order the grid of tv values
  tv.seq <- sort(tv.seq, decreasing = FALSE)

  # define the returned values
  dwithat.left  <- rep(NA, length(s))
  dwithat.right <- rep(NA, length(s))
  tvhat.left    <- rep(NA, length(s))
  tvhat.right   <- rep(NA, length(s))
  tvhat         <- rep(NA, length(s))

  # main loop over the different splits
  for (k in 1:length(s)) {

    # first we obtain the ranks from the rho
    if (ordering.type == "single") {

      ranks <- rank(rho, ties.method = "random")
    } else {

      ranks <- rank(rho[,k], ties.method = "random")
    }

    # then we cast everything in a data.table object
    ranks.table <- data.table::data.table(rank = ranks, t = t)

    # let us order by rank the table
    data.table::setkey(ranks.table, "rank")

    # getting the parameters m and n
    m <- nrow(ranks.table[t <= s[k]])
    n <- nrow(ranks.table) - m

    if (m == 0 || n == 0) {
      warning(paste0("m or n is zero for split ", s, ", please see ?dWit."))
      next
    }


    ## binomial estimator
    if (estimator.type == "binomial-test") {

      # compute the decisions
        decision <- ifelse(rho == 0.5,
                           sample(c(TRUE, FALSE), size = length(rho), replace = TRUE),
                           rho > threshold)


      # compute the rightly classified number of points per class
      obs0 <- sum(!decision & t == 0)
      obs1 <- sum(decision & t == 1)

      # probability for class 0
      phat0 <- stats::binom.test(x = obs0,
                                 n = sum(t == 0),
                                 p = 0.5,
                                 alternative = "greater",
                                 conf.level = 1 - (alpha / 2))$conf.int[1]

      # probability for class 1
      phat1 <- stats::binom.test(x = obs1,
                                 n = sum(t == 1),
                                 p = 0.5,
                                 alternative = "greater",
                                 conf.level = 1 - (alpha / 2))$conf.int[1]

      # compute the tvhat value
      tvhat <- max(phat0 + phat1 - 1, 0)

    }



    ## binomial estimator
    if (estimator.type == "asymptotic-binomial-test") {

      # compute the decisions
      decision <- ifelse(rho == 0.5,
                         sample(c(TRUE, FALSE), size = length(rho), replace = TRUE),
                         rho > threshold)


      # probability for class 0
      Ahat0<-sum(!decision & t == 0)/sum(t==0)

      # probability for class 1
      Ahat1<-sum(decision & t == 1)/sum(t==1)

      # compute the Variance
      sigma= sqrt(Ahat0 * (1-Ahat0)/sum(t==0)  +  Ahat1 * (1-Ahat1)/sum(t==1) )



      # compute the tvhat value
      tvhat <- max(Ahat0 + Ahat1 - 1 -  qnorm(1-alpha)*sigma, 0)

    }



    ## branching on estimators using the V-function

    # define the V function corresponding to m and z (and n and z respectively)
    V.left  <- cumsum(ranks.table$t <= s[k])
    V.right <- cumsum(rev(ranks.table$t > s[k]))

    if (estimator.type == "asymptotic-tv-search") {

      # recursive search
      max.stat <- 1
      step     <- 1

      while (max.stat > 0) {

        # search updates
        tv.cur <- tv.seq[step]
        step   <- step + 1

        # upper-bound dWit left and right
        dWit.left  <- stats::qbinom(p = 1 - (alpha / 3), size = m, prob = tv.cur)
        dWit.right <- stats::qbinom(p = 1 - (alpha / 3), size = n, prob = tv.cur)

        # left over hypergeometric parameters m0 and n0 in the middle
        m0 <- m - dWit.left
        n0 <- n - dWit.right

        if ((m0 <= 2) | (n0 <= 2)) {
          warning(paste0("m0 or n0 smaller or equal to 2 for split: ", s[k], ", asymtotic regime may not hold, see ?dWit."))
          break
        }

        # asymptotic parameters
        x.alpha <- -log(-log(1-alpha/3)/2)
        beta.left <- sqrt(2 * log(log(m0))) + (log(log(log(m0)))-log(pi)+2*x.alpha) / (2 * sqrt(2 * log(log(m0))))
        beta.right <- sqrt(2 * log(log(n0))) + (log(log(log(n0)))-log(pi)+2*x.alpha) / (2 * sqrt(2 * log(log(n0))))

        # mean and sd under the null (asymptotics)
        sd0 <- sapply(1:(m0+n0), function(z) sqrt(z * (m0 /(m0+n0)) * (n0 /(m0+n0)) * ((m0+n0-z)/(m0+n0-1))))
        mean0 <- sapply(1:(m0+n0), function(z) z * (m0 /(m0+n0)))

        # creating the bounding function
        Q <- rep(0, m + n)

        if (dWit.left != 0) {
          Q[1:dWit.left] <- 1
          Q <- cumsum(Q)
        }

        # hypergeometric in the middle
        ind.hypergeom <- c((dWit.left+1):(m+n-dWit.right))
        Q[ind.hypergeom] <- Q[ind.hypergeom] + mean0 + beta.left * sd0

        # witnesses right
        Q[(m+n-dWit.right):(m+n)] <- m
        Q <- pmin(Q, m)

        max.stat <- max(V.left-Q)

        if (verbose.plot && (tv.cur == 0)) {
          plot(V.left - mean0, type="l",
               main="Observed centered V-function with asymptotic bounding function", xlab="z", ylab="V(z)")
          lines(beta.left * sd0, col="blue")
        }
      }

      # return the tv estimate
      tvhat[k] <- tv.cur

    } else if (estimator.type %in% c("custom-tv-search", "empirical-tv-search")) {

      # recursive search
      max.stat <- 1
      step     <- 1

      while (max.stat > 0 & step != length(tv.seq)) {

        tv.cur <- tv.seq[step]
        max.stat <- max(V.left - custom.bounding.seq[[step]])
        step <- step + 1

        if (verbose.plot && (tv.cur == 0)) {
          plot(V.left - mean0, type="l",
               main="Observed centered V-function with asymptotic bounding function",xlab="z",ylab="V(z)")
          lines(beta.left * sd0, col="blue")
        }
      }

      # return the tv estimate
      tvhat[k] <- tv.cur

    } else if (estimator.type == "hypergeometric-test") {

      # recursive search
      max.stat <- 1
      step     <- 1

      while (max.stat > 0) {

        # search updates
        tv.cur <- tv.seq[step]
        step   <- step + 1

        # constructing witness overestimations
        lambda.left  <- stats::qbinom(p = 1 - (alpha / 3), size = m, prob = tv.cur)
        lambda.right <- stats::qbinom(p = 1 - (alpha / 3), size = n, prob = tv.cur)

        # left over hypergeometric parameters m0 and n0 in the center
        m0 <- m - lambda.left
        n0 <- n - lambda.right

        # constructing hypergeometric overestimation
        q.left <- stats::qhyper(p = 1 - (alpha / 3), m = m0, n = n0, k = z[k] - dwit.left)

        # overestimation
        Q <- dwit.left + q.left

        # new observed statistic
        max.stat <- V.left[z[k]] - Q
      }

      # return the tv estimate
      tvhat[k] <- tv.cur

    } else if (estimator.type == "asymptotic-dwit-search") {

      if (direction[k] == "left") {

        m <- nrow(ranks.table[t <= s[k]])
        n <- nrow(ranks.table) - m

      } else if (direction[k] == "right") {

        n <- nrow(ranks.table[t <= s[k]])
        m <- nrow(ranks.table) - n
      }

      # recursive search default values
      # TODO: implement a faster search (dyadic)
      max.stat <- 1
      step     <- 1
      dwit.left.grid <- c(0:m)

      while (max.stat > 0) {

        # search updates
        dwit.left <- dwit.left.grid[step]
        step      <- step + 1

        # overestimating the tv for the number of right dwit
        tv.cur <- stats::binom.test(x = dwit.left,
                                    n = m,
                                    p = 0.5,
                                    alternative = "less",
                                    conf.level = 1 - (alpha / 3))$conf.int[2]

        # dwit right
        dwit.right <- stats::qbinom(p    = 1 - (alpha / 3),
                                    size = n,
                                    prob = tv.cur)

        # define mo and no for the hypergeometric part
        m0 <- m - dwit.left
        n0 <- n - dwit.right

        # asymptotic values
        x.alpha <- -log(-log(1-alpha/3)/2)
        beta.left <- sqrt(2 * log(log(m0))) + (log(log(log(m0)))-log(pi)+2*x.alpha) / (2 * sqrt(2 * log(log(m0))))
        beta.right <- sqrt(2 * log(log(n0))) + (log(log(log(n0)))-log(pi)+2*x.alpha) / (2 * sqrt(2 * log(log(n0))))

        # mean and sd under the null (asymptotics)
        sd0 <- sapply(1:(m0+n0), function(z) sqrt(z * (m0 /(m0+n0)) * (n0 /(m0+n0)) * ((m0+n0-z)/(m0+n0-1))))
        mean0 <- sapply(1:(m0+n0), function(z) z * (m0 /(m0+n0)))

        # creating the bounding function
        Q <- rep(0, m + n)

        # adding witnesses left
        if (dwit.left != 0) {
          Q[1:dwit.left] <- 1
          Q <- cumsum(Q)
        }

        # hypergeometric filling in the middle
        ind.hypergeom <- c((dwit.left+1):(m+n-dwit.right))
        Q[ind.hypergeom] <- Q[ind.hypergeom] + mean0 + beta.left * sd0

        # adding witnesses right
        Q[(m+n-dwit.right):(m+n)] <- m
        Q <- pmin(Q, m)

        if (direction[k] == "left") {
          max.stat <- max(V.left - Q)
        } else if (direction[k] == "right") {
          max.stat <- max(V.right - Q)
        }

      }

      # updating the value for the split depending on the direction
      if (direction[k] == "left") {
        dwithat.left[k] <- dwit.left
      } else if (direction[k] == "right") {
        dwithat.right[k] <- dwit.left
      }

    }
  }

  # return a list with relevant estimates
  return(Filter(x = list(dwithat.left  = dwithat.left,
                         dwithat.right = dwithat.right,
                         tvhat.left    = tvhat.left,
                         tvhat.right   = tvhat.right,
                         tvhat         = tvhat),
                f = function(ll) !any(is.na(ll))))
}
