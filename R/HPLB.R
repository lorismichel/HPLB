#' High Probability Lower Bounds (HPLB) for the Total Variation Distance (TV) Based on Finite Samples
#'
#' Implementations of different HPLBs for TV as described in (Michel et al., 2020).
#'
#' @param t a numeric vector value corresponding to a natural ordering of the observations. For a two-sample
#'   test 0-1 numeric values values should be provided.
#' @param rho a numeric vector value providing an ordering. This could be
#'   a binary classifier, a regressor, a witness function from a MMD kernel or anything else that would witness a distributional difference.
#' @param s a numeric vector value giving split points on t.
#' @param estimator.type a character value indicating which estimator to use. One option out of:
#'   \itemize{
#'   \item{\code{adapt}:}{adaptive binary classification estimator (asymptotic bounding function)}
#'   \item{\code{bayes}:}{binary classification estimator}
#'   \item{\code{bayes_finite_sample}:}{binary classification finite sample estimator}
#'   \item{\code{adapt_empirical}:}{adaptive binary classification estimator (simulation-based bounding function)}
#'   \item{\code{adapt_custom}:}{adaptive binary classificatrion estimator (user-defined bounding function)}
#'   \item{\code{adapt_dwit}:}{adaptive binary classificatrion estimator (for distributional witnesses estimation)}
#'   }
#' @param alpha a numeric value giving the overall type-I error control level.
#' @param tv.seq a sequence of values between 0 and 1 used as the grid search for the total variation distance in case of tv-search.
#' @param custom.bounding.seq a list of bounding functions respecting the order of tv.seq used in case of estimator.type "custom-tv-search".
#' @param direction a character vector value made of "left" or "right" giving which distribution witness count to estimate (t<=s or t>s?).
#' @param cutoff a numeric value. This is the cutoff used if bayes estimators are used. The theory suggests to use 1/2 but this can be changed.
#' @param verbose.plot a boolean value for additional plots.
#' @param seed an integer value. The seed for reproducibility.
#' @param ... additional parameters for the function \code{empiricalBF}.
#'
#' @return a \code{list} containing the relevant lower bounds estimates. For the total variation distance the relevant entry is \code{tvhat}.
#'
#' @import data.table
#'
#' @author Loris Michel, Jeffrey Naef
#'
#' @references
#' L. Michel, J. Naef and N. Meinshausen (2020). High-Probability Lower Bounds for the Total Variation Distance \cr
#'
#' @examples
#' ## libs
#' library(HPLB)
#' library(ranger)
#' library(distrEx)
#'
#' ## reproducibility
#' set.seed(0)
#'
#' ## Example 1: TV lower bound based on two samples (bayes estimator), Gaussian mean-shift example
#'
#' n <- 100
#' means <- rep(c(0,2), each = n / 2)
#' x <- stats::rnorm(n, mean = means)
#' t <- rep(c(0,1), each = n / 2)
#'
#' bayesRate <- function(x) {
#'   return(stats::dnorm(x, mean = 2) /
#'     (stats::dnorm(x, mean = 2) + stats::dnorm(x, mean = 0)))
#' }
#'
#' # estimated HPLB
#' tvhat <- HPLB(t = t, rho = bayesRate(x), estimator.type = "bayes")
#' # true TV
#' TotalVarDist(e1 = Norm(2,1), e2 = Norm(0,1))
#'
#' ## Example 2: optimal mixture detection (adapt estimator), Gaussian mean-shift example
#'
#' n <- 100
#' mean.shift <- 2
#' t.train <- runif(n, 0 ,1)
#' x.train <- ifelse(t.train>0.5, stats::rnorm(n, mean.shift), stats::rnorm(n))
#' rf <- ranger::ranger(t~x, data.frame(t=t.train,x=x.train))
#'
#' n <- 100
#' t.test <- runif(n, 0 ,1)
#'x.test <- ifelse(t.test>0.5, stats::rnorm(n, mean.shift), stats::rnorm(n))
#' rho <- predict(rf, data.frame(t=t.test,x=x.test))$predictions
#'
#' ## out-of-sample
#' tv.oos <- HPLB(t = t.test, rho = rho, s = seq(0.1,0.9,0.1), estimator.type = "adapt")
#'
#' \donttest{
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
#'  if (s < 0.5) {
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
#' oldpar <- par(no.readonly =TRUE)
#' par(mfrow=c(2,1))
#' plot(t.test,x.test,pch=19,xlab="t",ylab="x")
#' plot(seq(0.1,0.9,0.1), tv.oos$tvhat,type="l",ylim=c(0,1),xlab="t", ylab="TV")
#' lines(seq(0.1,0.9,0.1), tv, col="red",type="l")
#' par(oldpar)
#' }
#' @export
HPLB <- function(t,
                 rho,
                 s                   = 0.5,
                 estimator.type      = "adapt",
                 alpha               = 0.05,
                 tv.seq              = seq(from = 0, to = 1, by = 1/length(t)),
                 custom.bounding.seq = NULL,
                 direction           = rep("left", length(s)),
                 cutoff              = 0.5,
                 verbose.plot        = FALSE,
                 seed                = 0,
                 ...) {

  ## parameters checks

  # seed check
  if (!is.null(seed) && (!is.numeric(seed) || any(!is.finite(seed))  || length(seed) != 1)) {
    stop("Invalid seed, please see ?HPLB.")
  }

  # setting seed for reproducility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # estimator.type check
  if (!is.character(estimator.type) || !(estimator.type %in% c("adapt",
                                                               "bayes",
                                                               "bayes_finite_sample",
                                                               "adapt_empirical",
                                                               "adapt_dwit"
  ))) {
    stop("Invalid estimator type, please see ?HPLB.")
  }

  # t check
  if (!is.numeric(t) || any(!is.finite(t)) || is.matrix(t)) {
    stop("Invalid values or type for t, please see ?HPLB.")
  }

  # s check
  if (!is.numeric(s) || any(!is.finite(s)) || any(s > max(t)) || any(s < min(t)) || is.matrix(s)) {
    stop("Invalid values or type for s, please see ?HPLB.")
  }

  # rho check
  if (!is.numeric(rho) || any(!is.finite(rho))) {
    stop("Invalid values for rho, please see ?HPLB.")
  }

  # two cases depending on the dimension of rho
  if (is.matrix(rho)) {

    if (nrow(rho) != length(t)) {
      stop("Incompatible dimensions between t and rho, please see ?HPLB.")
    }
    if (ncol(rho) != length(s)) {
      stop("Incompatible dimensions between s and rho, please see ?HPLB.")
    }
    ordering.type = "multiple"

  } else {
    if (length(rho) != length(t)) {
      stop("Incompatible dimensions between t and rho, please see ?HPLB.")
    }
    ordering.type = "single"
  }

  # ordering.type check
  if (estimator.type == "bayes" && ordering.type == "multiple") {
    stop("Invalid ordering type for binomial-test estimator, please see ?HPLB.")
  }

  # alpha check
  if (!is.numeric(alpha) || (length(alpha) != 1) || (alpha > 1) || (alpha <= 0)) {
    stop("Invalid type I error level alpha, please see ?HPLB.")
  }

  # tv.seq check
  if (!is.numeric(tv.seq) || any(!is.finite(tv.seq)) || any(tv.seq > 1) || any(tv.seq < 0) || any(duplicated(tv.seq))) {
    stop("Invalid grid for tv, please see ?HPLB.")
  }

  # direction check
  if ((estimator.type == "adapt_dwit") && (!is.character(direction) || !(direction %in% c("left", "right")) || (length(s) != length(direction)))) {
    stop("Invalid direction, please see ?HPLB.")
  }

  # check that with user-defined bounding sequence or empirical-tv-seach s should be uni-dim
  if (estimator.type %in% c("adapt_empirical", "adapt_custom")) {
    if (length(s) != 1) {
      stop("Incompatible value of s and estimator.type, please see ?HPLB.")
    }
    if (estimator.type == "adapt_custom") {
      if (is.null(custom.bounding.seq)) {
        stop("Invalid custom.bounding.seq, please see ?HPLB.")
      } else {
        if (length(custom.bounding.seq) != length(tv.seq)) {
          stop("Incompatible dimensions between custom.bounding.seq and tv.seq, please see ?HPLB.")
        }
      }
    }

    if (estimator.type %in% c("adapt_empirical")) {
      custom.bounding.seq <- empiricalBF(tv.seq = tv.seq,
                                         alpha = alpha,
                                         m = sum(t <= s),
                                         n = sum(t > s), ...)
    }

    # dimensions of s, k and tresholds check
    if (length(s) != length(k)) {
      stop("Incompatible dimensions between s, and k, please see ?HPLB.")
    }
  }


  ## core of the code

  # sort in increasing order the grid of tv value candidates
  tv.seq <- sort(tv.seq, decreasing = FALSE)

  # define the returned values
  dwithat.left  <- rep(NA, length(s))
  dwithat.right <- rep(NA, length(s))
  tvhat.left    <- rep(NA, length(s))
  tvhat.right   <- rep(NA, length(s))
  tvhat         <- rep(NA, length(s))

  # main loop, iteration over the different splits
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
      warning(paste0("m or n is zero for split ", s, ", please see ?HPLB."))
      next
    }


    ## binary classifer (finite sample version of bayes)
    if (estimator.type == "bayes_finite_sample") {

      # compute the decisions with randomization when rho is exactly 0.5
        decision <- ifelse(rho == 0.5,
                           sample(c(TRUE, FALSE),
                                  size = length(rho), replace = TRUE),
                           rho > cutoff)


      # compute the rightly classified number of points per class
      obs0 <- sum(!decision & t == 0)
      obs1 <- sum(decision & t == 1)

      # accuracy lower bound for class 0
      phat0 <- stats::binom.test(x = obs0,
                                 n = sum(t == 0),
                                 p = 0.5,
                                 alternative = "greater",
                                 conf.level = 1 - (alpha / 2))$conf.int[1]

      # accuracy lower bound for class 1
      phat1 <- stats::binom.test(x = obs1,
                                 n = sum(t == 1),
                                 p = 0.5,
                                 alternative = "greater",
                                 conf.level = 1 - (alpha / 2))$conf.int[1]

      # compute the estimated tvhat LB value
      tvhat <- max(phat0 + phat1 - 1, 0)

    }


    ## binary classifier (bayes from the paper)
    if (estimator.type == "bayes") {

      # compute the decisions
      decision <- ifelse(rho == 0.5,
                         sample(c(TRUE, FALSE),
                                size = length(rho), replace = TRUE),
                         rho > cutoff)


      # accuracy for class 0
      Ahat0 <- sum(!decision & t == 0) / sum(t==0)

      # accuracy for class 1
      Ahat1 <- sum(decision & t == 1) / sum(t==1)

      # compute the asymptotic standard deviation
      sigma <- sqrt(Ahat0 * (1-Ahat0) / sum(t==0)  +  Ahat1 * (1-Ahat1) / sum(t==1))

      # compute the tvhat value
      tvhat <- max(Ahat0 + Ahat1 - 1 -  stats::qnorm(1-alpha) * sigma, 0)

    }

    ## estimators based on adaptive binary classification (adapt)

    # define the V function corresponding to m and z (and n and z respectively)
    V.left  <- cumsum(ranks.table$t <= s[k])
    V.right <- cumsum(rev(ranks.table$t > s[k]))

    if (estimator.type == "adapt") {

      # recursive search
      max.stat <- 1
      step     <- 1

      while (max.stat > 0 ) {

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
          warning(paste0("m0 or n0 smaller or equal to 2 for split: ",
                         s[k], ",
                         asymtotic regime may not hold, see ?HPLB."))
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
          graphics::plot(V.left - mean0, type="l",
               main="Observed centered V-function with asymptotic bounding function", xlab="z", ylab="V(z)",
               font.main=1,font.lab=1,font.axis=1)
          graphics::lines(beta.left * sd0, col="blue")
        }

        if (tv.cur == max(tv.seq)){
          break
        }

     }

      # return the tv LB estimate
      tvhat[k] <- tv.cur

    } else if (estimator.type %in% c("adapt_empirical", "adapt_custom")) {

      # recursive search
      max.stat <- 1
      step     <- 1

      while (max.stat > 0 & step != length(tv.seq)) {

        tv.cur <- tv.seq[step]
        max.stat <- max(V.left - custom.bounding.seq[[step]])
        step <- step + 1

        if (verbose.plot && (tv.cur == 0)) {
          graphics::plot(V.left - mean0, type="l",
               main="Observed centered V-function with custom bounding function",xlab="z",ylab="V(z)",
               font.main=1,font.lab=1,font.axis=1)
          graphics::lines(beta.left * sd0, col="blue")
        }
      }

      # return the tv LB estimate
      tvhat[k] <- tv.cur

    } else if (estimator.type == "adapt_wit") {

      if (direction[k] == "left") {

        m <- nrow(ranks.table[t <= s[k]])
        n <- nrow(ranks.table) - m

      } else if (direction[k] == "right") {

        n <- nrow(ranks.table[t <= s[k]])
        m <- nrow(ranks.table) - n
      }

      # recursive search default values
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

  # return a final list with relevant LB estimates
  return(Filter(x = list(dwithat.left  = dwithat.left,
                         dwithat.right = dwithat.right,
                         tvhat.left    = tvhat.left,
                         tvhat.right   = tvhat.right,
                         tvhat         = tvhat),
                f = function(ll) !any(is.na(ll))))
}
