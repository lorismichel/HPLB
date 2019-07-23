#' Lower-bound estimate of the number of distributional witness & total variation distance
#' for mixture models
#' @param t a numeric value corresponding to an ordering of the observations. For two-sample
#'   test times could be 0-1 numeric values values.
#' @param rho a numeric value giving an ordering function. This could be
#'   a classifier, a regressor, a witness function from an MMD kernel or anything else.
#' @param s a numeric value giving split points on t.
#' @param alpha the level of the type 1 error control.
#' @param estimator.type either "basic", "tv-search", "wit-search", "binomial"
#' @param tv.grid a grid of values for the search.
#' @param direction which witness to estimate.
#' @param threshold this is the threshold used if missclassification error in used.
#' @import data.table
#' @export
dWit <- function(t,
                 rho,
                 s,
                 alpha              = 0.05,
                 estimator.type     = "basic",
                 tv.grid            = seq(from = 0,
                                          to   = 1,
                                          by   = 0.01),
                 direction          = rep("left", length(s)),
                 threshold          = 0.5) {

  # parameters checks

  # test for t
  if (!is.numeric(t) || any(!is.finite(t)) || is.matrix(t)) {
    stop("Invalid values or type for t, please see ?dWit.")
  }

  # test for s
  if (!is.numeric(s) || any(!is.finite(s)) || any(s > max(t)) || any(s < min(t)) || is.matrix(s)) {
    stop("Invalid values or type for s, please see ?dWit.")
  }

  # test for rho
  if (!is.numeric(rho) || any(!is.finite(rho))) {
    stop("Invalid values for rho, please see ?dWit.")
  }

  # two cases
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

  # test for the estimator type
  if (!is.character(estimator.type) || !(estimator.type %in% c("basic", "tv-search", "wit-search", "binomial"))) {
    stop("Invalid estimator type, please see ?dWit.")
  }

  # test for the type I error level
  if (!is.numeric(alpha) || (alpha > 1) || (alpha <= 0)) {
    stop("Invalid type I error level alpha, please see ?dWit.")
  }

  # test for the tv.grid
  if (!is.numeric(tv.grid) || any(tv.grid > 1) || any(tv.grid < 0) || any(duplicated(tv.grid))) {
    stop("Invalid grid for tv, please see ?dWit.")
  }

  # check direction
  if (!is.character(direction) || !(direction %in% c("left", "right"))) {
    stop("Invalid direction, please see ?dWit.")
  }

  # sort in increasing order the grid of tv values
  tv.grid <- sort(tv.grid, decreasing = FALSE)

  # define the returned values
  lambdahat.Fs <- rep(NA, length(s))
  lambdahat.Gs <- rep(NA, length(s))
  tvhat.Fs <- rep(NA, length(s))
  tvhat.Gs <- rep(NA, length(s))
  tvhat <- rep(NA, length(s))

  # main loop over splits
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
      stop("m or n is zero, please see ?dWit.")
    }

    # define the function function V_z corresponding to m and z
    V_zFs <- cumsum(ranks.table$t <= s[k])
    V_zGs <- cumsum(rev(ranks.table$t > s[k]))
    #V_zFs <- sapply(1:nrow(ranks.table), function(z) nrow(ranks.table[t <= s[k] & rank <= z]))
    #V_zGs <- sapply(nrow(ranks.table):1, function(z) nrow(ranks.table[t > s[k] & rank >= z]))

    # branching on estimators
    if (estimator.type == "basic") {

      # define the beta bounding sequence from the asymptotic result & the standard deviation
      x_alpha <- -log(-log(1-alpha/2)/2)
      betaFs <- sqrt(2 * log(log(m))) + (log(log(log(m)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(m))))
      betaGs <- sqrt(2 * log(log(n))) + (log(log(log(n)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(n))))

      sd0 <- sapply(1:nrow(ranks.table), function(z) sqrt(z * (m /(m+n)) * (n /(m+n)) * ((m+n-z)/(m+n-1))))

      # centered stats
      centered.statFs <- (V_zFs - sapply(1:nrow(ranks.table), function(z) (z*m)/(m+n)))
      centered.statGs <- (V_zGs - sapply(1:nrow(ranks.table), function(z) (z*n)/(m+n)))

      # estimators construction
      lambdahat.Fs[k] <- max(max(centered.statFs - betaFs * sd0), 0)
      lambdahat.Gs[k] <- max(max(centered.statGs - betaGs * sd0), 0)

      tvhat.Fs[k] <- invertBinMeanTest(n.success     = lambdahat.Fs[k],
                                       n.trial       = m,
                                       alpha         = alpha/2,
                                       rule.of.three = FALSE)
      tvhat.Gs[k] <- invertBinMeanTest(n.success     = lambdahat.Gs[k],
                                       n.trial       = n,
                                       alpha         = alpha/2,
                                       rule.of.three = FALSE)

    } else if (estimator.type == "wit-search") {

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
      wit.left.grid <- c(0:m)

      while (max.stat > 0) {

        # 1) get a candidate nb witness left
        lambda.left <- wit.left.grid[step]

        # 2) get an estimate of tv by inverting binomial test
        tv.cur <- invertBinMeanTest(n.success     = lambda.left,
                                    n.trial       = m,
                                    alpha         = 1 - (alpha / 3),
                                    rule.of.three = FALSE)
        step <- step + 1

        # 3) take a high quantile of binomial with tv
        lambda.right <- qbinom(p = 1-alpha/3, size = n, prob = tv.cur)

        # define mo and no
        m0 <- m - lambda.left
        n0 <- n - lambda.right

        # asymptotic values
        x_alpha <- -log(-log(1-alpha/3)/2)
        betaFs <- sqrt(2 * log(log(m0))) + (log(log(log(m0)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(m0))))
        betaGs <- sqrt(2 * log(log(n0))) + (log(log(log(n0)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(n0))))

        # mean and sd under the null
        sd0 <- sapply(1:(m0+n0), function(z) sqrt(z * (m0 /(m0+n0)) * (n0 /(m0+n0)) * ((m0+n0-z)/(m0+n0-1))))
        mean0 <- sapply(1:(m0+n0), function(z) z * (m0 /(m0+n0)))

        # creating the quantile bands
        Q <- rep(0, m + n)

        # adding witnesses left
        if (lambda.left != 0) {
          Q[1:lambda.left] <- 1
          Q <- cumsum(Q)
        }

        # hypergeom filling
        ind.hypergeom <- c((lambda.left+1):(m+n-lambda.right))
        Q[ind.hypergeom] <- Q[ind.hypergeom] + mean0 + betaFs * sd0

        # witness right filling
        Q[(m+n-lambda.right):(m+n)] <- m
        Q <- pmin(Q, m)

        if (direction[k] == "left") {
          max.stat <- max(V_zFs - Q)
        } else if (direction[k] == "right") {
          max.stat <- max(V_zGs - Q)
        }

      }

      # updating the value for the split
      if (direction[k] == "left") {
        lambdahat.Fs[k] <- lambda.left
      } else if (direction[k] == "right") {
        lambdahat.Gs[k] <- lambda.left
      }
    } else if (estimator.type == "tv-search") {


      # recursive search
      max.stat <- 1
      step <- 1

      while (max.stat > 0) {

        tv.cur <- tv.grid[step]
        step <- step + 1

        lambda.left <- qbinom(p = 1-alpha/3, size = m, prob = tv.cur)
        lambda.right <- qbinom(p = 1-alpha/3, size = n, prob = tv.cur)

        m0 <- m - lambda.left
        n0 <- n - lambda.right

        # asymptotic values
        x_alpha <- -log(-log(1-alpha/3)/2)
        betaFs <- sqrt(2 * log(log(m0))) + (log(log(log(m0)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(m0))))
        betaGs <- sqrt(2 * log(log(n0))) + (log(log(log(n0)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(n0))))

        # mean and sd under the null
        sd0 <- sapply(1:(m0+n0), function(z) sqrt(z * (m0 /(m0+n0)) * (n0 /(m0+n0)) * ((m0+n0-z)/(m0+n0-1))))
        mean0 <- sapply(1:(m0+n0), function(z) z * (m0 /(m0+n0)))

        # creating the bands
        Q <- rep(0, m + n)

        if (lambda.left != 0) {
          Q[1:lambda.left] <- 1
          Q <- cumsum(Q)
        }

        # hypergeom filling
        ind.hypergeom <- c((lambda.left+1):(m+n-lambda.right))
        Q[ind.hypergeom] <- Q[ind.hypergeom] + mean0 + betaFs * sd0

        # witness right filling
        Q[(m+n-lambda.right):(m+n)] <- m
        Q <- pmin(Q, m)

        max.stat <- max(V_zFs-Q)
      }

      # return the tv estimate
      tvhat[k] <- tv.cur
    } else if (estimator.type == "binomial") {
      obs <- sum((rho > threshold & t == 0) | (rho <= threshold & t == 1))
      phat <- invertBinMeanTest(n.success = obs, n.trial = length(t),
                                alpha = 1-alpha,
                                rule.of.three = FALSE)
      tvhat <- max(1-2*phat,0)
    }
  }

  # return list
  return(list(lambdahat.Fs = lambdahat.Fs,
              lambdahat.Gs = rev(lambdahat.Gs),
              tvhat.Fs     = tvhat.Fs,
              tvhat.Gs     = tvhat.Gs,
              tvhat        = tvhat))
}
