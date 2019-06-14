#' Lower-bound estimate of the number of distributional witness & total variation distance
#' for mixture models
#' @param times a numeric value corresponding to the times of the observations. For two-sample
#'   test times could be 0-1 values.
#' @param preds a numeric value giving the ordering of the observations. This could come either from
#'   a classification model or a regression model.
#' @param s the time value at which we consider a split and test for difference in distribution between the
#'   the mixture sitting at the left of s and the one at the right of s.
#' @param verbose.plot a logical value indicating if we want some summary plots.
#' @param alpha the level of the type 1 error control.
#' @param type should we use the base estimator?
#' @param use.correction.term do we want to use correction term?
#' @param tv.grid a grid of values for the search.
#' @import data.table
#' @export
dwlb <- function(times,
                 preds,
                 s,
                 verbose.plot = FALSE,
                 alpha = 0.05,
                 type = "global null",
                 use.correction.term = FALSE,
                 tv.grid = seq(0.05, 0.95, 0.01)) {

  # sort in decreasing order the grid of tv values
  tv.grid <- sort(tv.grid, decreasing = FALSE)

  # first we obtain the ranks from the predictions
  ranks <- rank(preds, ties.method = "random")

  # then we cast everything in a data.table object
  ranks.table <- data.table::data.table(rank = ranks, t = times)

  # do we plot the time against timehat scatter plot
  if (verbose.plot) {
    plot(times, preds, pch=19, xlab="t", ylab="that")
  }

  # let us order by rank the table
  data.table::setkey(ranks.table, "rank")

  TVhat_asymptoticFs<-c()
  TVhat_asymptoticGs<-c()
  lambdahat_asymptoticFs<-c()
  lambdahat_asymptoticGs<-c()
  i <- 0


  for (ss in s) {

    i <- i + 1

    if (i > 1) {
      print(paste(i,"out of", length(s)))
    }

    # getting the parameters m and n
    m <- nrow(ranks.table[t <= ss])
    n <- nrow(ranks.table) - m

    if (m == 0 || n==0) {
      stop("m or n is zero.")
    }

    # define the function function V_z corresponding to m and z
    V_zFs <- sapply(1:nrow(ranks.table), function(z) nrow(ranks.table[t <= ss & rank <= z]))
    #V_zGs <- sapply(1:nrow(ranks.table), function(z) nrow(ranks.table[t > ss & rank <= z]))
    V_zGs <- sapply(nrow(ranks.table):1, function(z) nrow(ranks.table[t > ss & rank >= z]))

    if (type == "global null") {
    # define the beta bounding sequence from the asymptotic result & the standard deviation
    # TODO: check these values
    x_alpha <- -log(-log(1-alpha)/2)
    betaFs <- sqrt(2 * log(log(m))) + (log(log(log(m)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(m))))
    betaGs <- sqrt(2 * log(log(n))) + (log(log(log(n)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(n))))

    sd0 <- sapply(1:nrow(ranks.table), function(z) sqrt(z * (m /(m+n)) * (n /(m+n)) * ((m+n-z)/(m+n-1))))

    # centered stats
    centered.statFs <- (V_zFs - sapply(1:nrow(ranks.table), function(z) (z*m)/(m+n)))
    centered.statGs <- (V_zGs - sapply(1:nrow(ranks.table), function(z) (z*n)/(m+n)))


    # summary plot
    if (verbose.plot) {
      plot(1:nrow(ranks.table),centered.statFs,type="l",ylab="V_z - (zm)/(m+n)",ylim=c(min(c(betaFs*sd0,centered.statFs)), max(c(betaFs*sd0,centered.statFs))),
           font.main=3,font.lab=3,
           xlab="z (order statistics of combined sample)",main="Observed Maximal Discrepancies")
      lines(1:nrow(ranks.table), betaFs*sd0,col="black",lty=2)
      lines(1:nrow(ranks.table), centered.statGs,col="brown")
      lines(1:nrow(ranks.table), betaGs*sd0,col="brown",lty=2)
    }


    # correction term applied
    lambdahat_asymptoticFs[i] <- max(max(centered.statFs - betaFs * sd0), 0) / (if (use.correction.term) (n/(m+n)) else 1)
    lambdahat_asymptoticGs[i] <- max(max(centered.statGs - betaGs * sd0), 0) / (if (use.correction.term) (m/(m+n)) else 1)

    TVhat_asymptoticFs[i]<-invertBinMeanTest(n.success = lambdahat_asymptoticFs[i],
                                             n.trial = m,alpha = alpha,
                                             rule.of.three = T)
    TVhat_asymptoticGs[i]<-invertBinMeanTest(n.success = lambdahat_asymptoticGs[i],
                                             n.trial = n,alpha = alpha,
                                             rule.of.three = T)
    } else {
      # recursive search
      max.stat <- 1
      step <- 1

      while(max.stat > 0) {

        tv.cur <- tv.grid[step]
        step <- step + 1

        lambda_left <- qbinom(p = 1-alpha/3, size = m, prob = tv.cur)
        lambda_right <- qbinom(p = 1-alpha/3, size = n, prob = tv.cur)
        m0 <- m - lambda_left
        n0 <- n - lambda_right

        # asymptotic values
        x_alpha <- -log(-log(1-alpha/3)/2)
        betaFs <- sqrt(2 * log(log(m0))) + (log(log(log(m0)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(m0))))
        betaGs <- sqrt(2 * log(log(n0))) + (log(log(log(n0)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(n0))))

        # mean and sd under the null
        sd0 <- sapply(1:(m0+n0), function(z) sqrt(z * (m0 /(m0+n0)) * (n0 /(m0+n0)) * ((m0+n0-z)/(m0+n0-1))))
        mean0 <- sapply(1:(m0+n0), function(z) z * (m0 /(m0+n0)))

        # creating the bands
        Vz_bands <- rep(0, m + n)

        if (lambda_left != 0) {
          Vz_bands[1:lambda_left] <- 1
        }
        Vz_bands <- cumsum(Vz_bands)

        if (lambda_left != 0) {
          ind.0 <- which(!c(1:(m+n)) %in% c(1:lambda_left, (m+n):(m+n-lambda_right+1)))
        } else {
          ind.0 <- 1:(m+n)
        }

        Vz_bands[ind.0] <- Vz_bands[ind.0] + mean0 + betaFs * sd0
        Vz_bands[(m+n):(m+n-lambda_right+1)] <- m
        Vz_bands <- pmin(Vz_bands, m)

        max.stat <- max(V_zFs-Vz_bands)
      }

      # updating the value for the split
      lambdahat_asymptoticFs[i] <- qbinom(p = tv.cur, size = m, prob = alpha)
      lambdahat_asymptoticGs[i] <- qbinom(p = tv.cur, size = n, prob = alpha)
      TVhat_asymptoticFs[i] <- tv.cur
      TVhat_asymptoticGs[i] <- tv.cur
    }
  }

  return(list(lambdahat_asymptoticFs = lambdahat_asymptoticFs,
              lambdahat_asymptoticGs = lambdahat_asymptoticGs,
              TVhat_asymptoticFs = TVhat_asymptoticFs,
              TVhat_asymptoticGs = TVhat_asymptoticGs))
}

