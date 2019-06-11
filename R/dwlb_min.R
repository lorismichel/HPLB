#' Lower-bound estimate of the number of distributional witness & total variation distance
#' for mixture models or two-sample tests
#' @param times a numeric value corresponding to the times of the observations. For two-sample
#'   test times could be 0-1 values.
#' @param preds a numeric value giving the ordering of the observations. This could come either from
#'   a classification model or a regression model.
#' @param s the time value at which we consider a split and test for difference in distribution between the
#'   the mixture sitting at the left of s and the one at the right of s.
#' @param verbose.plot a logical value indicating if we want some summary plots.
#' @param alpha the level of the type 1 error control.
#' @import data.table
#' @export
dwlb_min <- function(times,
                 preds,
                 s,
                 verbose.plot = FALSE,
                 alpha = 0.05) {

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
  lambdahat_asymptoticFs <- rep(2, length(s))

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

    # define the beta bounding sequence from the asymptotic result & the standard deviation
    # TODO: check these values
    m.loc <- m
    while (lambdahat_asymptoticFs[i] > 0) {
      x_alpha <- -log(-log(1-alpha)/2)
      betaFs <- sqrt(2 * log(log(m.loc))) + (log(log(log(m.loc)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(m.loc))))
      centered.statFs <- (V_zFs - sapply(1:nrow(ranks.table), function(z) (z*m.loc)/(m.loc+n)))
      sd0 <- sapply(1:nrow(ranks.table), function(z) sqrt(z * (m.loc /(m+n)) * (n /(m.loc+n)) * ((m.loc+n-z)/(m.loc+n-1))))
      lambdahat_asymptoticFs[i] <- max(max(centered.statFs - betaFs * sd0), 0)
      m.loc <- m.loc - 1
    }

    #betaGs <- sqrt(2 * log(log(n))) + (log(log(log(n)))-log(pi)+2*x_alpha) / (2 * sqrt(2 * log(log(n))))

    TVhat_asymptoticFs[i]<-invertBinMeanTest(n.success = lambdahat_asymptoticFs[i],
                                             n.trial = m+n,alpha = alpha,
                                             rule.of.three = T)/ss
  }

  return(list(lambdahat_asymptoticFs = lambdahat_asymptoticFs,
              TVhat_asymptoticFs = TVhat_asymptoticFs))
}

