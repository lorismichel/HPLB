#' Inverting a Binomial Test
#' @param n.success number of successes.
#' @param n.trial number of trials.
#' @param alpha the level.
#' @param rule.of.three should we apply the rule of three?
invertBinMeanTest <- function(n.success,
                              n.trial,
                              alpha,
                              rule.of.three = TRUE) {

  if (n.success == 0) {
    if (rule.of.three)
      return(3/n.trial)
    else return(0)
  }

  phat <- n.success/n.trial
  by <- phat / max(100000, n.trial)
  p.candidates <- seq(0, 1, by = by)
  q <- sapply(p.candidates, FUN = function(y) qbinom(p = 1-alpha, size = n.trial, prob = y))

  return(p.candidates[tail(which(q <= n.success), 1)])
}
