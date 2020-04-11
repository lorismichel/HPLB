#' Bounding Operation
#'
#' @param v a a numeric value giving an ordering permutation of 1 to m+n.
#' @param left a numeric value giving the number of witnesses left.
#' @param right a numeric value giving the number of witnesses right.
#' @param m a numeric value, the number of observations left.
#' @param n a numeric value, the number of observations right.
#'
#' @return a cumulative counting function represented as a numeric vector.
boundingOperation <- function(v,
                              left,
                              right,
                              m,
                              n) {
  l <- length(v)

  if (left != 0) {
    v <- setdiff(v,c(1:left))
  }
  if (right != 0) {
    v <- setdiff(v,c((l-right+1):l))
  }

  return(cumsum(c(rep(1, left), c(rep(1,m),rep(0, n))[v], rep(0, right))))
}



#' Empirical Bounding Functions
#'
#' @param tv.seq a vector of total variation values between 0 and 1.
#' @param nsim a numeric value giving the number of repetitions.
#' @param m a numeric value, the number of observations left.
#' @param n a numeric value, the number of observations right.
#' @param alpha a numeric value giving the type-I error level.
#' @return a list of empirical bounding functions indexed by the tv.seq (in the respective order).
#'
#' @export
empiricalBF <- function(tv.seq,
                        nsim  = 1000,
                        m     = 100,
                        n     = 100,
                        alpha = 0.05) {

  l <- lapply(1:nsim, function(i) {v <- sample(1:(m+n)); lapply(tv.seq, function(tv) boundingOperation(v = v,
                                                                     left = stats::qbinom(prob = tv, size = m, p = 1-alpha/2),
                                                                     right = stats::qbinom(prob = tv, size = n, p = 1-alpha/2),
                                                                     m = m, n = n))})

  ll <- Reduce(x = l, f = function(x, y) mapply(x, y, FUN = function(xx,yy) pmax(xx,yy), SIMPLIFY = FALSE))

  return(ll)
}
