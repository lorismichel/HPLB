#' Bounding operation
boundingOperation <- function(v, left, right, m, n) {
  l <- length(v)
  if (left != 0) {
    v <- setdiff(v,c(1:left))
  }
  if (right != 0) {
    v <- setdiff(v,c((l-right+1):l))
  }
  return(cumsum(c(rep(1, left), c(rep(1,m),rep(0, n))[v], rep(0, right))))
}



#'obtain empirical bounding functions
empiricalBF <- function(tv.seq, nrep = 1000, m = 100, n = 100, alpha = 0.05) {


  l <- lapply(1:nrep, function(i) {v <- sample(1:(m+n)); lapply(tv.seq, function(tv) boundingOperation(v = v,
                                                                     left = qbinom(prob = tv, size = m, p = 1-alpha/2),
                                                                     right = qbinom(prob = tv, size = n, p = 1-alpha/2),
                                                                     m = m, n = n))})
  ll <- Reduce(x = l, f = function(x, y) mapply(x,y,FUN = function(xx,yy) pmax(xx,yy), SIMPLIFY = FALSE))
}
