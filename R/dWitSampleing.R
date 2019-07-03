#' Distributional witnesses sampling
#'@export
dWitSampling <- function(x, f, g, seed = 0) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  u <- runif(length(x))
  wF <- as.numeric(u <= (pmax(f(x)- g(x),0) / f(x)))
  u <- runif(length(x))
  wG <- as.numeric(u <= (pmax(g(x)- f(x),0) / g(x)))

  return(list(wF=wF, wG=wG,pF=(pmax(f(x)- g(x),0) / f(x)), pG=(pmax(g(x)- f(x),0) / g(x))))
}
