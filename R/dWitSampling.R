#' Distributional Witnesses Sampling
#'
#' @param x a numeric matrix of observations as rows
#' @param f a function, the density for witnesses left
#' @param g a function, the density for witnesses right
#' @param seed a numeric value for reproducibility
#'
#' @return a list with values:
#'   \itemize{-wF, wG the witness indicators left and right respectively
#'            -pF, pG the corresponding probabilities}
#' @export
dWitSampling <- function(x,
                         f,
                         g,
                         seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # left
  u1 <- runif(nrow(x))
  wF <- as.numeric(u1 <= apply(x, 1, function(xx) (pmax(f(xx) - g(xx), 0) / f(xx))))

  # right
  u2 <- runif(nrow(x))
  wG <- as.numeric(u1 <= apply(x, 1, function(xx) (pmax(g(xx) - f(xx), 0) / g(xx))))

  return(list(wF = wF,
              wG = wG,
              pF = apply(x, 1, function(xx) (pmax(f(xx) - g(xx), 0) / f(xx))),
              pG = apply(x, 1, function(xx) (pmax(g(xx) - f(xx), 0) / g(xx)))))
}
