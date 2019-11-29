#' Distributional Witness Sampling
#'
#' @param x a numeric matrix of observations as rows.
#' @param f a function giving the density left.
#' @param g a function giving the density right.
#' @param seed a numeric value for reproducibility.
#'
#' @return a list with values:
#'   \itemize{-dwit.f, dwit.g the witness indicators left and right respectively.
#'            -pdwit.f, pdwit.g the corresponding probabilities.}
dWitSampling <- function(x,
                         f,
                         g,
                         seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # left
  u1 <- runif(nrow(x))
  dwit.f <- as.numeric(u1 <= (pdwit.f <- apply(x, 1, function(xx) (pmax(f(xx) - g(xx), 0) / f(xx)))))

  # right
  u2 <- runif(nrow(x))
  dwit.g <- as.numeric(u1 <= (pdwit.g <- apply(x, 1, function(xx) (pmax(g(xx) - f(xx), 0) / g(xx)))))

  return(list(dwit.f = dwit.f,
              dwit.g = dwit.g,
              pdwit.f = pdwit.f,
              pdwit.g = pdwit.g))
}
