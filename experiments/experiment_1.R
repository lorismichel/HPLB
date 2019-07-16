# low detection rate
require(dWit)
grid.N <- c(100, 500, 1000, 5000, 10000, 50000, 100000)
gam <- c(-0.5, -0.7, -0.9)

tv.matrix <- matrix(nrow=length(grid.N), ncol=length(gam))

for (i in 1:length(grid.N)) {
  N <- grid.N[i]
  m <- N/2
  n <- N/2
  x1 <- ifelse(runif(m) > m^{gam[1]}, rexp(m), -rexp(m))
  x2 <- rexp(n)
  t <- c(rep(0, m), rep(1, n))
  b <- function(x) {dexp(x) / (dexp(x) + (1-m^{gam[1]})*dexp(x) + m^{gam[1]}*dexp(-x))}
  preds <- b(c(x1,x2))

  tv <- dWit(t = t, rho = preds, s = 0.5, estimator.type = "tv-search")$tvhat
  tv.matrix[i,1] <- tv
  print(tv.matrix)
}
