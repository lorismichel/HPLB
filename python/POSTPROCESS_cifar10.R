library(igraph)
gg <- graph.ring(4)
ll =matrix(c(0,0,0,1,0,3,0,5),ncol=2,byrow=TRUE)
plot(gg,layout=ll)


nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
load("~/phdWork/dWit/python/distance.Rdata")
dist_m <- mat
dist_m
dist_mi <- 1/dist_m
library(qgraph)
qgraph(dist_mi, layout = "spring",vsize=10)
attr(mat, 'Labels') <- colnames(mat)

"coldiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
{
  require(gclus)

  if (max(D)>1) D <- D/max(D)

  if (byrank) {
    spe.color = dmat.color(1-D, cm.colors(nc))
  }
  else {
    spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(nc))
  }

  spe.o = order.single(1-D)
  speo.color = spe.color[spe.o,spe.o]

  op = par(mfrow=c(1,1), pty="s")

  if (diag) {
    plotcolors(spe.color, rlabels=attributes(D)$Labels,
               main="Dissimilarity Matrix",
               dlabels=attributes(D)$Labels)
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o],
               main="Ordered Dissimilarity Matrix",
               dlabels=attributes(D)$Labels[spe.o])
  }
  else {
   # plotcolors(spe.color, rlabels=attributes(D)$Labels,
   #            main="Dissimilarity Matrix")
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o],
               dlabels=attributes(D)$Labels[spe.o],
               #clabels=attributes(D)$Labels[spe.o],
               main="Total Variation Lower-bound Matrix Between Classes")
  }

  par(op)
}


coldiss(mat)
