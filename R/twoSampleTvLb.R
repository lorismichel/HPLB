#' Two Sample Total Variance Lower Bound Estimate
#' @param labels a factor value with levels 0 and 1 encoding the membership to the two samples.
#' @param preds a numeric value between 0 and 1 giving the probability predictions for class 1.
#' @param alpha a value giving the type 1 error control.
#' @return the lower bound estimate
#' @export
twoSampleTvLb <- function(labels,
                          preds,
                          alpha = 0.05) {

  if (!is.factor(labels) || !(levels(labels) %in% c(0,1))) {
    stop("labels is not of correct format.")
  }
  if (!is.numeric(preds)) {
    stop("preds should be numeric and should lie in [0,1].")
  }
  if (length(labels) != length(preds)) {
    stop("incompatibles lengths of labels and preds.")
  }

  unique.labels <- unique(labels)
  tab.labels <- table(labels) / sum(table(labels))

  small.labels <- names(which.min(tab.labels))

  times <- ifelse(labels == small.labels, runif(length(labels), min = 0, max = min(tab.labels)),
                  runif(length(labels), min = min(tab.labels), max = 1))

  res <- dwlb(times = times, preds = if(small.labels=="0") preds else 1-preds, s = min(tab.labels), alpha = alpha, verbose.plot = FALSE)

  tvlb <- res$TVhat_asymptoticFs
  names(tvlb) <- c("total variation distance lower bound")
  return(tvlb)
}


