# Sequence of permutations interesting to study
permLabels <- function(labels,
                       preds = rep(0.5, length(labels)),
                       prob = 0.5,
                       u = 1,
                       d = 0,
                       perm.inside = TRUE) {

  if (perm.inside) {
    labels.perm <- ifelse(preds <= u & preds >= d,
                          ifelse(rbinom(n = length(labels), size = 1, prob = prob)==0,
                                 labels,
                                 1-labels), labels)
  } else {
    labels.perm <- ifelse(preds >= u | preds <= d,
                          ifelse(rbinom(n = length(labels), size = 1, prob = prob)==0,
                                 labels,
                                 1-labels), labels)
  }

  return(labels.perm)
}


