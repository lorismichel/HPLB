# SIMUlATION_5
# type: checking power by a loglog plot of sample size against contamination/tv
#       comparison with missclassification error rate test.

## saving paths
PATH.SAVE = "../Data/"

## requirements
require(dWit)
require(parallel)
require(data.table)
require(stats)
require(ranger)
require(MASS)
require(mlbench)
require(pROC)

# source
source("./helpers-data.R")


## permutation of rho
permuteRho <- function(rho, u, d) {
  ind <- which(rho <= u & rho >= d)
  rho.perm <- rho
  rho.perm[ind] <- sample(rho.perm[ind])
  return(rho.perm)
}

# core analysis
runPermAnalysis <- function(dataset="Boston") {


  # running simulations on clusters, params
  grid.p <- rep(seq(0, 0.5, by = 0.01), each=5)
  grid   <- grid.p


  # running simulations
  set.seed(123, "L'Ecuyer")
  res <- mcmapply(FUN = function(p) {

    d <- getDataset(dataset = dataset)
    y.train <- d$y.train
    x.train <- d$x.train
    y.test  <- d$y.test
    x.test  <- d$x.test

    # 1) fit a forest RF
    rf <- ranger::ranger(y~., data = data.frame(y = y.train, x.train),classification = TRUE, probability = TRUE)

    # 2) get predictions on test
    rho <- predict(rf, data = data.frame(x.test))$predictions[,"1"]

    # 3) get quantile predictions
    dq <- quantile(rho, p)
    uq <- quantile(rho, 1-p)

    # 4) permute labels
    rho.perm <- permuteRho(rho = rho, u = uq, d = dq)

    # 5) evaluate
    tvhat1 <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5,estimator.type = "asymptotic-tv-search")$tvhat
    tvhat2 <- dWit(t = as.numeric(levels(y.test))[y.test], rho = rho.perm, s = 0.5, estimator.type = "hypergeometric-test", z = sum(y.test==0),  threshold = mean(rho.perm))$tvhat
    auc <- auc(as.numeric(levels(y.test))[y.test], rho.perm)
    c(tvhat1, tvhat2, auc)}, grid, mc.cores = 10)


# gathering data
 return(power.table <- data.table(dataset        = dataset,
                                  p              = grid,
                                  tvhat_search   = res[1,],
                                  tvhat_binomial = res[2,],
                                  auc            = res[3,],
                                  auc_search     = 0.5 + 0.5*res[1,]^2,
                                  auc_binomial   = 0.5 + 0.5*res[2,]^2))
}

# dataset names
dataset.names <- c("Boston", "titanic", "BreastCancer", "Ionosphere", "abalone",# "adult",
                   "banknotes", "Default", "credit")

two.classes.datasets <- c("acute-inflammation", "acute-nephritis", "balloons", "bank", "blood",
                          "breast-cancer", "breast-cancer-wisc", "breast-cancer-wisc-diag", "breast-cancer-wisc-prog", "chess-krvkp", "congressional-voting",
                          "conn-bench-sonar-mines-rocks", "connect-4", "credit-approval", "cylinder-bands", "echocardiogram", "fertility", "haberman-survival",
                          "heart-hungarian", "heart-switzerland", "hepatitis", "hill-valley", "horse-colic", "ilpd-indian-liver", "ionosphere",
                          "magic", "mammographic", "miniboone", "molec-biol-promoter", "monks-1", "monks-2", "monks-3", "mushroom", "musk-1", "musk-2",
                          "oocytes_merluccius_nucleus_4d", "oocytes_merluccius_states_2f", "ozone", "parkinsons", "pima", "pittsburg-bridges-T-OR-D",
                          "planning", "ringnorm", "spambase", "spect", "spectf", "statlog-australian-credit", "statlog-german-credit", "statlog-heart",
                          "tic-tac-toe", "titanic", "trains", "twonorm", "vc-2classes")

# running the sims
res <- data.table()
for (n in two.classes.datasets) {
  res <- rbind(res, runPermAnalysis(dataset = n))
  print(n)
}



save(res, file = paste0(PATH.SAVE, "DATA_SIMULATION_5.Rdata"))
