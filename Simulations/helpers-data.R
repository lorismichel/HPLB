getClassDataset <- function(dataset="Boston",
                       p.train = 0.5,
                       PATH.UCI.DATA = "~/Downloads/data 2/") {

  if (dataset=="Boston") {

    # libs
    require(MASS)

    # define the labels
    crime01 = rep(0, length(Boston$crim))
    crime01[Boston$crim > median(Boston$crim)] = 1

    # sample
    Boston = data.frame(Boston)
    Boston <- Boston[sample(1:nrow(Boston), size = nrow(Boston), replace = FALSE),]

    # half train set size
    train = 1:(dim(Boston)[1]*p.train)
    test = setdiff(1:dim(Boston)[1], train)

    Boston.train = Boston[train, -c(1)]
    Boston.test = Boston[test, -c(1)]

    crim01.test = crime01[test]
    crim01.train = crime01[train]

    crim01.test = factor(rep(crim01.test,each=1))

    crim01.train = factor(rep(crim01.train,each=1))

    return(list(x.train = Boston.train, x.test = Boston.test,
                y.train = crim01.train, y.test = crim01.test))

  } else if (dataset == "Titanic") {

    titanic <- read.csv(url("https://web.stanford.edu/class/archive/cs/cs109/cs109.1166/stuff/titanic.csv"),
                        header = T)

    # labels
    survived = titanic$Survived

    # sample
    titanic = data.frame(titanic)
    titanic <- titanic[sample(1:nrow(titanic), replace = FALSE, size = nrow(titanic)),]

    # preprocess the labels
    titanic$Sex <- ifelse(titanic$Sex=="male",0,1)
    titanic <- titanic[,-c(3)]

    train = 1:(dim(titanic)[1]*p.train)
    test = setdiff(1:dim(titanic)[1], train)

    titanic.train = titanic[train, -c(1)]
    titanic.test = titanic[test, -c(1)]

    survived.train = factor(survived[train])
    survived.test = factor(survived[test])

    return(list(x.train = titanic.train, x.test = titanic.test,
                y.train = survived.train, y.test = survived.test))
  } else if (dataset == "BreastCancer") {

    # libs
    require(mlbench)

    # remove na
    data("BreastCancer")
    breast <- na.omit(BreastCancer)

    # sample
    breast <- breast[sample(1:nrow(breast), replace = FALSE, size = nrow(breast)),]

    # labels
    response <- factor(breast$Class)
    levels(response) <- c(0,1)

    train = 1:(dim(breast)[1]*p.train)
    test = setdiff(1:dim(breast)[1], train)

    breast.train = breast[train, -c(1,11)]
    breast.test = breast[test, -c(1,11)]

    response.train = factor(response[train])
    response.test = factor(response[test])

    return(list(x.train = breast.train, x.test = breast.test,
                y.train = response.train, y.test = response.test))

  } else if (dataset == "Ionosphere") {

    # libs
    require(mlbench)

    # remove na
    data("Ionosphere")
    iono <- na.omit(Ionosphere)

    # labels
    response <- factor(iono$Class)
    levels(response) <- c(0,1)

    # sample
    iono <- iono[sample(1:nrow(iono), replace = FALSE, size = nrow(iono)),]


    train = 1:(dim(iono)[1]*p.train)
    test = setdiff(1:dim(iono)[1], train)

    iono.train = iono[train, -c(2,35)]
    iono.test = iono[test, -c(2,35)]

    response.test = factor(response[test])
    response.train = factor(response[train])

    return(list(x.train = iono.train, x.test = iono.test,
                y.train = response.train, y.test = response.test))

  } else if (dataset == "Abalone") {

    abalone.cols = c("sex", "length", "diameter", "height", "whole.wt",
                     "shucked.wt", "viscera.wt", "shell.wt", "rings")

    url <- 'http://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data'
    abalone <- read.table(url, sep=",", row.names=NULL, col.names=abalone.cols,
                          nrows=4177)
    colnames(abalone) <- abalone.cols

    # sample
    abalone <- abalone[sample(1:nrow(abalone), replace = FALSE, size = nrow(abalone)),]


    train = 1:(dim(abalone)[1]*p.train)
    test = setdiff(1:dim(abalone)[1], train)
    abalone.train = abalone[train, -c(9)]
    abalone.test = abalone[test, -c(9)]

    # young and adult vs old, labels
    response <- ifelse(abalone$rings <=5, 0, ifelse(abalone$rings <= 13, 0, 1))
    response.test = factor(response[test])
    response.train = factor(response[train])

    return(list(x.train = abalone.train, x.test = abalone.test,
                y.train = response.train, y.test = response.test))

  } else if (dataset == "Banknotes") {

    banknotes <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/00267/data_banknote_authentication.txt',
                            sep = ',', fill = F, strip.white = T)

    # sample
    banknotes <- banknotes[sample(1:nrow(banknotes), size = nrow(banknotes), replace = FALSE),]
    colnames(banknotes) <- c("X1","X2","X3","X4","y")
    banknotes <- banknotes[sample(1:nrow(banknotes), replace = FALSE, size = nrow(banknotes)),]


    train = 1:(dim(banknotes)[1]*p.train)
    test = setdiff(1:dim(banknotes)[1], train)
    banknotes.train = banknotes[train, -c(5)]
    banknotes.test = banknotes[test, -c(5)]

    response <- banknotes$y
    response.test = factor(response[test])
    response.train = factor(response[train])

    return(list(x.train = banknotes.train, x.test = banknotes.test,
                y.train = response.train, y.test = response.test))

  } else if (dataset == "Default") {

    require(ISLR)
    default <- ISLR::Default

    # sample
    default <- default[sample(1:nrow(default), replace = FALSE, size = nrow(default)),]

    train = 1:(dim(default)[1]*p.train)
    test = setdiff(1:dim(default)[1], train)
    default.train = default[train, -c(1)]
    default.test = default[test, -c(1)]

    # labels
    response <- ifelse(default$default=="Yes",1,0)
    response.test = factor(response[test])
    response.train = factor(response[train])


    return(list(x.train = default.train, x.test = default.test,
                y.train = response.train, y.test = response.test))

  } else if (dataset == "Credit") {

    # credit
    url="http://freakonometrics.free.fr/german_credit.csv"
    credit=read.csv(url, header = TRUE, sep = ",")
    credit <- credit[sample(1:nrow(credit), replace = FALSE, size = nrow(credit)),]

    train = 1:(dim(credit)[1]*p.train)
    test = setdiff(1:dim(credit)[1], train)
    credit.train = credit[train, -c(1)]
    credit.test = credit[test, -c(1)]

    # labels
    response <- credit$Creditability
    response.test = factor(response[test])
    response.train = factor(response[train])



    return(list(x.train = credit.train, x.test = credit.test,
                y.train = response.train, y.test = response.test))

  } else if (dataset == "Adult") {
     adult <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data',
                         sep = ',', fill = F, strip.white = T)
     colnames(adult) <- c('age', 'workclass', 'fnlwgt', 'educatoin',
                          'educatoin_num', 'marital_status', 'occupation', 'relationship', 'race', 'sex',
                          'capital_gain', 'capital_loss', 'hours_per_week', 'native_country', 'income')
     levels(adult$income) <- c(0,1)
     adult <- adult[sample(1:nrow(adult), replace = FALSE, size = nrow(adult)),]


     train = 1:(dim(adult)[1]*p.train)
     test = setdiff(1:dim(adult)[1], train)
     adult.train = adult[train, -c(15)]
     adult.test = adult[test, -c(15)]

     response <- adult$income
     response.test = factor(response[test])
     response.train = factor(response[train])

     return(list(x.train = adult.train, x.test = adult.test,
                 y.train = response.train, y.test = response.test))
  } else if (dataset == "ringnorm2000") {

    d <- mlbench:::mlbench.ringnorm(n = 2000)

    return(list(x.train = d$x[1:1000,], x.test = d$x[1001:2000,],
                y.train = d$labels[1:1000,], y.test = d$labels[1001:2000,]))

  } else if (dataset == "ringnorm500") {

    d <- mlbench:::mlbench.ringnorm(n = 500)

    return(list(x.train = d$x[1:250,], x.test = d$x[251:500,],
                y.train = d$labels[1:250,], y.test = d$labels[251:500,]))

  } else if (dataset == "ringnorm10000") {

    d <- mlbench:::mlbench.ringnorm(n = 10000)

    return(list(x.train = d$x[1:5000,], x.test = d$x[5001:10000,],
                y.train = d$labels[1:5000,], y.test = d$labels[5001:10000,]))

  } else if (dataset %in% c("annealing")) {

    d.train <- read.table(file = paste0(PATH.UCI.DATA, dataset, "/", dataset, "_train_R.dat"))
    d.test <- read.table(file = paste0(PATH.UCI.DATA, dataset, "/", dataset, "_test_R.dat"))


    return(list(x.train = d.train[,-ncol(d.train)], x.test = d.test[,-ncol(d.test)],
                y.train = factor(d.train[,ncol(d.train)]), y.test = factor(d.test[,ncol(d.test)])))

  } else {
      d <- read.table(file = paste0(PATH.UCI.DATA, dataset, "/", dataset, "_R.dat"))

      ids.train <- sample(1:nrow(d), size = p.train * nrow(d), replace = FALSE)

      return(list(x.train = d[ids.train,-ncol(d)], x.test = d[-ids.train,-ncol(d)],
                  y.train = factor(d[ids.train,ncol(d)]), y.test = factor(d[-ids.train,ncol(d)])))
  }
}



# two.classes.datasets <- c("acute-inflammation", "acute-nephritis", "balloons", "bank", "blood",
#                             "breast-cancer", "breast-cancer-wisc", "breast-cancer-wisc-diag", "breast-cancer-wisc-prog", "chess-krvkp", "congressional-voting",
#                             "conn-bench-sonar-mines-rocks", "connect-4", "credit-approval", "cylinder-bands", "echocardiogram", "fertility", "haberman-survival",
#                             "heart-hungarian", "heart-switzerland", "hepatitis", "hill-valley", "horse-colic", "ilpd-indian-liver", "ionosphere",
#                             "magic", "mammographic", "miniboone", "molec-biol-promoter", "monks-1", "monks-2", "monks-3", "mushroom", "musk-1", "musk-2",
#                             "oocytes_merluccius_nucleus_4d", "oocytes_merluccius_states_2f", "ozone", "parkinsons", "pima", "pittsburg-bridges-T-OR-D",
#                             "planning", "ringnorm", "spambase", "spect", "spectf", "statlog-australian-credit", "statlog-german-credit", "statlog-heart",
#                             "tic-tac-toe", "titanic", "trains", "twonorm", "vc-2classes")
