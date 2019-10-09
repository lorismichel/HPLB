
# source
source("./climate-preprocessing.R")

# libs
require(dWit)

# prepro
d <- climatePrePro()

# analysis
dat <- data.frame(class=d$tenyears.ind,
                  x = cbind(t(d$air),t(d$mslp),t(d$prate),t(d$shum)))

# splits train-test
ind.train <- unlist(lapply(lapply(unique(dat$class), function(i) which(dat$class==i)), function(l) l[l > quantile(l,0.25) & l <= quantile(l,0.75)]))
ind.test <- unlist(lapply(lapply(unique(dat$class), function(i) which(dat$class==i)), function(l) l[l <= quantile(l,0.25) | l > quantile(l,0.75)]))

dat.train <- dat[ind.train,]
dat.test <- dat[ind.test,]

# source
source("./helpers-data.R")



"~/Downloads/data 2/steel-plat/"

datasets.multi <- c("annealing", "arrhythmia", "audiology-std", "balance-scale", "breast-tissue", "car",
                    "cardiotocography-10clases", "cardiotocography-3clases", "chess-krvk", "conn-bench-vowel-deterding",
                    "contrac", "dermatology", "ecoli",
                    "energy-y1", "energy-y2", "flags", "glass", "hayes-roth", "heart-cleveland", "heart-va",
                    "image-segmentation", "led-display", "lenses", "letter", "libras", "low-res-spect",
                    "lung-cancer", "lymphography", "molec-biol-splice",
                    "nursery", "oocytes_merluccius_states_2f", "oocytes_trisopterus_states_5b",
                    "optical", "page-blocks", "pendigits", "pittsburg-bridges-MATERIAL", "pittsburg-bridges-REL-L",
                    "pittsburg-bridges-SPAN", "pittsburg-bridges-TYPE", "plant-margin",
                    "plant-shape", "plant-texture", "primary-tumor", "seeds", "semeion",
                    "soybean", "st-image", "statlog-image", "statlog-landsat", "statlog-shuttle",
                    "st-vehicle", "steel-plates", "synthetic-control", "teaching", "thyroid", "vc-3classes",
                    "wall-following", "waveform", "waveform-noise", "wine", "wine-quality-red", "wine-quality-white",
                    "yeast", "zoo")


for (n in datasets.multi) {

  # make sure at least two points by class in test
  d <- getClassDataset(dataset = n, p.train = 0.5, PATH.UCI.DATA = "~/Downloads/data 2/")
  while(any(table(d$y.test)<=1)) {
    d <- getClassDataset(dataset = n, p.train = 0.5, PATH.UCI.DATA = "~/Downloads/data 2/")
  }
  # orderings
  ordering.array <- array(dim=c(length(levels(d$y.train)),length(levels(d$y.train)),nrow(d$x.test)))

  # fit a forest
  mRF <- ranger::ranger(class~., data = data.frame(d$x.train, class = d$y.train), probability = TRUE)

  # build the rho functions
  preds <- predict(mRF, data = data.frame(d$x.test))$predictions

  for (i in 1:length(levels(d$y.train))) {
    for (j in 1:length(levels(d$y.train))) {
      ordering.array[i,j,] <- preds[,j]-preds[,i]
    }
  }

  # get the tv mat
  tv.mat <- getTvLbDistanceMatrix(labels = d$y.test, ordering.array = ordering.array)
  print(tv.mat)
}











