## Load data
require(reticulate)


py_run_string("
import pickle
with open('./cifar10/predictions_cifar10.pkl', 'rb') as f:
    data = pickle.load(f)
    y = data['labels']
    preds = data['predictions']")

print(table(py$y))
print(dim(py$preds))

ordering.array <- array(dim = c(10,10,nrow(py$preds)))
for (i in 1:10) {
  for (j in 1:10) {
    ordering.array[i,j,] <- py$preds[,j]-py$preds[,i]
  }
}

mat <- dWit::getTvLbDistanceMatrix(labels = py$y, ordering.array = ordering.array)
rownames(mat) <- colnames(mat) <- c('plane', 'car', 'bird', 'cat',
           'deer', 'dog', 'frog', 'horse', 'ship', 'truck')
print(mat)

save(mat, file="./distance.Rdata")
