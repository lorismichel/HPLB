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

~                                  
