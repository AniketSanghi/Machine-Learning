import utils
import predict
import time as tm
import numpy as np
from sklearn import model_selection

d = 16385
L = 3400
(X, y) = utils.load_data( "data", d, L )

X_train, X_test, Y_train, Y_test = model_selection.train_test_split(X, y, test_size=0.20, random_state=100)

utils.dump_data(X_train, Y_train, "train")
utils.dump_data(X_test, Y_test, "test")