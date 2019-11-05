import numpy as np
from numpy import random as rand
import os
from sklearn.datasets import load_svmlight_file
from sklearn.datasets import dump_svmlight_file
from scipy import sparse as sps

# DO NOT CHANGE THE NAME OF THIS METHOD OR ITS INPUT OUTPUT BEHAVIOR

# INPUT CONVENTION
# X: n x d matrix in csr_matrix format containing d-dim (sparse) features for n test data points
# k: the number of recommendations to return for each test data point in ranked order

# OUTPUT CONVENTION
# The method must return an n x k numpy nd-array (not numpy matrix or scipy matrix) of labels with the i-th row 
# containing k labels which it thinks are most appropriate for the i-th test point. Labels must be returned in 
# ranked order i.e. the label yPred[i][0] must be considered most appropriate followed by yPred[i][1] and so on

# CAUTION: Make sure that you return (yPred below) an n x k numpy (nd) array and not a numpy/scipy/sparse matrix
# The returned matrix will always be a dense matrix and it terribly slows things down to store it in csr format
# The evaluation code may misbehave and give unexpected results if an nd-array is not returned

def load_data( filename, L ):
    # X, _ = load_svmlight_file( "%s.X" % filename, multilabel = True, n_features = d, offset = 1 )
    y, _ = load_svmlight_file( filename, multilabel = True, n_features = L, offset = 1 )
    return y

def dump_data( X, filename ):
    (n, d) = X.shape
    dummy = sps.csr_matrix( (n, 1) )
    dump_svmlight_file( X, dummy, filename, multilabel = True, zero_based = True, comment = "%d %d" % (n, d) )

def getReco( X, k ):
    # Find out how many data points we have
    n = X.shape[0]
    # Load and unpack the dummy model
    # The dummy model simply stores the labels in decreasing order of their popularity
    # npzModel = np.load( "model.npz" )
    # model = npzModel[npzModel.files[0]]
    # Let us predict a random subset of the 2k most popular labels no matter what the test point
    # shortList = model[0:2*k]
    # Make sure we are returning a numpy nd-array and not a numpy matrix or a scipy sparse matrix
    dump_data(X, "test.txt")

    os.system('./fastXML_predict test.txt score.txt model')

    yPred = np.zeros( (n, k) )

    # Open the score file so created
    i = -1
    file = open("score.txt", "r")

    # Read File
    for line in file:
        # Ignore the first line
        if i==-1:
            i += 1
            continue

        # Remove the collons in the data
        Line = [x.split(":") for x in str(line).split(" ")]

        # Sort it based on the ranking so received
        Line.sort(reverse = True, key = lambda x: x[1])

        # Update our np.array
        yPred[i,:] = [int(x[0]) for x in Line[:k]]

        i += 1
  
    file.close()

    return yPred