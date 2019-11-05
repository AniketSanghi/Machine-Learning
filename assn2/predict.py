import numpy as np
from numpy import random as rand
import os

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
    os.system('../Tree_Extreme_Classifiers/Parabel/parabel_predict dataX.txt ../Tree_Extreme_Classifiers/Parabel/model score.txt')

    yPred = np.zeros( (n, k) )

    i = -1
    f = open("score.txt", "r")
    for temp in f:
        if i==-1:
            i += 1
            continue
        lol = [x.split(":") for x in str(temp).split(" ")]
        lol.sort(reverse = True, key = lambda x: x[1])
        lol = [int(x[0]) for x in lol]
        yPred[i,:] = lol[0:5]
        i += 1
    f.close()
    # for i in range( n ):
        # yPred[i,:] = rand.permutation( shortList )[0:k]
    return yPred