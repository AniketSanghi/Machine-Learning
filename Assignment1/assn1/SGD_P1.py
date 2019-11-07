import numpy as np
import random as rnd
import time as tm
from matplotlib import pyplot as plt
import math

#  These are for 10 seconds
#  eta = 0.0007 , sq : 5400 , random
#  0.0008 , sq, 5491 , random
#  0.0008 , sq, 5415 , cyclic
#  0.0007 , sq, 5458 , cyclic
#  0.0009 , sq, 5384 , cyclic
#  0.00095 , sq, 5366 , cyclic
#  0.001 , sq, 5361 , cyclic
#  0.002 , sq, 5436 , cyclic
#  0.0015 , sq, 5307 , cyclic
#  0.0016 , sq, 5289 , cyclic
#  0.00165 , sq, 5283 , cyclic  //Not working now
#  0.00165 , sq, 5279 , random
#  Random is sometimes giving good answers but sometimes it is horrible
#  0.00003, -, 5252, cyclic
#  0.00005. -, 5234, cyclic
#  0.00007. -, 5224, cyclic
#  0.000072. -, 5224, cyclic


#  0.00165 , sq, 5250 , random cyclic
#  0.00007 , - , 5228 , random cyclic


def getCyclicCoord( currentCoord, d ):
	if currentCoord >= d-1:
		return 0
	else:
		return currentCoord + 1

def getRandCoord( d ):
	return rnd.randint( 0, d-1 )

randpermInner = -1

def getRandpermCoord( currentCoord, d ):
# 	# samples = rnd.sample( range(0, d), B )
# 	# return samples
	global randperm, randpermInner
	if randpermInner >= d-1 or randpermInner < 0 or currentCoord < 0:
		randpermInner = 0
		randperm = np.random.permutation( d )
		return randperm[randpermInner]
	else:
		randpermInner = randpermInner + 1
		return randperm[randpermInner]

def batch_grad(theta, C, X, y, j):
	(n, d) = X.shape
	X_ = X
	y_ = y
	discriminant = np.multiply((X_.dot(theta)), y_)
	g = np.zeros( (n,) )
	g[discriminant < 1] = -1
	return theta[j] + C * 2 * (g * (X_.T[j,:])).dot(np.multiply(y_, (1 - discriminant)))

def getObj( X, y, theta, C):
	w = theta[0:-1]
	hingeLoss = np.maximum( 1 - np.multiply( (X.dot( theta )), y ), 0 )
	return 0.5 * w.dot( w ) + C * hingeLoss.dot( hingeLoss )


def solver( X, y, C, timeout, spacing ):
	
	(n, d) = X.shape
	X = np.c_[X, np.ones(n)]
	d = d + 1
	t = 0
	totTime = 0
	totalTime = 0

	# w is the normal vector and b is the bias
	# These are the variables that will get returned once timeout happens
	# w = np.zeros( (d-1,) )
	# b = 0
	tic = tm.perf_counter()

	theta = np.zeros( (d,) )
	cumulative = theta
	w = theta[0:-1]
	b = theta[-1]
	eta = 0.00007
	# B = 5
	j = -1
	# global randpermInner randperm
	
	# randperm = np.random.permutation( d )

	obj_SGD = np.array([getObj(X, y, theta, C)])
	time_SGD = np.array([0])

	tic1 = tm.perf_counter()
	while True:
		t = t + 1
		if t % spacing == 0:
			toc = tm.perf_counter()
			totTime = totTime + (toc - tic)
			if totTime > timeout:

		
				print(getObj(X, y, theta, C), t)
				# print(t)
				return (w, b, totTime)
			else:
				tic = tm.perf_counter()
		

		thetanew = theta
		j = getCyclicCoord(j, d)
		# j = getRandCoord(d)
		# j = getRandpermCoord(j, d)

		thetanew[j] = theta[j] - (batch_grad(theta, C, X, y, j) * (eta)) 
		theta = thetanew

		toc1 = tm.perf_counter()

		# cumulative = cumulative + theta
		w = theta[0:-1]
		b = theta[-1]

	return (w, b, totTime) # This return statement will never be reached

# C = 1