import numpy as np
import random as rnd
import time as tm
from matplotlib import pyplot as plt


def getObj( X, y, w, b, C):
	hingeLoss = np.maximum( 1 - np.multiply( (X.dot( w ) + b), y ), 0 )
	return 0.5 * w.dot( w ) + C * hingeLoss.dot( hingeLoss )

def getRandomPerm( curr , n):
	global randperm, randindex
	if randindex >= n-1 or randindex < 0 or curr < 0:
		randindex = 0
		randperm = np.random.permutation( n )
		return randperm[randindex]
	else:
		randindex = randindex + 1
		return randperm[randindex]

def getRandomCoord( curr, n):
	return rnd.randint(0, n-1)

def getCyclic(curr, n):
	if curr<0 or curr>=n-1:
		return 0
	else:
		return curr+1

def solver( X, y, C, timeout, spacing ):
	(n, d) = X.shape
	t = 0
	totTime = 0
	
	w = np.zeros( (d,) )
	b = 0
	tic = tm.perf_counter()

	# Hide the bias term
	Xtemp = np.ones( (n,d+1))
	Xtemp[:,0:-1] = X
	X = Xtemp
	w_ = np.zeros((d+1,))

	alpha = np.zeros( (n) )
	sq_norm_x = np.square( np.linalg.norm (X, axis = 1))
	w_ = X.T.dot(np.multiply(alpha, y))
	
	
	obj_SGD = np.array(getObj(X[:,0:-1], y, w, b, C))
	time_SGD = np.array([0])
	total_time = 0
	tic1 = tm.perf_counter()

	global randperm, randindex

	i = -1
	randindex = 0
	randperm = np.random.permutation(n)

	while True:
		t = t + 1
		if t % spacing == 0:
			toc = tm.perf_counter()
			totTime = totTime + (toc - tic)
			if totTime > timeout:
				# plt.plot( time_SGD, obj_SGD, color = 'r', linestyle = '-', label = "SGD" )
				# plt.xlabel( "Elapsed time (sec)" )
				# plt.ylabel( "C-SVM Objective value" )
				# plt.legend()
				# plt.ylim( 0, 20000 )
				# plt.xlim( 0, timeout )
				# plt.show()
				# print(getObj(X[:,0:-1], y, w, b, C), t)
				print(t)
				return (w, b, totTime)
			else:
				tic = tm.perf_counter()

		# i = getRandomCoord(i, n)
		i = getRandomPerm(i, n)
		# i = getCyclic(i, n)

		x = X[i,:]
		q = sq_norm_x[i]
		p = y[i]*(x.dot(w_) - alpha[i]*y[i]*q)
		
		alpha_new_i = 2*C*(1 - p)/(2*C*q + 1)

		if alpha_new_i < 0:
			alpha_new_i = 0

		w_ = w_ + (alpha_new_i - alpha[i])*y[i]*x
		w = w_[0:-1]
		b = w_[-1]
		alpha[i] = alpha_new_i

		# print(getObj(X[:,0:-1], y, w, b, C))
		toc1 = tm.perf_counter()

		# total_time = total_time + (toc1 - tic1)
		# obj_SGD = np.append(obj_SGD, getObj(X[:,0:-1], y, w, b, C))
		# time_SGD = np.append(time_SGD, toc1-tic1)



		
		
	return (w, b, totTime) # This return statement will never be reached
# 