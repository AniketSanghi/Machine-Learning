import numpy as np
import random as rnd
import time as tm
from matplotlib import pyplot as plt
import math

# 5400 eta = 0.000001, B = 5, eta/1, no. of Trials = 5
# 5433 eta = 0.000005, B = 5, eta/1, 
# 5600 eta = 0.0000005, B = 5, eta/1,
# 5400 eta = 0.0000009, B = 5, eta/1, 
# 5348 eta = 0.000002, B = 5, eta/1, 
# 5500 eta = 0.000003, B = 5, eta/1, 
# 5400 eta = 0.0000025, B = 5, eta/1, 
# 5292 eta = 0.000003, B = 10, eta/1, 
# 5327 eta = 0.000002, B = 10, eta/1, 
# 5403 eta = 0.000004, B = 10, eta/1, 
# 5327 eta = 0.0000025, B = 10, eta/1, 
# 5261 eta = 0.000003, B = 20, eta/1, 
# 5271 eta = 0.000004, B = 20, eta/1, 
# 5321 eta = 0.000002, B = 20, eta/1, 
# 5264 eta = 0.000003, B = 30, eta/1, 
# 5277 eta = 0.000004, B = 30, eta/1, 
# 5280 eta = 0.000002, B = 30, eta/1,
# 5280 eta = 0.000003, B = 40, eta/1, 
# 5248 eta = 0.000004, B = 40, eta/1, 
# 5326 eta = 0.000002, B = 40, eta/1,
# 5244 eta = 0.000003, B = 100, eta/1,
# 5241 eta = 0.000004, B = 50, eta/1,
# 5237 eta = 0.000003, B = 50, eta/1, 10s == 30s
# 5248 eta = 0.000004, B = 200, eta/1,
# 5480 eta = 0.000001, B = 500, eta/1,
# 5258 eta = 0.000003, B = 500, eta/1,


def getObj( X, y, w, b, C):
	hingeLoss = np.maximum( 1 - np.multiply( (X.dot( w ) + b), y ), 0 )
	return 0.5 * w.dot( w ) + C * hingeLoss.dot( hingeLoss )

def gradientSGD(X, y, C, w, B, n):

	if B <= n:
		random_indices = rnd.sample( range(0,n), B)
		Xnew = X[random_indices,:]
		ynew = y[random_indices]
	else:
		Xnew = X
		ynew = y

	discriminant = np.multiply( Xnew.dot(w) , ynew)
	g = np.zeros( (B,) )
	g[ discriminant < 1] = -1

	grad = w + C * (n/B) * 2 * (g * Xnew.T).dot(np.multiply(ynew, (1 - discriminant)))
	return grad


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

	eta = 0.0001
	B = 50
	
	obj_SGD = np.array(getObj(X[:,0:-1], y, w, b, C))
	time_SGD = np.array([0])
	total_time = 0
	tic1 = tm.perf_counter()

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
				return (w, b, totTime)
			else:
				tic = tm.perf_counter()

		
		delta = gradientSGD(X, y, C, w_, B, n)
		w_ = w_ - (eta/math.sqrt(t))*delta
		w = w_[0:-1]
		b = w_[-1]
		# print(getObj(X[:,0:-1], y, w, b, C))
		toc1 = tm.perf_counter()

		# total_time = total_time + (toc1 - tic1)
		# obj_SGD = np.append(obj_SGD, getObj(X[:,0:-1], y, w, b, C))
		# time_SGD = np.append(time_SGD, toc1-tic1)



		
		
	return (w, b, totTime) # This return statement will never be reached
