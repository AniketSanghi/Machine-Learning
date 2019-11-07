import numpy as np
import random as rnd
import time as tm

# You may define any new functions, variables, classes here
# For example, functions to calculate next coordinate or step length

# Function to get next coordinate from a randomly generated permutation of [n]
def getRandPerm( currIndex, n ):
	
	global index, randperm
	
	if index < 0 or index > n-2 or currIndex == -1:
		index = 0
		randperm = np.random.permutation(n)
	else:
		index = index + 1

	return randperm[index]

################################
# Non Editable Region Starting #
################################
def solver( X, y, C, timeout, spacing ):
	(n, d) = X.shape
	t = 0
	totTime = 0
	
	# w is the normal vector and b is the bias
	# These are the variables that will get returned once timeout happens
	w = np.zeros( (d,) )
	b = 0
	tic = tm.perf_counter()
################################
#  Non Editable Region Ending  #
################################

	# You may reinitialize w, b to your liking here
	# You may also define new variables here e.g. eta, B etc

	# Hide the bias term within w_b
	# And add an extra column of 1 in X
	X_temp = np.ones((n,d+1))
	X_temp[:,0:-1] = X
	X = X_temp
	w_b = np.zeros((d+1,))

	# Initialising alpha with zeros
	# And storing the square norm of x before hand
	# Initialing w with its initial start value
	alpha = np.zeros( (n) )
	sq_norm_x = np.square( np.linalg.norm (X, axis = 1))
	w_b = X.T.dot(np.multiply(alpha, y))

	# Initialing our parameters to find index from random permutation everytime
	global randperm, index

	i = -1
	index = 0
	randperm = np.random.permutation(n)

################################
# Non Editable Region Starting #
################################
	while True:
		t = t + 1
		if t % spacing == 0:
			toc = tm.perf_counter()
			totTime = totTime + (toc - tic)
			if totTime > timeout:
				return (w, b, totTime)
			else:
				tic = tm.perf_counter()
################################
#  Non Editable Region Ending  #
################################

		# Write all code to perform your method updates here within the infinite while loop
		# The infinite loop will terminate once timeout is reached
		# Do not try to bypass the timer check e.g. by using continue
		# It is very easy for us to detect such bypasses - severe penalties await
		
		# Please note that once timeout is reached, the code will simply return w, b
		# Thus, if you wish to return the average model (as we did for GD), you need to
		# make sure that w, b store the averages at all times
		# One way to do so is to define two new "running" variables w_run and b_run
		# Make all GD updates to w_run and b_run e.g. w_run = w_run - step * delw
		# Then use a running average formula to update w and b
		# w = (w * (t-1) + w_run)/t
		# b = (b * (t-1) + b_run)/t
		# This way, w and b will always store the average and can be returned at any time
		# w, b play the role of the "cumulative" variable in the lecture notebook
		# w_run, b_run play the role of the "theta" variable in the lecture notebook

		# Get a coordinate from our random permutation
		i = getRandPerm(i, n)

		# Based on our mathematical derivation
		# as mentioned in pdf, calculating the 
		# values of variables in max O(d) :)
		x_i = X[i,:]
		q = sq_norm_x[i]
		p = y[i]*(x_i.dot(w_b) - alpha[i]*y[i]*q)

		# C multiplied on top instead of 1/2C in denominator
		# To avoid divide by zero exception
		alpha_i = 2*C*(1 - p)/(2*C*q + 1)

		# Taking the projection following the constraint
		if alpha_i < 0:
			alpha_i = 0

		# Updating the w for the new alpha_i and update alpha[i]
		w_b = w_b + (alpha_i - alpha[i])*y[i]*x_i
		alpha[i] = alpha_i

		# update w and b, in case the loop terminates after this
		w = w_b[0:-1]
		b = w_b[-1]
		
	return (w, b, totTime) # This return statement will never be reached