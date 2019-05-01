import numpy as np
from scipy.special import erf



def lenCurve_Exc(x, ke, a, c):
	'''
	Length Tuning Curve presuming single excitatory region

	Returns the integration of a gaussian from -length/2 to length/2 using ERF
	
	Takes an array of lengths and returns results of same shape

	Recommended bounds for fitting:
	([0,0,-np.inf],[np.inf,np.inf,np.inf])
	'''

	int_vals = np.array([-1, 1]) * x[:,np.newaxis] /2 / a

	vals = 0.5 * np.pi**0.5 * a * erf(int_vals)

	result = ke * np.diff(vals) + c
	
	return np.squeeze(result)



def lenCurve_Exc_Inh(x, ke, a, ki, b, oi, c):
	'''
	Length Tuning Curve presuming excitatory and inhibitory regions


	Returns the integration of a gaussian from -length/2 to length/2 using ERF
	
	Takes an array of lengths and returns results of same shape
	([0,0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
	'''

	int_vals = np.array([-1, 1]) * x[:,np.newaxis] /2
	
	exc_vals = 0.5 * np.pi**0.5 * a * erf(int_vals / a)
	exc_result = ke * np.diff(exc_vals)
	
	inh_vals = 0.5 * np.pi**0.5 * b * erf( (int_vals - oi) / b)
	inh_result = ki * np.diff(inh_vals)
	
	
	result = np.squeeze(exc_result - inh_result + c)
	result[result <0] = 0


	
	return result



def lenCurve_2Exc_Inh(x, ke, a, ke2, a2, oe, ki, b, oi, c):
	'''
	Length Tuning Curve presuming excitatory and inhibitory regions


	Returns the integration of a gaussian from -length/2 to length/2 using ERF
	
	Takes an array of lengths and returns results of same shape
	([0,0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
	'''

	# Parameter constraint logic
	if not (0 < oe < oi):
		return np.ones_like(x) * 1e6

	# if not (a2 < a):
	# 	return np.ones_like(x) * 1e6

	int_vals = np.array([-1, 1]) * x[:,np.newaxis] /2
	
	exc_vals = 0.5 * np.pi**0.5 * a * erf(int_vals / a)
	exc_result = ke * np.diff(exc_vals)

	exc_vals2 = 0.5 * np.pi**0.5 * a2 * erf( (int_vals - oe) / a2)
	exc_result2 = ke2 * np.diff(exc_vals2)
	


	inh_vals = 0.5 * np.pi**0.5 * b * erf( (int_vals - oi) / b)
	inh_result = ki * np.diff(inh_vals)
	
	
	result = np.squeeze(exc_result + exc_result2 - inh_result + c)
	result[result <0] = 0


	
	return result


def lenCurve_2Exc(x, ke, a, ke2, a2, oe, c):
	'''
	Length Tuning Curve presuming excitatory and inhibitory regions


	Returns the integration of a gaussian from -length/2 to length/2 using ERF
	
	Takes an array of lengths and returns results of same shape
	([0,0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
	'''

	# Parameter constraint logic

	# if not (a2 < a):
	# 	return np.ones_like(x) * 1e6

	int_vals = np.array([-1, 1]) * x[:,np.newaxis] /2
	
	exc_vals = 0.5 * np.pi**0.5 * a * erf(int_vals / a)
	exc_result = ke * np.diff(exc_vals)

	exc_vals2 = 0.5 * np.pi**0.5 * a2 * erf( (int_vals - oe) / a2)
	exc_result2 = ke2 * np.diff(exc_vals2)
	


	
	result = np.squeeze(exc_result + exc_result2  + c)
	result[result <0] = 0


	
	return result




def genMultiLin_3Seg(x, slopes, incpts, breaks):
	
	if (breaks[0] == x[0]) and (len(breaks) == (3+1)): # if start point in breaks too
		breaks = breaks[1:] # drop first break point
	
	assert len(breaks) == 3, 'breaks must have 3 numbers, one for end, and two for separation points'
	
	assert (breaks[0] < breaks[1] < breaks[2]), 'breaks not sequential'
	
	lines =  np.piecewise(
		x,
		[x<breaks[0], (x>=breaks[0]) & (x<breaks[1]), (x>=breaks[1]) & (x<=breaks[2])],
		[lambda x: x*slopes[0] + incpts[0], lambda x: x*slopes[1] + incpts[1], lambda x: x*slopes[2] + incpts[2]]
		 )
	
	return lines