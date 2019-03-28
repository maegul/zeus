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