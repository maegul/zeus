from functools import partial

import numpy as np

from pyth.zeus import hermes, themis

def genPeakLen(x_curve, y_curve, thresh = 0.9, provide_arg = False):
	'''
	Returns values for when curve first gets to value that is
	equal to or greater than the max value times the threshold
	'''
	
	maxLen = np.max(y_curve)
	
	# Arg for first instance that passed condition
	# nonzero returns array for each dimension, if one, than a tuple with one element
	# Thus, two index extractions, one for typle, other for first arg of array
	peakLenArg = np.nonzero(y_curve > (thresh*maxLen))[0][0]

	peakLen = x_curve[peakLenArg]

	if provide_arg:
		return peakLenArg, peakLen
	else:
		return peakLen
	

def genGain(x_curve, y_curve):

	arg, peakLen = genPeakLen(x_curve, y_curve, provide_arg=True)
	
	peakResp = y_curve[arg]
	firstResp = y_curve[0]

	gain = (peakResp - firstResp) / peakResp

	return gain


def genESI(x_curve, y_curve):
	
	peakResp = np.max(y_curve)
	
#     arg, peakLen = genPeakLen(x_curve, y_curve)

#     peakResp = y_curve[arg]

	endResp = y_curve[-1]

	ESI = (peakResp-endResp) / peakResp

	return ESI


def analyseLengthData(proj, keys = None, cond_type = 'len', 
	n_interp = None, analysisFuncKeys = None):

	'''
	Returns dictionary of analysis data values for a proj object and
	set of analysis functions provided as a dict

	Parameters
	----
	keys : iterable, 1D (None)
		Run keys of the runs that are to be analysed.
		If None (default), then all run keys in proj will be analyse that
		match the cond_type specified.

	cond_type : str ('len')
		Cond type for filtering the run keys from proj if not keys are 
		specified.

	n_interp : int (None)
		number of points with which to generate the curve for a run key
		Passed to proj.getCurveData

	analysisFuncKeys : dict (None)
		Dictionary of functions to be used to analyse the data.
		Functions should be from this module, but need not be.

		Keys of this dict will be used for storing the result
		of the function in the output analysis data.

	'''

	if keys is None:
		# getting run_keys of the length tuning runs
		run_keys = proj.CellData.query('cond_type == @cond_type').index.get_level_values(1)

	else:
		assert hasattr(keys, '__iter__'), 'keys must be an iterable'

		run_keys = keys


	# Filter to those that have curves fit to them
	# run_keys = [rk for rk in run_keys if rk in proj.CurveFits]
	filtered_run_keys = []

	for rk in run_keys:
		if rk in proj.CurveFits:
			filtered_run_keys.append(rk)
		else:
			print(f'Run Key {rk} not in proj.CurveFits ... fit a curve?')

	if analysisFuncKeys is None:
		analysisFuncKeys = {
			'peak_length': genPeakLen,
			'ESI':genESI,
			'len_gain': genGain
		}

	analysisData = {}

	for rk in run_keys:

		analysisData[rk] = {}
		curve_data = proj.getCurveData(rk, n_interp=n_interp)

		for f in analysisFuncKeys.keys():
			analysisData[rk][f] = analysisFuncKeys[f](curve_data[0], curve_data[1])

	return analysisData






# ###
# ON OFF analysis
# ###

def get_on_off_resp(proj, run_key):
	'''
	Extracts the on and off responses from a themis object
	'''

	resp_data = proj.getData(key=run_key)

	# take max_responses out in order of conditions (just in case conditions are contrasts not 0,1 )
	off_max, on_max = resp_data[ 1, np.argsort(resp_data[0,:]) ]
	
	return on_max, off_max


def on_off_index(on_max, off_max):
	'''
	Returns the on off index
	Reflects whether overall response is dominated by ON or OFF
	or responses evenly to both

	1 -> ON cell
	0 -> OFF cell
	0.5 -> even

	Van Hooser et al (2013)
	'''
	
	return on_max / (on_max + off_max)


def sign_index(on_max, off_max):
	
	'''
	Returns sign index which reflects whether
	the cell responds only to one polarity or evenly
	to both

	0 -> responses evenly to both (not "signed")
	1 -> responds to one polarity ("signed")
	'''

	return np.abs(on_max-off_max) / (on_max + off_max)


def seg_index(dark_curve, light_curve):
	'''
	index of segregation to dark and light bars
	
	'''
	
	dl_diff = np.abs(dark_curve - light_curve)
	dl_sum = dark_curve + light_curve
	
	return dl_diff.sum() / dl_sum.sum()


def gen_on_off_index(proj, run_key):
	
	on_max, off_max = get_on_off_resp(proj, run_key)
	
	oo_index = on_off_index(on_max, off_max)
	
	return oo_index


def gen_sign_index(proj, run_key):
	
	on_max, off_max = get_on_off_resp(proj, run_key)
	
	s_index = sign_index(on_max, off_max)
	
	return s_index


def gen_seg_index(proj, run_key, resp='hist'):
	
	assert resp in ['hist', 'sdf'], f'resp "{resp}" not recognised'

	# Retrieving 
	cell_id = proj.getCellIdFromRunKey(run_key)
	themis_file_name = hermes.mk_themis_file_name(**cell_id)
	themis_obj = themis.load(themis_file_name)
	
	if resp == 'hist':
		resp_data = themis_obj.conditions_hist_mean
	elif resp == 'sdf':
		resp_data = themis_obj.spike_dens_func

	dark_idx = themis_obj.conditions.argmin()
	light_idx = themis_obj.conditions.argmax()

	if themis_obj.parameters['biphase_select_resp'] is not None:

		resp_sel = themis_obj.parameters['biphase_select_resp']

		half_mark = int(themis_obj.parameters['biphase_split_point'] * themis_obj.bins.size)
		dark_idx = themis_obj.conditions.argmin()
		light_idx = themis_obj.conditions.argmax()

		if resp_sel == 1:
			dark_curve = resp_data[dark_idx, :half_mark]
			light_curve = resp_data[light_idx, :half_mark]

		elif resp_sel == 2:
			dark_curve = resp_data[dark_idx, half_mark:]
			light_curve = resp_data[light_idx, half_mark:]
	else:

		dark_curve = themis_obj.conditions_hist_mean[dark_idx, :]
		light_curve = themis_obj.conditions_hist_mean[light_idx, :]

	si = seg_index(dark_curve, light_curve)
	
	return si


def analysePolarityData(proj, keys = None, cond_type = 'dl_bar',
	analysisFuncKeys = None, resp = 'hist'):

	'''
	Returns dictionary of analysis data values for a proj object and
	set of analysis functions provided as a dict

	Parameters
	----
	keys : iterable, 1D (None)
		Run keys of the runs that are to be analysed.
		If None (default), then all run keys in proj will be analyse that
		match the cond_type specified.

	cond_type : str ('len')
		Cond type for filtering the run keys from proj if not keys are 
		specified.

	analysisFuncKeys : dict (None)
		Dictionary of functions to be used to analyse the data.
		Functions should be from this module, but need not be.

		Keys of this dict will be used for storing the result
		of the function in the output analysis data.

	resp : 'str ('hist')
		Argument for gen_seg_index().
		Defines which response curve to use.
		Options: 'hist' (PSTH respone), 'sdf' (Spike density function)

	'''


	if keys is None:
		# getting run_keys of the length tuning runs
		run_keys = proj.CellData.query('cond_type == @cond_type').index.get_level_values(1)

	else:
		assert hasattr(keys, '__iter__'), 'keys must be an iterable'

		run_keys = keys


	if analysisFuncKeys is None:
		analysisFuncKeys = {
			'seg_idx': partial(gen_seg_index, resp=resp),
			'on_off_idx': gen_on_off_index,
			'sign_idx': gen_on_off_index
		}

	analysisData = {}

	for rk in run_keys:

		analysisData[rk] = {}

		for k, f in analysisFuncKeys.items():

			analysisData[rk][k] = f(proj, rk)

	return analysisData


