import pathlib as pthl
# import inspect
# from pprint import pprint

import math
import numpy as np
import scipy as sp
import h5py # for importing matlab files formatted after version 7.3


import scipy.stats as stats
from scipy.ndimage import gaussian_filter1d
import numpy.fft as fft



# import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator
plt.style.use('ggplot')

import pandas as pd

import os
import glob

import pickle

from . import hermes
from . import hephaistos as heph
# from . import athena
from .utils import funcTracker, getMethods



def load_data(data_source, mode, cell_no=None):
	'''
	Wrapper around loading functions for all loading
	'''		

	if mode == 'txt':
		return load_raw_data(data_source, matlab=False)
	elif mode == 'matlab':
		return load_raw_data(data_source, matlab=True)
	elif mode == 'heph':
		return load_hephaistos(data_source, cell_no=cell_no)

def load_raw_data(dir, matlab=False):
	
	"""
	Loads files from the directory provided.
	
	Loads only txt files.
	
	Sorts and assigns expecting a number of either markerss, spikess times and spikess shape.
	
	Returns a dictionary of the imported data types
	"""

	current_dir = os.getcwd()    
	
	os.chdir(dir)
	
	file_names = []
	data = {}
	
	
	## For text files
	if not matlab:
		files = glob.glob('*.txt')
		
		assert len(files) > 0, 'No *.txt files found!'

		if len(glob.glob('*.mat')) > 0:
			print('WARNING: matlab files also found in directory: \t%s'%dir)
		
		for f in files:
			f_name = f.lower()
		
			if f_name.find('mark') > -1:
				data['markers'] = np.loadtxt(f_name, skiprows=1)
				file_names.append(f)
			
			elif f_name.find('spike') > -1:
				data['spikes'] = np.loadtxt(f_name, skiprows=1)
				file_names.append(f)
			
			elif f_name.find('shape') > -1:
				data['shape'] = np.loadtxt(f_name, skiprows=1)
				file_names.append(f)
	

	## For matlab files
	# These matlab files have more useful data than is extracted here.
	elif matlab:
		files = glob.glob('*.mat')
		
		assert len(files) > 0, 'No matlab files found!'
		
		if len(glob.glob('*.txt')) > 0:
			print('WARNING: text files also found in directory: \t%s' %dir)

		for f in files:
			f_name = f.lower()
			
			
			if f_name.find('mark') > -1:
				
				mark_file = h5py.File(f) # Loads hfd5 file
				mark_key = mark_file.keys()[0] # Gets name of relevant file for extract
				
				# Extract times of the markers
				data['markers'] = np.array(mark_file['%s/times' %mark_key])
				data['markers'] = np.reshape(data['markers'], -1) # turn to 1D array, as first axis redundant
				
				# Extract the numerical codes of the markers, which are listed one-to-one
				# with the times extracted above.  Useful for an integrity check.
				# Zero index necessary as marker codes has three empty columns
				data['marker_codes'] = np.array(mark_file['%s/codes' %mark_key][0])
				data['marker_codes'] = np.reshape(data['marker_codes'], -1) # turn to 1D array, as first axis redundant
				file_names.append(f)

			elif f_name.find('spike') > -1:

				spike_file = h5py.File(f) # Loads hfd5 file
				spike_key = spike_file.keys()[0] # Gets name of relevant file for extract
				
				# Extract times of the spikes
				data['spikes'] = np.array(spike_file['%s/times' %spike_key])
				data['spikes'] = np.reshape(data['spikes'], -1) # turn to 1D array, as first axis redundant


				#Extract trace for each spike. First Dim-trace, second-spikes.
				spike_traces = np.array(spike_file['%s/values' %spike_key])
				
				# Calculate Average shape (for all templates, which are coded in '/codes')
				avg_spike_trace = np.mean(spike_traces, axis=1)
				sem_avg_spike_trace = stats.sem(spike_traces, axis=1, ddof=1)
				
				data['shape'] = avg_spike_trace
				data['shape_SEM'] = sem_avg_spike_trace
				file_names.append(f)                
				
						
	os.chdir(current_dir)

			
	if len(data.keys()) != len(files):
		mesg = 'Not all of your file names are recognised; they may not have been imported appropriately.'
		mesg2 = 'File names must contain the key words "mark", "spike" and/or "shape."'
		print(mesg)
		print(mesg2)
		print('\nFollowing files loaded successfully:\n')
		for i in file_names: print(i)
		return data

	
	elif len(data.keys()) == len(files):
		print('All files imported and assigned')
		print('\nFollowing files loaded successfully:\n')
		for i in file_names: print(i)
		return data
		
		

def load_hephaistos(heph_unit, cell_no=None):

	# return heph native Spikes, Template, avg and SEM, markers and markerCodes

	# Extract 1d arrays for spike times, avg, SEM, marker times and marker Codes 



	data = {}

	data['markers'] = heph_unit.MarkTimes
	data['marker_codes'] = heph_unit.MarkCodes

	cluster_labels = heph_unit.Spikes.cluster_label.unique()

	if cell_no is None:

		assert cluster_labels.size == 1, \
			f'Hephaistos spike data has more than one cell ({cluster_labels}).  \nMUST PROVIDE Cell/cluster Number'

		spike_data = heph_unit.Spikes

	else:

		assert cell_no in cluster_labels, f'Cell no ({cell_no}) not in unit cluter labels ({cluster_labels})'

		spike_data = heph_unit.Spikes.query('cluster_label == @cell_no')


	data['spikes'] = spike_data.time.values

	data['shape'] =  heph_unit.SpikeAvgs.loc[cell_no, :].values
	data['shape_SEM'] = heph_unit.SpikeSem.loc[cell_no, :].values


	return data
		

def load(file_path):
	"""
	Un-Pickles a data object
	
	Parameters
	__________
	file_path : string
		File path of the pickled file to be unpickled
	"""
	assert type(file_path) == str, 'File_path must be a string'
	
	with open(file_path, 'rb') as f:
		return pickle.load(f)
		
		


#def plotform(ax, tickint=False):
#    
#
#
#    # remove top and right axes
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('left')
#    
#    # make ticks outward
#    ax.tick_params('both', direction='out')
#    
#    
##    fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
##    'weight' : 'normal', 'size' : 12}
#    
##    if tickint:
##        p.set_xticklabels(p.get_xticks().astype('int'), fontProperties)
##        p.set_yticklabels(p.get_yticks().astype('int'), fontProperties)
##    else:
##        p.set_xticklabels(p.get_xticks(), fontProperties)
##        p.set_yticklabels(p.get_yticks(), fontProperties)
#        
#        
#    #set the style of the major and minor grid lines, filled blocks
#    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
#    ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
#    ax.patch.set_facecolor('0.85')
#    ax.set_axisbelow(True)
#    
#    #set minor tick spacing to 1/2 of the major ticks
#    ax.xaxis.set_minor_locator(MultipleLocator( (plt.xticks()[0][1]-plt.xticks()[0][0]) / 2.0 ))
#    ax.yaxis.set_minor_locator(MultipleLocator( (plt.yticks()[0][1]-plt.yticks()[0][0]) / 2.0 ))
#    
#    #remove axis border
##    for child in ax.get_children():
##        if isinstance(child, matplotlib.spines.Spine):
##            child.set_alpha(0)
#       
#    #restyle the tick lines
#    for line in ax.get_xticklines() + ax.get_yticklines():
#        line.set_markersize(5)
#        line.set_color("gray")
#        line.set_markeredgewidth(1.4)
#    
#    #remove the minor tick lines    
#    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
#        line.set_markersize(0)
#
#    
#    
#    if ax.legend_ <> None:
#        lg = ax.legend_
#        lg.get_frame().set_linewidth(0)
#        lg.get_frame().set_alpha(0.5)
#        
#    plt.tight_layout()


def steppify(arr, axis='x'):
	"""
	expands datasets to create de facto step (step-post) plots
	"""
	
	if axis == 'x':
		newarr = np.r_[arr[0], np.dstack((arr[1:], arr[1:])).flatten()]
	
	elif axis == 'y':
		newarr = np.r_[np.dstack((arr[:-1], arr[:-1])).flatten(), arr[-1]]
	
	else:
		print('your axes in steppify are improperly identified')

	return newarr


#==============================================================================
# Generalised Analysis Functions
#==============================================================================

def fourier(data, temp_freq, axis, output = 'amplitude'):
	"""
	Generalised fourier analysis for drifting grating stimuli.
	
	
	Parameters
	__________
	
	axis : int
		Axis along which to perform fourier analysis
	
	"""
		
	
	# take largest possible multiple of F1 from PSTH.
	# Generate freq and fft
	# generate amplitude
	# return amplitude, F0, F1 and F2 values
	
	

def bootstrap(data, alpha=0.05, n_bootstrap = 2000, func=None, **func_args):
	
	"""Generalised bootstrapping function for data structured as con x trial x bins(time).
	
	Parameters
	__________
	data : array (3-D)
		Data to be bootstrapped.  Must be three dimensional with 
		trials along the second dimension (conditions on the first, time on the third)
	
	alpha : float
		Determines the confidence interval, where the probability of a type I error
		is `alpha`.  confidence intervals will be returned at the (100*alpha/2) and
		(100*(1-alpha/2)) percentiles.  Default is 0.05 (ie, 5%).
		
	n_bootstrap : int
		Number of iterations for the bootstrap distribution.  Default is 2000.
		Advisable that between 1000 and 10,000.
		
	func : function
		A function to operate on the bins or time dimension of the data before
		confidence intervals are calculated (philosophy being that the most
		non-normal data ought to be bootstrapped the most).
		
		Must accept 3-D data set, operate on the third dimension (axis=2),
		and return a 3-D data set.  If a pre-defined function is sought, wrapping
		in a custom function to meet these needs is advised.
		
	
	Returns
	_______
	out : tuple of arrays
		Tuple of 2 arrays, first the positive confidence interval, second, negative
		condidence interval, both at percentiles corresponding to a `alpha` chance
		of a type I error.  Each array is 2-D, with conditions on first axis, and bins
		(or time) on the second axis.
	"""
	
	assert data.ndim == 3, 'Data is not 3-dimensional.  Function only works for 3-D data.'    
	
	# Trials form the second dimension
	n_trials = data.shape[1]
	
	# generate randomised bootstrap resamples as random indices
	bootstrap_index = np.random.randint(0, n_trials, 
										(n_trials, n_bootstrap) )
	
	# For each bin in the histogram, randomly samples from the results
	# of each trial and repeats, effectively, n_bootstrap times                
	trials_bootstrap = data[:, bootstrap_index, :]
	
	# dimension one is the trials, zero is the conditions; this averaging 
	# goes across the trials creating a PSTH for each condition, and,
	# importantly, for each bootstrap resample
	avg_bootstrap = trials_bootstrap.mean(axis=1)
	
	if func:
		avg_bootstrap = func(avg_bootstrap, **func_args)
		
	# find percentile values for each bin along the bootstrap resamples,
	# which are on axis 1                                              
	CI_pos = np.percentile(avg_bootstrap, 100*(1 - (alpha/2.)), 
								axis=1)
	CI_neg = np.percentile(avg_bootstrap, 100*(alpha/2.), 
								axis=1)


	return CI_pos, CI_neg

	

## Class Definitions

class Themis:

	def __repr__(self):

		methods = getMethods(Themis, incl_hidden=True)

		# List of methods to report, in order of preferred/necessary execution
		method_report = [
			'__init__',		
			'_cell_ID',
			'_sort',
			'_conditions',
			'_analyse'


		]

		conditional_methods = {
			'__init__' : [
				'_assign_projects',
				'_stim_params',
				'_stim_params_print'
			],
			'_cell_ID' : '_cell_id_print',
			'_sort' : ['_marker_diag', '_plot_psth_flat'],
			'_conditions' : '_plot_psth',
			'_analyse' : '_plot_tuning'
		}

		rep = ''

		for m in method_report:

			if m in methods:

				rep += '{0}\t{1}()\n'.format(
					u"\u2705" if m in self._function_tracking else u"\u274c",
					m
					)

			if m in conditional_methods:
				condMethods = conditional_methods[m]
				if isinstance(condMethods, list):
					for cm in condMethods:
						rep += '\t\t{0}()\n'.format(cm)
				else:
					rep += '\t\t{0}()\n'.format(condMethods)


		return rep

	
	
	@funcTracker	
	def __init__(self, data_path = None, output_path = '.', 
		data_mode='matlab', cell=None):
		
		"""
		Instatiates the output of `zeus.load` as self.data.
		
		`zeus.load` argument, directory, is passed from class instantiation.
		
		For markers, spikes and shape, an attribute is directly .
		
		self.parameters instantiated, and directory argument appended.

		Parameters
		_________

		data_path : path / str / Hephaistos object
			Path to Hephaistos object or pickle, or
			Path to directory containing relevant data files 

		data_mode : str
			Relevant only if not using hephaistos objects or save files
			But instead using spike 2 matlab files or txt files directly

			'matlab' or 'txt', to define what kind of data to load
		
		"""

		# tracking function executions
		# Useful for repr and user guidance
		self._function_tracking = dict()

		# Load files and assign to attributes
		
		# Treat Hephaistos object differently 
		
		# Check if data_path is hephaistos object
		if data_path.__class__.__name__ == 'Hephaistos':
			# Take data directly from hephaistos object
			self.data = load_data(data_path, 'heph', cell_no=cell)

		# If not, treat as path / str
		else:
			self.Data_path = pthl.Path(data_path).absolute()
			assert self.Data_path.exists(), f'Does not exist; data_path: {self.Data_path}'

			# If path is to a pkl file, treat as a hephaistos save file
			if self.Data_path.suffix == '.pkl':
				unit = heph.load(self.Data_path)
				self.data = load_data(unit, 'heph', cell)

			# If a directory, then treat as directory of txt or matlab data files
			elif self.Data_path.is_dir():

				self.data = load_data(data_path, data_mode)



		self.Output_path = pthl.Path(output_path)

		self._absolute_paths = {}
		self._absolute_paths['Output_path'] = self.Output_path.absolute()

		if not self.Output_path.absolute().is_dir():
			self.Output_path.absolute().mkdir()




		
		try:
			self.markers = self.data['markers']
		except:
			pass

		try:
			self.marker_codes = self.data['marker_codes']
		except:
			pass


		try:
			self.spikes = self.data['spikes']
		except:
			pass

		try:
			self.shape = self.data['shape']
		except:
			pass
		
		try:
			self.shape_SEM = self.data['shape_SEM']
		except:
			pass


		self.parameters = {}
		


		# self.parameters['ops_directory'] = curr_dir

		# assert type(meta_data) == dict, 'meta data should be a ditionary from athena'
		
		
		print('\n\n*********\n')	
		print ('\nRun the following methods:'
				'\n self.sort() - to sort trials and generate PSTHs.'
				'\n self.conditions() - to generate and define condition descriptions and labels'
				'\nself._analyse() - to extract out the relevant measurement of response amplitude for each condition\n\n'
				'Saving and plotting functions may then be used.')

	def _assign_project(self, project):
		'''
		Accepts a path, or an Athena object to attach and assign a project
		to this themis object

		Parameters
		______
		project : Athena object | str (path)
			The athena object for the project, or a str providing the path to the 
			pickled object.
		'''

		if project.__class__.__name__ == 'Athena':
			# pull out project ID

			self.PROJ_ID = project.PROJ_ID
			# take path of
			self.PROJ_PATH = project.path
			self.ATHENA_PATH = project.SavePath

			self._PROJ_PATH_ABSOLUTE = project._absolute_paths['path']
			self._ATHENA_PATH_ABSOLUTE = project._absolute_paths['SavePath']

			project.save()


		# elif isinstance(project, str):
		else:
			try:

				athena_path = pthl.Path(project)

				with open(athena_path, 'rb') as f:
					project = pickle.load(f)

				self.PROJ_ID = project.PROJ_ID
				# take path of
				self.PROJ_PATH = project.path
				self.ATHENA_PATH = project.SavePath

				self._PROJ_PATH_ABSOLUTE = project._absolute_paths['path']
				self._ATHENA_PATH_ABSOLUTE = project._absolute_paths['SavePath']

				project.save()

			# Catching all exceptions :(
			except:
				print((
					f'Could not treat project as path or string\n'
					f'Type of project {type(project)}\n'
					f'No project assigned'
					))



		# Check if athena and themis are pointing to same directory
		# Using absolute of proj path, as paths may simply be '.', which is ok for
		# internal use, but, as athena and its path should already have been established
		# we use it as the reference and test for having the same output path is it

		if self.Output_path.absolute() != self._PROJ_PATH_ABSOLUTE:
			print((f'Themis output path ("{self.Output_path.absolute()}")'
				f'is not same as project path ("{self._PROJ_PATH_ABSOLUTE}")'
				))

		assert pthl.Path.cwd() == self.Output_path.absolute(), (
			f'Current path ("{pthl.Path.cwd()}") is not same as'
			f'Output path ("{self.Output_path.absolute()}")'
			'You should be working in the output path once coupling with a project!'
			)


	@funcTracker
	def _cell_ID(self, force = False, 
		experiment = None, unit = None, cell = None, run = None):


		cell_id = hermes.mk_cell_ID(
			experiment = experiment, unit=unit, cell=cell, run = run
			)

		self.CELL_KEY = f'{experiment}u{unit}c{cell}r{run}'


		if hasattr(self, 'PROJ_ID'):
			with open(self.ATHENA_PATH, 'rb') as f:
				athena_obj = pickle.load(f) 

			# Check if stim and cell dataset exists
			if hasattr(athena_obj, 'CellData'):

				# Check for repeats
				checkResult = hermes.checkCellData(cell_id, athena_obj.CellData, splitDataSet=True)

				if (force or checkResult):
					self.CELL_ID = cell_id

				else:
					print('potential redundancies, stim params NOT assigned')
					
			else:
				print('Athena has not cell data yet; stim params assigned')
				self.CELL_ID = cell_id
		else:
			print('No project assigned, stim params assigned')
			self.CELL_ID = cell_id


		# Assign cell ID
		self.CELL_ID = hermes.mk_cell_ID(
			experiment = experiment, unit=unit, cell=cell, run = run
			)



	def _stim_params(self, force=False, **kwargs):
		'''
		Update stimulus meta data in addition to the conditions of each test
		kwargs are a intended to be an exhaustive list of options
		'''
		# args = inspect.getargvalues(inspect.currentframe())

		stim_params = hermes.mk_stim_params(**kwargs)

		# Assigned projected?
		if hasattr(self, 'PROJ_ID'):
			with open(self.ATHENA_PATH, 'rb') as f:
				athena_obj = pickle.load(f) 

			# Check if stim and cell dataset exists
			if hasattr(athena_obj, 'CellData'):

				# Check for repeats
				checkResult = hermes.checkStimParams(stim_params, athena_obj.CellData, splitDataSet=True)

				if (force or checkResult):
					self.STIM_PARAMS = stim_params
					
			else:
				print('Athena has no cell data yet')
				self.STIM_PARAMS = stim_params
		else:
			print('No project assigned')
			self.STIM_PARAMS = stim_params





	def _stim_params_print(self):

		if hasattr(self, 'STIM_PARAMS'):

			print('Stimulus Parameters\n')

			hermes.show_info(self.STIM_PARAMS)

		else:
			print('Stim params not set yet ... use self._stim_params')


	def _cell_id_print(self):

		if hasattr(self, 'CELL_ID'):

			print('Cell and Experiment ID:\n')

			hermes.show_info(self.CELL_ID)

		else:
			print('No cell ID!!')



	@funcTracker	
	def _sort(self, use_codes = False, conditions = 9, trials = 10, stim_len = None,
			  bin_width=0.02, auto_spont=False, sigma=3, alpha=0.05, 
			  n_bootstrap=2000):
		
		"""Sorts spikes into conditions and trials on the basis of the markers
		loaded by `zeus.__init__` and the expected amount of conditions and
		trials defined by the passed arguments.  Then, generates histograms, 
		or PSTHS, for each of the conditions, averaged across the trials, 
		as well as filtered spike density functions, both with confidence intervals.
		
		Parameters
		__________
		use_codes : boolean
			Instead of identifyin amount of conditions and trials, use the marker 
			codes to derive these numbers.
			Must have marker codes from the data source
			If True, "conditions" and "trials" are ignored

		conditions : int
			Amount of stimulus conditions in the experiment.
			
		trials : int
			Amount of trials of each condition in the experiment.
			Each condition must have the same amount of trials.
			
		stim_len : float, optional
			Default is `None`.
			Defines the length of the stimulus for a single trial.
			If not defined, it is calculated as the maximum difference in time
			between markers.
			
		bin_width : float, default = 0.02
			Determines the width of the bins in histograms or PSTHs
			
		auto_spont : Boolean
			Not yet complete.  Determines whether to calculate a spontaneous rate automatically.
			Currently, only sure method for this is to count spikes before the first stimulus.
			More methods are likelty to be used later.
			
		sigma : float
			Defines the standard deviation, in bins, of the gaussian kernel applied
			to time-series data to smoothen the data.  To obtain the standard 
			deviation in units of absolute time -> `self.bin_width * sigma`.
			
		alpha : float
			To be passed to global method `bootstrap()`
			Determines the confidence interval, where the probability of a type I error
			is `alpha`.  confidence intervals will be returned at the (100*alpha/2) and
			(100*(1-alpha/2)) percentiles.  Default is 0.05 (ie, 5%).
		
		n_bootstrap : int
			To be passed to global method `bootstrap()`
			Number of iterations for the bootstrap distribution.  Default is 2000.
			Advisable that between 1000 and 10,000.

		
		Returns
		_______
		
		self.conditions_trials : dict of lists
			Three dimensional array containing the REAL spike times of 
			each trial of each condition.
			Keys are the codes of the actual conditions
			
		self.conditions_trials_zeroed : dict of lists
			Three dimensional array containing the spike times, relative to the
			commencement of the relevant trial, of each trial of each 
			condition.
			Keys are the codes of the actual conditions
			
		self.parameters : dict, pre-instantiated
			'conditions' (Number of conditions)
			'trials' (Number)

		self.bin_width : float
			Width of the bins used in generating the histograms.  Defines in seconds.
			
		self.bins : array (1-D)
			The bind boundaries used in generating the histograms
			
		self.conditions_trials_hist : array (3-D)
			PSTH/histograms for each trial for each condition.  
			Provides spike counts for each trial
			First dimension (axis=0) - conditions (numerical: range(conditions.size))
			Second dimension (axis=1) - trials
			Third dimension (axis=2) - bins (or time)
			
		self.conditions_hist_mean : array (2-D)
			As `self.conditions_trials_hist` but averaged across trials, 
			provides the mean spike count.
			Second dimension (axis=1) - bins (or time),
			
		self.conditions_hist_stderr : array (2-D)
			As `self.conditions_hist_mean` but with standard error of the mean rather than the mean
			itself provided.
			
		self.conditions_trials_sdf : array (3-D)
			As `self.conditions_trials_hist` but with gaussian smoothing applied
			to each trial. 
			This gaussian is 5-sigmas wide on each side (10 total) and has 
			`sigma` standard deviation. It is the same gaussian applied for 
			the spike density function.
			
		self.spike_density_function : array (2-D)
			As `self.conitions_hist_mean` but with gaussian smoothing applied to
			the mean values.
			
		self.CI_pos : array (2-D)
			As `self.conditions_hist_stderr` but with bootstrap derived confidence
			intervals, for the positive side.  
			Accompanies `self.spike_density_function`, as the intervals are derived
			from a population sampled from smoothened values.  
			Generated by global method `bootstrap()`.

		"""
		# Loading parameters to parameters dictionary

		if use_codes:
			assert self.marker_codes.size, 'No marker codes from data source, cannot use codes'

			conds, cond_counts = np.unique(self.marker_codes, return_counts=True)

			assert np.all(cond_counts == cond_counts[0]), (
				'Not all conditions have the same amount of trials according to marker codes'
				)

			conditions = conds.size
			trials = cond_counts[0]

		else:
			# useful for iterating later in the function
			# start at 1, +1 to get to last condition
			conds = np.arange(1, conditions+1)


		self.parameters['conditions'] = conditions
		self.parameters['trials'] = trials

		self.parameters['sdf_alpha'] = alpha
		self.parameters['sdf_number_bootstrap'] = n_bootstrap     
		self.parameters['sdf_sigma'] = sigma       

		self.parameters['bin_width'] = bin_width
		
#==============================================================================
# Post-marker buffer & Spontaneous calculation needed.
#==============================================================================
		
		assert self.markers.size == (conditions * trials), (
			'The number of recorded markers' 
			'does not match your dataset specifications.'
			f'\n#Conditions ("{conditions}") * #trials ("{trials}") != #markers ("{self.markers.size}").\n'
			'Markers will probably have to be inserted.  Run self._marker_diag() to diagnose')

		# If marker codes generated by data source
		if self.marker_codes.size:
			assert self.markers.size == self.marker_codes.size, (
				f'Number of markers ("{self.markers.size}") != number of marker codes ("{self.marker_codes.size}") '	
				)
		
		
		# determine stimulus length
		if stim_len:
			self.parameters['stimulus_length'] = float(stim_len)
		else:
			stim_diffs = np.diff(self.markers)
			stim_len = np.max(stim_diffs)
			self.parameters['stimulus_length'] = stim_len
			
			print('\nstimulus length has been calculated from the set of markers; ' 
				  'taking the maximum distance between all adjacent pairs of markers')
			print(f'Calculated Stim Len (max): {stim_len}\nMean Stim Len: {np.mean(stim_diffs)}\nMedian Stim Len: {np.median(stim_diffs)}')
			
		

		# If no markcodes, make synthetic marker codes
		if not self.marker_codes.size:

			self.marker_codes = (

				# reshape to enable broadcasting for repeats
				conds.reshape(1, conditions)
				* np.ones((trials, conditions)) # trials dimension will have the conditions from the arange broadcast along it (ie, repeated)
				
				).flatten().astype('uint8')


		# Initialising 

		self.conditions_trials = {}
		self.conditions_trials_zeroed = {}


		for c in conds:
			# these dictionaries are now keyed by the actual codes for the conditions
			self.conditions_trials[c] = []
			self.conditions_trials_zeroed[c] = []

		cond_range = range(conds.size)
		# as conds are the actual codes used, this dictionary converts
		# from numerical index (from zero upward) to actual code used
		self.condition_code_idx_convert = {cond_range[i]:conds[i] for i in range(conds.size)}

		# Iterate through each marker
		for i, t in enumerate(self.markers):

			# extract spike times for each trial
			trial = np.extract( (self.spikes >= t) & (self.spikes < t+stim_len), 
								self.spikes )
			# Also zero to beginning of trial
			trial_zeroed = trial - t

			# condition identify from marker codes
			cond = self.marker_codes[i]

			# Append both trials and zeroed trials to conditions dictionary
			self.conditions_trials[cond].append(trial)
			self.conditions_trials_zeroed[cond].append(trial_zeroed)
			
#==============================================================================
# Generating PSTHs
#==============================================================================
	
		self.bin_width = bin_width
		
		# bin boundaries, in time.  Fundamental to basic PSTH analysis and so made an attribute
		self.bins = np.arange(0, stim_len, bin_width)
		
		
		# Initialising.  Will contain, for each condition, a PSTH for each trial
		self.conditions_trials_hist = np.zeros((conditions, trials, self.bins.size - 1))
		
		# One condition, and trial, at a time
		for ic,c in enumerate(conds):
			
			trials_hist = np.zeros((trials, self.bins.size - 1))
			
			for t in range(trials):
				
				trials_hist[t] = np.histogram(self.conditions_trials_zeroed[c][t], self.bins)[0]
				
			# PSTH for each trial added for the condition    
			self.conditions_trials_hist[ic] = trials_hist
			
		# Average and Standard Error, for each bin, calculated across the trials (second axis)    
		self.conditions_hist_mean = np.mean(self.conditions_trials_hist, axis=1)
		self.conditions_hist_stderr = stats.sem(self.conditions_trials_hist, axis= 1, ddof=0)


#==============================================================================
# Create more robust spontaneous function
#==============================================================================
		# Calculate spontaneous rate from unstimulated portion of the record
		if auto_spont:
			index = self.spikes < self.markers[0]
			self.spont_rate = self.spikes[index].size / self.markers[0]



#==============================================================================
# Spike Density Function generation - previously _sdf()            
#==============================================================================
	   
		def gauss_smooth(data, sigma):
			"""
			Wrapper for guassian smooting, to be passed to bootstrap function
			"""
			# make the kernel 5 sigmas wide in each direction
			kernel = stats.norm.pdf(np.arange(-5*sigma, (5*sigma)+1), scale=sigma)
			
			return sp.ndimage.convolve1d(data, kernel, axis=2)


		# Gaussian Filter
		kernel = stats.norm.pdf(np.arange(-5*sigma, (5*sigma)+1), scale=sigma)

		# Confidence intervals derived from bootstrapping                                      
		self.CI_pos, self.CI_neg = bootstrap(self.conditions_trials_hist, alpha=alpha,
											 func = gauss_smooth, **{'sigma':sigma})
							   
		# Gaussian smoothing for each trial
		self.conditions_trials_sdf = sp.ndimage.convolve1d(self.conditions_trials_hist,
															kernel, axis=2)
		# Gaussian smoothing over condition averages, for a spike density function
		self.spike_dens_func = sp.ndimage.convolve1d(self.conditions_hist_mean, 
													 kernel, axis=1)
			
			
	def _marker_diag(self, auto_replace = False):
		"""Attempts to pinpoint where markers have been lost or missed in the marker stream.
		
		Parameters
		__________
		auto_replace : Boolean
			Not yet active.  Intention is to allow auto-replacenemtn of missing markers
			
		Returns
		_______
		table : dict
			Dictionary containing table of where markers are missing and key numbers regarding
			how many markers are missing.
		"""
		

		try:
			self.parameters['conditions']
		except:
			print('THere are no sorting parameters ... you must run self._sort()')
			return

		num_missing_markers = (self.parameters['conditions']*self.parameters['trials']) - self.markers.size        
		
		if self.parameters['conditions']*self.parameters['trials'] == (self.markers.size):
			print(
				'Stimulus parameters appear to be appropriate.  Running this function is '
				'unlikely to be necessary.  If anomalies are found, timing errors during the '
				'experimen are likely to be the culprit, and perhaps irretrivable.')
		
		output = {}

		marker_diff = np.r_[0.0, ( np.subtract( self.markers[1:], self.markers[:-1]))]
		
		multiple = np.around((marker_diff / float(np.median(marker_diff))) - 1, 1)
	
		anomaly = np.greater(marker_diff, np.ones_like(marker_diff) * 1.2 * np.median(marker_diff))
		
		table = pd.DataFrame({'markers': self.markers, 'multiple': multiple, 
							  'difference': marker_diff, 'anomalous': anomaly})

		
		# Anomalous markers, require insertion
		
		anomalies = table[table['anomalous'] == True]
		output['bad_marks'] = anomalies
		output['bad_marks_index'] = anomalies.index
		output['num_mark_found'] = np.around(anomalies['multiple'].sum(), 1)
		output['num_missing_mark'] = num_missing_markers
		output['num_missing_mark_at_beg/end'] = max(output['num_missing_mark']-output['num_mark_found'], 0)
		output['maximum_marer_diff'] = marker_diff.max()
		output['median_marker_diff'] = np.median(marker_diff)

		if output['num_missing_mark_at_beg/end'] > 0:
			output['spikes_before_beg'] = self.spikes[self.spikes < self.markers[0]]

			output['spikes_after_end'] = self.spikes[self.spikes > self.markers[-1]]
			
#==============================================================================
# Add:
#            
#            Auto Replace function.  Create a parameter reflecting, and do not touch the original
#            file.
#            Capacity to tell whether marker missing from front or end of experiment.
#            Use average spike rates to test, or require feedback from teh original data file.        
#==============================================================================
			
		return output
			
			 
	@funcTracker
	def _conditions(self, beg=-90, intvl=20, con_type='ori', stim='bar', 
					biphasic=True, unit='deg', con_list=[], temp_freq = 2):
						
		"""Generates condition and stimulus descriptions both as numbers and strings.
		
		Parameters
		__________
		
		beg : float
			First condition parameter/description.  Presumes a linear series for
			list of conditions.  If not true, use `list` parameter to provide 
			conditions list manually
			
		intvl : float
			Increment by which condition series increases from `beg`.  Increment 
			will be iterated as many times as there are conditions (derived from 
			`self.parameters`)
			
		con_type : string
			Condition type.  Stimulus parameter varying in experiment.  
			Predefined defaults are - `orientation` - `spat_freq` - `temporal freq` - `chromatic`.  
			Default is `orientation`.
			If a predifined type is not provided, it is presumed to describe a
			linear range of conditions according to `beg` and `intvl`.
			
		stim : string
			Stimulus type.  
			Select either - `bar` - `grating`.
			Default is `bar`.
			
		biphasic : Boolean
			Whether bar is biphasic, that is, whether it reverses direction in the
			one trial.
			
		unit : string
			Units of the condition parameters.  Default is `deg`.
			
		con_list : list (of strings)
			List of condition descriptions.  Must be strings convertible to floats,
			and provide one for each condition.
			
		temp_freq : float
			temporal frequency of drifting grating.  Required if `stim = 'grating'`.
			
		
		Returns
		_______
		
		self.conditions : array (1-D)
			Condition descriptions in numbers
			
		self.conditions2 : array (1-D)
			For when `biphasic == True`.  Lists the secondary condition descriptions for the
			returning bar.
			
		self.cond_label : list (of strings)
			String descriptions of conditions, chiefly for plotting purposes.
		"""
		
		
		con_types = ['ori', 'spat_freq', 'temporal_freq', 'chromatic']
		stims = ['bar', 'grating']
		
		
		# Checking if condition and stimulus type recognised.  
		if not con_type.lower() in con_types:
			print('con_type not recognised. '  
					'Predefined options, if desired, are %s \n'%con_types
					)

		if not stim.lower() in stims:
			print('stimulus not recognised. '  
					'Predefined options, if desired, are %s \n'%con_types
					)


		
		n_con = self.parameters['conditions']
		
		self.parameters['condition_type'] = con_type.lower()
		self.parameters['condition_unit'] = unit.capitalize()
		self.parameters['stimulus'] = stim.lower()
		
		if stim.lower() == stims[1]:
			# Gratings are GENERALLY not biphasic
			self.parameters['biphasic'] = 'N/A'
		else:
			self.parameters['biphasic'] = biphasic
		
		# Address issue of whether the sampling rate suits teh temporal frequency of 
		# the grating for FFT analysis
		if stim.lower() == 'grating':
			self.parameters['temp_freq'] = float(temp_freq)
			
			# Sample rate must be a multiple of F1/temp_freq for it to be a frequency measured
			# in the FFT.
			samp_rate = 1/float(self.bin_width)
			
			
			assert samp_rate % temp_freq == 0., ('Bin_width (%s) is incompatible wih obtaining' 
												 'an FFT containing the specified temp_freq (%s). '
												 'The sampling frequency (1/bin_width) must be a'
												 'multiple of the temp_freq. \n\n Try as a' 
												 'bin_width %s and rerun self._sort().'
												 % (self.bin_width, temp_freq, 
													1/(np.ceil(samp_rate/float(temp_freq))*temp_freq)))
		
		self.cond_label = []

		
		def circ(ori, bound = 360):
			"""Func that Ensures all orientation values are between 0 and 360 degrees.
			"""
			# ori[ori<-360] += 720
			# ori[ori<0] += 360
			# ori[ori>360] -= 360
			# ori[ori>720] -= 720


			return ori % bound

		# if list of conditions provided directly
		if len(con_list) > 0:
			
			# Must match number of conditions
			assert len(con_list) == n_con, ('the number of labels provided '
										'manually (%s) does not match the '
										'number of conditions (%s).' % 
										(len(con_list), n_con))
			 
			# Must all be strings                           
			assert all(isinstance(l, str) for l in con_list), ('not all the '
														   'labels provided '
														   'are strings')
																										 
			# List of conditions as strings
			self.cond_label = con_list
			
			# Convert to floats
			# Relying on numpy conversion error should list be unable to convert to float.
			self.conditions = np.array(con_list).astype('float')
			
			
			if biphasic:
				

				# self.conditions has been defined as an np.ndarray
				self.conditions2 = self.conditions 
# 
#                 # Generate list of strings or labels
#                 for c in range(n_con):
#                     label = '%s / %s %s' %(self.conditions[c], self.conditions2[c],
#                                            self.parameters['condition_unit'])
#                     self.cond_label.append(label)
# 
#             else:
#                 for c in range(n_con):
#                     
#                     label = '%s %s' %(self.conditions[c],
#                                       self.parameters['condition_unit'])
#                     self.cond_label.append(label)

				
		
		# if condition tpye is orientation
		elif con_type.lower() == con_types[0]:
			
			# Generate full range of conditions
			self.conditions = circ(np.arange(beg, beg+(n_con*intvl), intvl))
			
			assert len(self.conditions) == n_con, ('The amount of condition labels (%s) '
											'and conditions (%s) do not match; '
											'check your condition parameters' % 
											(self.cond_label.size, n_con))
			
			if biphasic:
				

				# self.conditions has been defined as an np.ndarray
				self.conditions2 = circ(self.conditions + 180) 

				# Generate list of strings or labels
				for c in range(n_con):
					label = '%s / %s %s' %(self.conditions[c], self.conditions2[c],
										   self.parameters['condition_unit'])
					self.cond_label.append(label)
			# Generate list of strings for non-biphasic.        
			else:
				
				for c in range(n_con):
					label = '%s %s' %(self.conditions[c],
									  self.parameters['condition_unit'])
					self.cond_label.append(label)
					
		# IF condition type is Spat Freq            
		elif con_type.lower() == con_types[1]:
			self.conditions = np.arange(beg, beg + (n_con*intvl), intvl)
			
			assert len(self.conditions) == n_con, ('The amount of condition labels (%s) '
											'and conditions (%s) do not match; '
											'check your condition parameters' % 
											(self.cond_label.size, n_con))

			for c in range(n_con):
				label = '%s %s' %(self.conditions[c], self.parameters['condition_unit'])
				self.cond_label.append(label)

					
		# if condition type is not predefined in this method, presume linear range           
		elif not con_type.lower() in con_types:
			
			self.conditions = np.arange(beg, beg+(n_con*intvl), intvl)


			if biphasic:
				

				# self.conditions has been defined as an np.ndarray
				self.conditions2 = self.conditions 

				# Generate list of strings or labels
				for c in range(n_con):
					label = '%s / %s %s' %(self.conditions[c], self.conditions2[c],
										   self.parameters['condition_unit'])
					self.cond_label.append(label)

			else:
				for c in range(n_con):
					
					label = '%s %s' %(self.conditions[c],
									  self.parameters['condition_unit'])
					self.cond_label.append(label)





	@funcTracker
	def _analyse(self, source='sdf', alpha = 0.05, n_bootstrap = 2000, 
		biphas_split_point = 0.5, biphase_select_resp = None):
		
		
		# Need to add capacity to handle two things:
		# Qualitative conditions ne
		# Conditions split - need to to deal with splitting a single dataset, where one part
		# is qualitiative and the other quantitative
		
		# For qualitative, the self.cond_tuning array is numerical.  Replace with record
		# see Initial Chrom Analysis.  Keep conditions as strings, and convert to numerical
		# for plotting (?).  Where qualitative, only use bar plot, where mixed, split.
					
		
		"""
		Extracts relevant response amplitude for each condition.

		Relies on previously defined stimulus type in deciding how to extract response
		amplitude.  
		
		If it is a bar, the maximum response is taken.  If the bar stimulus was biphasic,
		the stimulus is treated as being two adjacent conditions accordingly.
		
		If the stimulus is a grating, an FFT is taken along with bootstrapped confidence intervals.
		
		
		
		Provides arrays and Pandas Dataframes of compiled response amplitudes with confidence
		intervals.
		
		Parameters
		__________
		
		source : str
			The spike data to be used for ensuing analysis.
			Select one of - `sdf`, `mov_avg` or `raw`.
		
		alpha : float
			For bootstrap confidence intervals generated for the fourier analysis.
			Defines the confidence interval as 100*(1-alpha).
		
		n_bootstrap : int
			For bootstrap confidence intervals generated for the fourier analysis.
			Defines how many trials are re-sampled.

		biphas_split_point : float
			If biphasic, what point to split the PSTH.
			0.5 -> half, 0.25 -> at quarter point, etc.

		biphase_select_resp : int
			If stimuli is biphasic, but desire only one set of responses
			Biphase ids start are 1 and 2 (for first and second response)
			Provide integer for which response is to be recorded, exclusively,
			in cond_tuning and cond_tuning_pd

			
		
		Returns
		_______
		
		self.cond_tuning : array(2-D)
			Response to each condition with positive and negativeconfidence intervals.
			First row (axis 0 - [0,:]) - conditions
			Second row - response amplitude
			Third row - negative confidence interval, the 95% lower limit of the amplitude
			Fourth row - postive confidence interval
		
		self.cond_tuning_pd : Pandas Dataframe
			Pandas Dataframe of self.cond_tuning.
			Transposed such that conditions are listed in the first column rather than row, and
			so on.
			
		
		
		"""
		## Add parameters to parameters dictionary
		
		# Organising source selection - raw and mov_avg not develoepd fully yet.
		sources = {'sdf': (self.spike_dens_func, self.CI_pos, self.CI_neg), 
				   'mov_avg': 'doesnt exist yet, call it self.spike_mov_avg', 
				   'raw': (self.conditions_hist_mean, 
						   self.conditions_hist_mean + 2*self.conditions_hist_stderr, 
						   self.conditions_hist_mean - 2*self.conditions_hist_stderr)}
				   
		assert source.lower() in sources.keys(), ('Tuning source data "%s" is invalid '
													'select one of %s' %(source, sources.keys()))      
	   
		## Need to expand this functionality to the mean and CI_pos and CI_neg.  Doing so for
	   # raw and moving average is not a priority, using sdf and bootstrap is pretty good.
	   # overall aim is to clean this function up to accomadte a number of tuning functinos
	   # in a clear and easy to use fasion.
	   
		n_con = self.parameters['conditions']
		
		# values for transient bar responses
		if self.parameters['stimulus'] == 'bar':
			
			resp, CI_pos, CI_neg = sources[source.lower()]
			
			
			if self.parameters['biphasic']:
				
				# Take max response for each half of each PSTH, including Conf Intvls
				half = int(self.bins.size * biphas_split_point)

				max_val_arg = (resp[:, :half].argmax(axis=1),
							   resp[:, half:].argmax(axis=1)+half)
									
				max_val = (resp[:, :half].max(axis=1),
						   resp[:, half:].max(axis=1))
						   
							   
				max_val_CI_neg = (CI_neg[np.arange(n_con), max_val_arg[0]],
								  CI_neg[np.arange(n_con), max_val_arg[1]])
								  
				max_val_CI_pos = (CI_pos[np.arange(n_con), max_val_arg[0]],
								  CI_pos[np.arange(n_con), max_val_arg[1]])

				# encode which of the two responses the data is attached to
				biphas_id = np.zeros_like(np.hstack((self.conditions, 
														self.conditions2)))
				biphas_id[:self.conditions.size] = 1
				biphas_id[self.conditions2.size:] = 2


								  
				self.cond_tuning = np.vstack((np.hstack((self.conditions, 
														self.conditions2)),
											 np.hstack(max_val),
											 np.hstack(max_val_CI_neg),
											 np.hstack(max_val_CI_pos),
											 biphas_id))
											 
				# Convert to Hertz - design choice is to keep all PSTH datasets as raw average spike
				# counts, with easy option of seeing frequency in the plotting, but converting to 
				# Hertz for all condition tuning data.
				self.cond_tuning[1:-1,:] *= (1/self.bin_width)
											 
							  
				# Column labels for pd.dataframe of tuning data
				# Percentage of confidence intervals
				# ci_perc = (100 * (1 - self.parameters['sdf_alpha']))
				
				# Labels
				idx = ['condition', 'max_resp', 'neg_CI', 'pos_CI', 'biphas_id']
				
				# Pandas object, with transpose of tuning array to data frame object       
				self.cond_tuning_pd = pd.DataFrame(self.cond_tuning.transpose(), columns=idx)
		
			#non biphasic version of above
			if not self.parameters['biphasic']:

				max_val_arg = resp[:, :].argmax(axis=1)
									
				max_val = resp[:, :].max(axis=1)
						   
							   
				max_val_CI_neg = CI_neg[np.arange(n_con), max_val_arg]
								  
				max_val_CI_pos = CI_pos[np.arange(n_con), max_val_arg]
								  
				self.cond_tuning = np.vstack((self.conditions,
											 max_val,
											 max_val_CI_neg,
											 max_val_CI_pos))
											 
				# Convert to Hertz - design choice is to keep all PSTH datasets as raw average spike
				# counts, with easy option of seeing frequency in the plotting, but converting to 
				# Hertz for all condition tuning data.
				self.cond_tuning[1:,:] *= (1/self.bin_width)
				
				
				# Column labels for pd.dataframe of tuning data
				# ci_perc = (100 * (1 - self.parameters['sdf_alpha']))

				idx = ['condition', 'max_resp', 'neg_CI', 'pos_CI']

					   
				# transpose of tuning array to data frame object       
				self.cond_tuning_pd = pd.DataFrame(self.cond_tuning.transpose(), columns=idx)
		
		
		# values for sinusoids/gratings
		## Note issue of temporal frequency tuning - need variable tf.
		if self.parameters['stimulus'] == 'grating':
			
			self.parameters['fft_alpha'] = alpha
			self.parameters['fft_number_bootstrap'] = n_bootstrap
			
			if source == 'sdf':
				print ('WARNING, using a smoothed/filtered dataset will artificially increase'
					   'the amplitude of the DC component and decrease that of the F1') 
			
			sources = {'sdf': self.conditions_trials_sdf,
					   'mov_avg': "doesn't exist yet (?)",
					   'raw': self.conditions_trials_hist}
			
			resp = sources[source]
			
			temp_freq = self.parameters['temp_freq']
			stim_len = self.parameters['stimulus_length']
			
			# ensuring that the temp_freq is measured in the FFT whilst taking the maximum time.
			# on the basis of delt-f = 1 / n*del-t; stim_len*F1=factor; 1/(bin_width*F1)=min bins
			# number times greater than minimum can fit in stim_length            
			factor = np.floor(stim_len * temp_freq).astype('int')
			
			# number of bins to take - the window size necessary for temp_freq to be measured
			bins_take = np.floor(factor / (self.bin_width * temp_freq)).astype('int')

			# Frequency axis generation
			self.freq = fft.rfftfreq(bins_take, self.bin_width)
			
			#Checkign whether the temp_freq is in the FFT.
			assert self.freq[factor] == temp_freq, ('The calculated FFT F1 frequency (%s)'
													   'does not equal the Stimulus temp_freq (%s)'
													   %(self.freq[bins_take], temp_freq))

			# Fourier Transform
			self.conditions_trials_fourier = fft.rfft(resp[:,:,:bins_take], axis=2)
			
			# Amplitude (peak-to-peak)
			self.conditions_trials_ampl = np.abs(self.conditions_trials_fourier)
			
			# normalising to dataset size, except the DC.
			self.conditions_trials_ampl[:,:,0] *= 1 / float(bins_take)
			self.conditions_trials_ampl[:,:,1:] *= 2 / float(bins_take)
			
			
			# Mean amplitudes and bootstrapped CI_intervals            
			self.conditions_ampl_mean = np.mean(self.conditions_trials_ampl, axis=1)
			
			CI_pos, CI_neg = bootstrap(self.conditions_trials_ampl, alpha=alpha, 
									   n_bootstrap=n_bootstrap)
			self.conditions_ampl_CI_pos, self.conditions_ampl_CI_neg = CI_pos, CI_neg
			
			# isolating F0, F1, and F2 responses and compiling into a single table.
			conditions_f0 = self.conditions_ampl_mean[:,0]
			conditions_f1 = self.conditions_ampl_mean[:,factor]
			conditions_f2 = self.conditions_ampl_mean[:,2*factor]
			
			# Condition Tuning array
			self.cond_tuning = np.vstack((self.conditions,
										 conditions_f0, CI_pos[:,0], CI_neg[:,0],
										 conditions_f1, CI_pos[:,factor], CI_neg[:,factor],
										 conditions_f2, CI_pos[:,2*factor], CI_neg[:,2*factor],
										 conditions_f1/conditions_f0))
			
			# Convert to Hertz - design choice is to keep all PSTH datasets as raw average spike
			# counts, with easy option of seeing frequency in the plotting, but converting to 
			# Hertz for all condition tuning data.
			
			self.cond_tuning[1:-1,:] *= (1/self.bin_width)
			
			# Column labels for pd.dataframe of tuning data
			# ci_perc = (100 * (1 - self.parameters['fft_alpha']))
			idx = ['conditions', 
				   'F0', 'F0_pos_CI', 'F0_neg_CI', 
				   'F1', 'F1_pos_CI', 'F1_neg_CI',
				   'F2', 'F2_pos_CI', 'F2_neg_CI',
				   'F1/F0_ratio']
			# transpose of tuning array to data frame object       
			self.cond_tuning_pd = pd.DataFrame(self.cond_tuning.transpose(), columns=idx)
		
		
			 
		# for orientation data, the orientation angles can get scrambled due to the circ() function
		# rotating the angles around.  This orders them numerically in the final cond_tuning
		
		if self.parameters['condition_type'] == 'orientation':
			self.cond_tuning = self.cond_tuning[:,self.cond_tuning[0].argsort()]
			self.cond_tuning_pd.sort_values(self.cond_tuning_pd.columns[0], inplace=True)


		# 
		# cond_tuning cleaning up
		# 

		if biphase_select_resp is not None:
			assert isinstance(biphase_select_resp, int) and biphase_select_resp in [1,2], \
			f'biphase_select_resp ({biphase_select_resp}) must be an integer of 1 or 2'

			assert self.parameters['biphasic'], 'Stimulus not analysed as biphasic'

			# cond tuning array
			cond_tuning_biphase_mask = self.cond_tuning[4,:] == biphase_select_resp
			self.cond_tuning = self.cond_tuning[:, cond_tuning_biphase_mask]

			# cond tuning pandas dataframe
			self.cond_tuning_pd = self.cond_tuning_pd.query('biphas_id == @biphase_select_resp')



		assert hasattr(self, 'CELL_ID'), 'Make Cell ID first'

		self.cond_tuning_pd.insert(0, 'key', self.CELL_KEY)
		self.cond_tuning_pd.set_index('key', inplace=True)

		self.cond_tuning_pd.insert(0, 'cond_type', self.parameters['condition_type'])
		self.cond_tuning_pd.insert(1, 'cond_unit', self.parameters['condition_unit'])


		
		
	def _out(self):
		'''
		Saves key data to CSV files in the directory of the original data files.
		
		Saves experimental parameters and the condition tuning dataset, in Pandas form.
		'''
		directory = self.Output_path # should be pathlib path
		data_label = self.parameters['data_label']
		con_type = self.parameters['condition_type']

		parameters_path = directory / (data_label + '_parameters.csv')
		tuning_path = directory / (data_label + '_' + con_type + '_tuning.csv')
		
		# Save parameters to csv
		param = pd.Series(self.parameters)
		with open(parameters_path, 'w') as f:
			param.to_csv(f)
		
		
		with open(tuning_path, 'w') as f:
			self.cond_tuning_pd.to_csv(f)
			
		## Other files to save etc
			# PSTH trial and means and SDF - pd with multiindexing and excel write?
			# make convenience functions for reading particular files?


	def _save(self):

		"""
		Pickles a class instantiation of `Zeus`
		"""
		
		directory = self.Output_path

		file_name = f'Themis_{self.CELL_ID["experiment"]}_u{self.CELL_ID["unit"]}_c{self.CELL_ID["cell"]}_r{self.CELL_ID["run"]}.pkl'

		save_path = directory / file_name

		# Atomic saving (helpful?)
		temp_path = save_path.with_suffix(save_path.suffix + '.tmp')
		
		self.SavePath = save_path

		
		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(save_path)



			
	def _plot_psth(self, plot_type = 'hist', frequency = True, 
				   smooth = False, smooth_type = 'gauss', sigma = 3, 
				   mov_avg_window = 3, figsize=(15, 8), format=True):

		"""Plots the PSTHs (histograms) for each of the conditions
		
		Parameters
		__________
		plot_type : string
			Type of plot for hte PSTHs.  Either 'hist', 'line', 'sdf'.
			
		frequency : Boolean
			Convert spike counts to spike frequency in Hz
			
		smooth : Boolean
			Plot, on top, a spike density curve derived from gaussian kernel convolution of
			the PSTH bins.
			Applies to both 'hist' and 'line' plot types.  Redundant for 'sdf'.
			
		smooth_type : string
			type of smoothing to employ over the PSTH data.
			'gauss' employs a gaussian kernel of `sigma` standard deviation (in bins)
			'mov_avg' employs a wimple square kernal of `sigma` width (in bins)
			
		sigma : float
			If `smooth_type == "gauss"`: standard deviation of the gaussian kernel 
			used to produce the smoothed curve.
			If `smooth_type == "mov_avg"`: width of the square kernel.
			This parameter is dimesionless.  It refers to the number of bins.  The standard
			deviation or window width in units of time can be calculated from 
			`sigma * self.bin_width`.
			
		mov_avg : Boolean
			Not developed yet.  Alternative to density curve is a simple moving average.
			
		
		Returns
		_______
		Plots : function
			Passes an `ax.plot()` function call.
		
		
		"""
		plot_types = ['hist', 'line', 'sdf']
		assert plot_type.lower() in plot_types, ('plot_type invalid,\n select from %s' %plot_types)
		
		assert smooth_type.lower() in ['gauss', 'mov_avg'], ('Smoothing type invalid, select '
															'either "gauss" or "mov_avg"')
		# Generate smoothed curves                                        
		def mkSmooth():

			if smooth_type.lower() == 'gauss':
				# Gaussian smoothing 
				spd = gaussian_filter1d(self.conditions_hist_mean, sigma, axis=1)
				
			if smooth_type.lower() == 'mov_avg':
				# Moving average smoothing
				spd = sp.ndimage.convolve1d(self.conditions_hist_mean,
											np.ones((sigma))/float(sigma), axis=1) 

			return spd

			

		n_con = self.parameters['conditions']
		bin_width = self.bin_width
		
		cols = 3
		rows = math.trunc(n_con/3.) + (n_con%3.)        
		
				
		fig = plt.figure(figsize = figsize)

				 
		if plot_type.lower() == plot_types[0]:
			
			# stepping the dataset for a fill_between plot by stacking and 
			# flattening to double up the values and bin boundaries
			# this is substantially faster than the bar plot for dataset this
			# large
			
			step_bins = np.column_stack((self.bins, self.bins)).flatten()[1:-1]
			step_hist_mean = np.column_stack((self.conditions_hist_mean.flatten(),
											  self.conditions_hist_mean.flatten()))
			step_hist_mean = step_hist_mean.flatten()
			step_hist_mean = step_hist_mean.reshape((n_con, (self.bins.size-1)*2))
			
			for c in range(n_con):
				ax = fig.add_subplot(rows, cols, c+1)
				
				ax.set_ylim(0, self.conditions_hist_mean.max() * 1.14)       
				ax.set_xlim(0, self.bins[-1])
				
				ax.fill_between(step_bins, step_hist_mean[c], lw=0, 
								facecolor='0.3', zorder=2)
								
				ax.set_title(self.cond_label[c])
				
				

				if smooth:
					spd = mkSmooth()
					ax.plot(self.bins[:-1] + 0.5*bin_width, spd[c], 
							linestyle='-', color='Coral', linewidth=4, alpha=0.8, 
							zorder=3)
					
			# convert ylabl to frequency units
					
				if frequency:
					freq_label = np.round(ax.get_yticks() * (1 / bin_width), 
										  decimals=1)
					ax.set_yticklabels( freq_label)
				
					# ylabel for all left most subplots
				for sub_plt in np.arange(1, rows*cols, cols):
					if sub_plt == (c+1):
						if frequency:
							ax.set_ylabel('Frequency')
						else:
							ax.set_ylabel('Average count')
				
#                plotform(ax)



		elif plot_type.lower() == plot_types[1]:
			
			for c in range(n_con):
				
				pos_stderr = self.conditions_hist_mean + 2*self.conditions_hist_stderr
				neg_stderr = self.conditions_hist_mean - 2*self.conditions_hist_stderr
	
				ax = fig.add_subplot(rows, cols, c+1)
				
				ax.set_ylim(0, pos_stderr.max() * 1.14)       
				ax.set_xlim(0, self.bins[-1])            
				
				
				ax.plot(self.bins[:-1] + 0.5*bin_width, self.conditions_hist_mean[c],
						color='0.28', linewidth=1)
	
				ax.fill_between(self.bins[:-1] + 0.5*bin_width, 
								pos_stderr[c], 
								neg_stderr[c],
								color='0.6', alpha=0.6)
								
				ax.set_title(self.cond_label[c])
				
	
	
				if smooth:
					spd = mkSmooth()
					ax.plot(self.bins[:-1] + 0.5*bin_width, spd[c], 
							linestyle='-', color='Coral', linewidth=3, alpha=0.8, 
							zorder=3)
	
	
					
				# convert ylabl to frequency units
						
				if frequency:
					freq_label = np.round(ax.get_yticks() * (1 / bin_width), 
										  decimals=1)
					ax.set_yticklabels( freq_label)
				
					# ylabel for all left most subplots
				for sub_plt in np.arange(1, rows*cols, cols):
					if sub_plt == (c+1):
						if frequency:
							ax.set_ylabel('Frequency')
						else:
							ax.set_ylabel('Average count')
				

#                plotform(ax)
				
				
				
		elif plot_type.lower() == plot_types[2]:
			
			cols = 3
			rows = math.trunc(n_con/3.) + (n_con%3.)        
			
			
			for c in range(n_con):
				ax = fig.add_subplot(rows, cols, c+1)
				
				ax.set_ylim(0, self.CI_pos.max() * 1.14)       
				ax.set_xlim(0, self.bins[-1])            
				
				
				ax.plot(self.bins[:-1] + 0.5*bin_width, self.spike_dens_func[c],
						lw=2, color='#036eb6', zorder=1)
	
				ax.fill_between(self.bins[:-1] + 0.5*bin_width, 
								self.CI_neg[c], self.CI_pos[c],
								color='#036eb6', alpha=0.3)
								
				ax.set_title(self.cond_label[c])
								
				if frequency:
					freq_label = np.round(ax.get_yticks() * (1 / bin_width), 
											  decimals=1)
					ax.set_yticklabels( freq_label)
					
				for sub_plt in np.arange(1, rows*cols, cols):
					if sub_plt == (c+1):
						if frequency:
							ax.set_ylabel('Frequency')
						else:
							ax.set_ylabel('Average count')
	
		  # bug with this and teh macosx backend      
#        plt.tight_layout()
		plt.subplots_adjust(hspace=0.45)

	def _plot_psth_flat(self, sigma=5, figsize = (15, 8)):
		""" Do all conditions side by side with a small filter
		"""
	
		gaus_filt = sp.ndimage.gaussian_filter1d
		all_resp = gaus_filt(self.conditions_hist_mean.flatten(), sigma)
		
		fig = plt.figure(figsize=figsize)
		ax = fig.add_subplot(1, 1, 1)
		
		ax.plot(all_resp, linestyle='-', color='0.28')
		
		n_con = self.parameters['conditions']
		con_mark = np.arange(0, (self.bins.size -1) * n_con, self.bins.size -1)
				
		ax.xaxis.set_ticks(con_mark)

		try:
			ax.xaxis.set_ticklabels(self.cond_label)
		except:
			ax.xaxis.set_ticklabels(np.unique(self.marker_codes))
		
		freq_label = np.round(ax.get_yticks() * (1/self.bin_width),
							  decimals = 1)
		ax.set_yticklabels(freq_label)
		ax.set_ylabel('Frequency')
		
		for label in ax.xaxis.get_majorticklabels():
			label.set_horizontalalignment('left')
			
		ax.set_xlim(0, (self.bins.size -1) * n_con)
		
		# bug with macosx backend
#        plt.tight_layout()
		plt.subplots_adjust(hspace=0.45)
		
			
	
	def _plot_tuning(self, plot_type = 'linear', frequency=True, modulation='both',
					 xaxis_scale = 'linear', xaxis_base = '2',
					 figsize = (10, 6)):
		
		"""
		
		Parameters
		__________
		
		modulation : string
			For grating stimuli, which modulation to plot.
			Either 'F1', 'F0' or 'both'.  Default is 'both'.
			
		
		"""
		
		def wrap(array):
			""" 
			returns wrapped version of array, with the first element
			in each rwo appended to the end of the same rows.
			
			For polar plotting or angles/orientations
			"""
			
			assert array.ndim == 2, ('works on 2D arrays only; expecting '
									 'self.cond_tuning, which should be 2D')
											
			return np.column_stack((array, array[..., 0]))
		

		assert plot_type.lower() in ['polar', 'linear'], ('plot type unrecognised.  Use either '
														' "polar" or "linear"')
														
		assert xaxis_scale.lower() in ['linear', 'log'], ('xaxis_scale must be either "linear" or "log"')                                                

		if plot_type.lower() == 'polar':
			assert self.conditions.min()>=0. and self.conditions.max()<=360. , ('Conditions do not'
																			'fall within a circle (in degs)')
																		
			assert xaxis_scale.lower() == 'linear', ('Log scale for a polar plot ...?')

															
		
		xaxis_scale_kw = {'basex':float(xaxis_base)}
		
		
		#==============================================================================
		# Extracting relevant data on the basis of stimulus type (bar v grating)  
		#==============================================================================
		# Defining data set according to which harmonics are desired
		if self.parameters['stimulus'] == 'grating':

			print('grat')
			
			if modulation.lower() == 'f0':
				data = self.cond_tuning[(0, 1, 2, 3),:]
			
			if modulation.lower() == 'f1':
				data = self.cond_tuning[(0, 4, 5, 6),:]

			if modulation.lower() == 'both':
				data = self.cond_tuning[:7,:]
		
		# take whole data set for bar stimuli        
		if self.parameters['stimulus'] == 'bar':
			data = self.cond_tuning


		#==============================================================================
		# Plot Orientation tuning
		#==============================================================================
				
		if self.parameters['condition_type'] == 'orientation':

		
			if plot_type.lower() == 'polar':
				
				fig = plt.figure(figsize=figsize)                
				ax=fig.add_subplot(111, polar=True)
				
				data = wrap(data)
				
				ax.plot(np.radians(data[0]), data[1], color='#507df5', lw=2)
				
				ax.fill_between(np.radians(data[0]), data[2], data[3], 
										   color='#507df5', alpha=0.25, lw=0)
															  
			
			elif plot_type.lower() == 'linear':
				
				fig = plt.figure(figsize=figsize)
				ax = fig.add_subplot(111)
								
				
				ax.plot(data[0], data[1], color='#507df5', lw=2)
				ax.fill_between(data[0], data[2], data[3], color='#507df5',
								alpha=0.25, lw=0)
								
				ax.set_xlabel(self.parameters['condition_type'].capitalize()+' '+self.parameters['condition_unit'])
					

		#==============================================================================
		# Plot Chromatic Tuning
		#==============================================================================

		elif self.parameters['condition_type'] == 'chromatic':
		
			if plot_type.lower() == 'polar':
				
				fig = plt.figure(figsize=figsize)                
				ax=fig.add_subplot(111, polar=True)
				
				data = wrap(data)
				
				ax.plot(np.radians(data[0]), data[1], color='#507df5', lw=2)
				
				ax.fill_between(np.radians(data[0]), data[2], data[3], 
										   color='#507df5', alpha=0.25, lw=0)
															  
			
			elif plot_type.lower() == 'linear':
				
				fig = plt.figure(figsize=figsize)
				ax = fig.add_subplot(111)
								
				
				ax.plot(data[0], data[1], color='#507df5', lw=2)
				ax.fill_between(data[0], data[2], data[3], color='#507df5',
								alpha=0.25, lw=0)
								
				ax.set_xlabel(self.parameters['condition_type'].capitalize()+' '+self.parameters['condition_unit'])



		# ==============================================================================
		# Plot others, including grating fourier harmonics
		# ==============================================================================


				
		else:

			fig = plt.figure(figsize=figsize)
			ax = fig.add_subplot(111)
			
			if plot_type.lower() == 'polar':
				
				fig = plt.figure(figsize=figsize)                
				ax=fig.add_subplot(111, polar=True)
				
				data = wrap(data)
				
				ax.plot(np.radians(data[0]), data[1], color='#507df5', lw=2)
				
				ax.fill_between(np.radians(data[0]), data[2], data[3], 
										   color='#507df5', alpha=0.25, lw=0)
			
						
			
			elif data.shape[0] == 7:# ie, if fourier data with f1 and f0
				ax.plot(data[0], data[1], color='#507df5', lw=2, label='F0')
				ax.fill_between(data[0], data[2], data[3], color='#507df5',
								alpha=0.25, lw=0)
								
				ax.plot(data[0], data[4], color='#bc5692', lw=2, label='F1')
				ax.fill_between(data[0], data[5], data[6], color='#bc5692',
								alpha=0.25, lw=0)
				ax.legend()


			else:

				# for biphasic/bar stimuli that are not orietnation tuning curves - plot both directions
				if self.parameters['biphasic'] == True and self.parameters['condition_type'] != 'orientation':
					id1 = data[-1] == 1
					id2 = data[-1] == 2

					ax.plot(data[0][id1], data[1][id1], color='#507df5', lw=2, label='first')
					ax.fill_between(data[0][id1], data[2][id1], data[3][id1], color='#507df5',
								alpha=0.25, lw=0)


					ax.plot(data[0][id2], data[1][id2], color='#bc5692', lw=2, label='second')
					ax.fill_between(data[0][id2], data[2][id2], data[3][id2], color='#bc5692',
								alpha=0.25, lw=0)

					ax.legend()


				else:
					ax.plot(data[0], data[1], color='#507df5', lw=2)
					ax.fill_between(data[0], data[2], data[3], color='#507df5',
									alpha=0.25, lw=0)
								
			ax.set_xlabel(self.parameters['condition_type'].capitalize()+' '+self.parameters['condition_unit'])
			ax.set_xscale(xaxis_scale, **xaxis_scale_kw)
			
			
		if plot_type.lower() == 'linear':
			ax.set_ylabel('Response (Hz)')
			
			

			

				
## spontaneous rate
				# plotting - subtract or line?
				# handling a blank condition - presume always last?
				# Three ways - blank condition, or blank portion of stimulus, or, pre/post-conditions

## Post-Marker Buffer
	# Unsure of implementation.  Perhaps at analysis, just take only a certain portion of PSTH.
	# Good to keep all data that we can.                
				
## Multiple sources for tuning or just sdf?  If multiple, need cleaner means of management.
#                Options include, raw mean, moving average.  All should be dealt with in the same way.
				
## Saving
				
## Bokeh plotting for inside notebook
	# Option within standard plotting functions, or alternative functions?
	# Probably option within standard functions.  Just for PSTH and tuning.

## in relation to Pandas integration ... general I/O.

## convenience function for listing all attributes

## Conditions / Trials load (for when marker issue doesn't exist)
## Pandas Integration -> question of what to provide a pandas form for
# key advantage will be outputting to excel or csv easily
# possibly provide functions for viewing and excel purposes that provide
# returns and not attributes

## ATHENA!!! (analysis!!)
# Orientation Tuning  Basic scheme - olympos (project), Athena (unit), Zeus(particular run)
# Spike Shape (simple), maybe a zeus function

# np.convolve()