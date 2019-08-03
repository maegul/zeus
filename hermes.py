import numpy as np
import pandas as pd

from time import time
from fuzzywuzzy import fuzz, process

# from collections import OrderedDict as OD
import inspect
import pathlib as pathl
import pickle

# def sortDict(dictionary):
# 	return OD(sorted(dictionary.items()))


def dictFromArgs(args, includeNone = True):


	argsDict = {
		k : args.locals[k]
		for k in args.args
		if (args.locals[k] is not None) or includeNone
	}

	# Add kwargs
	# Check if kwargs provided for in function

	if (args.keywords is not None) and (len(args.locals[args.keywords]) > 0):
		for k,v in args.locals[args.keywords].items():
			argsDict[k] = v




	return argsDict


def mkUID():
	'''
	Generates Unique ID from system time
	Returns hex of system epoch time
	'''

	return hex(int(time()*1e7))[2:]


def mk_project_ID(name = None, description = None):

	assert isinstance(name, str) and isinstance(description, str), (
		'Must provide a name and a description as strings'
		)
	args = inspect.getargvalues(inspect.currentframe())

	proj_id = dictFromArgs(args)

	proj_id['ID'] = mkUID()

	return proj_id


def mk_themis_file_name(themis_obj = None, 
	experiment = None, unit = None, cell = None, run = None):

	if themis_obj is not None:

		assert hasattr(themis_obj, 'CELL_ID'), 'Themis object must have an instantiated CELL ID'

		file_name = f'Themis_{themis_obj.CELL_ID["experiment"]}_' +\
					f'u{themis_obj.CELL_ID["unit"]}_' +\
					f'c{themis_obj.CELL_ID["cell"]}_' +\
					f'r{themis_obj.CELL_ID["run"]}.pkl'

	else:
		assert None not in (experiment, unit, cell, run), 'Must pass experiment ... run variables if themis_obj is None'

		file_name = f'Themis_{experiment}_' +\
					f'u{unit}_' +\
					f'c{cell}_' +\
					f'r{run}.pkl'


	return file_name	


def mk_track_file_name(track_obj = None, experiment = None, track = None):

	if track_obj is not None:
		assert hasattr(track_obj, 'experiment'), 'Track object has no experiment attribute'
		assert hasattr(track_obj, 'track'), 'Track object has no track attribute'

		file_name = f'Track_{track_obj.experiment}_t{track_obj.track}.pkl'

	else:
		assert None not in (experiment, track), 'Must pass experiment and track variables if track_obj is None'

		file_name = f'Track_{experiment}_t{track}.pkl'

	return file_name



def show_info(d, spec_keys=None):
	'''
	spec_keys : list
		List of keys, in the order in which the dictionary values are to be printed
		Each key in this list must be in the keys of the provided dictionary
	'''

	if spec_keys is not None:

		for sk in spec_keys:
			assert sk in d.keys(), f'{sk} in spec_keys not in provided info dictionary'
		ks = spec_keys

	else:
		ks = sorted(d.keys())

	lengths = [len(k) for k in d.keys()]
	max_len = max(lengths)
	tab = 1 + max_len//8
	
	for k in ks:
		space = '\t' * (tab - len(k)//8)
		print(f'{k}:{space}{d[k]}')


def mk_cell_ID(experiment = None, unit = None, cell = None, run = None):

	args = inspect.getargvalues(inspect.currentframe())

	for k,v in args.locals.items():
		# Convert numbers - int and float - to strings
		if isinstance(v, (int, float)):
			v = str(v)
			args.locals[k] = v
		assert isinstance(v, str), f'Entry for {k} must be a string, not {v}'

	cell_id = dictFromArgs(args)

	return cell_id


def mk_cell_key(experiment = None, unit = None, cell = None, run = None,
				pureCellKey = False
	):

	assert type(pureCellKey) == bool, 'pureCellKey must be boolean'

	if not pureCellKey:
		key = f'{experiment}u{unit}c{cell}r{run}'

	elif pureCellKey:
		key = f'{experiment}u{unit}c{cell}'

	return key


def mk_cell_key_from_iterable(iterable, pureCellKey = False):
	'''
	Presumes iterable contains strings in order of exp, unit, cell, run
	'''

	if not pureCellKey:
		key = f'{iterable[0]}u{iterable[1]}c{iterable[2]}r{iterable[3]}'
	if pureCellKey:
		key = f'{iterable[0]}u{iterable[1]}c{iterable[2]}'

	return key


def mk_stim_params(**kwargs):

	args = inspect.getargvalues(inspect.currentframe())

	# Deal with defined args first
	for k in args.args:
		v = args.locals[k]
		# Convert numbers - int and float - to strings
		if isinstance(v, (int, float)):
			v = str(v)
			args.locals[k] = v
		assert isinstance(v, str), f'Entry for "{k}" must be a string, not "{v}"'

	if args.keywords is not None:
		for k,v in args.locals[args.keywords].items():
			# Convert numbers - int and float - to strings
			if isinstance(v, (int, float)):
				v = str(v)
				args.locals[args.keywords][k] = v
			assert isinstance(v, str), f'Entry for "{k}" must be a string, not "{v}"'


	stim_params = dictFromArgs(args, includeNone=False)

	return stim_params



# Init dataset - make a dataframe from a single entry
def initDataSet(cell_key, run_key, cell_data, stim_data, cond_type):

	cell_data = cell_data.copy()
	cell_data.update(stim_data.copy())
	cell_data.update({'cond_type': cond_type})

	midx = pd.MultiIndex.from_product([[cell_key],[run_key]],
		names = ['cell_key', 'run_key']
		)
	dataSet = pd.DataFrame(data=cell_data, index=midx)

	return dataSet


def separateDataKeys(dataSet):
	'''
	Generate keys from dataset for both cell data and stim data
	'''

	all_columns = set(
		dataSet.columns.tolist()
		)

	cell_data_keys = set(
		inspect.getcallargs(mk_cell_ID).keys()
		)

	stim_data_keys = all_columns - cell_data_keys

	return cell_data_keys, stim_data_keys



def appendDataSet(cell_key, run_key, cell_data, stim_data, cond_type, dataSet, 
	force=False):


	# Separate Cell data from rest for dataset
	cell_data_keys, stim_data_keys = separateDataKeys(dataSet)


	if not force:
		cell_data_allGood = checkCellData(cell_data, dataSet[list(cell_data_keys)])

		stim_data_allGood = checkStimParams(stim_data, dataSet[list(stim_data_keys)])

		if not (cell_data_allGood and stim_data_allGood):

			print('**** force addition? ****')

			return

	# Prevent mutation of themis attributes!
	cell_data = cell_data.copy()
	cell_data.update(stim_data.copy())
	cell_data.update({'cond_type':cond_type})

	midx = pd.MultiIndex.from_product([[cell_key],[run_key]])
	cell_data = pd.DataFrame(data=cell_data, index=midx)

	# Relying on integrity check to prevent data being added twice
	dataSet = dataSet.append(cell_data, sort=False, verify_integrity=True)

	
	return dataSet	


def initAnalysisData(analData, key):

	# IE, run_keys are keys of the dictionary
	if key is None:
		analDataDF = pd.DataFrame.from_dict(analData, orient='index')

	if isinstance(key, str):
		analDataDF = pd.DataFrame(data=analData, index=[key])

	return analDataDF



def appendAnalysisData(analData, key, dataset):

	if key is None:
		analData = pd.DataFrame.from_dict(analData, orient='index')

	if isinstance(key, str):
		analData = pd.DataFrame(data=analData, index=[key])


	dataset = dataset.append(analData, sort=False, verify_integrity=True)

	return dataset


def checkEntries(query, choices):
	'''
	General wrapper for using fuzzy search of an entry on a list of choices

	Parameters
	----
	query : str
		Entry that is intended to either be new to choices or to match one value
	choices : list
		Previous entries to compare query against

	Returns
	----
	Process_results : list
		list of tuples (2 items).
		Each tuple contains the result of the match and its relative score
		[('match', score), ('match2', score2), ...]
	'''

	proc_results = process.extract(query, choices, scorer = fuzz.partial_ratio)

	return proc_results


# Check entries for cell id
def checkCellData(cell_data, cellDataSet, splitDataSet = False):
	'''
	Checking the entries of the cell data for uniqueness and redundancy
	'''

	if splitDataSet:
		cell_data_keys, _ = separateDataKeys(cellDataSet)

		cellDataSet = cellDataSet[list(cell_data_keys)]


	allGood = True

	for k,v in cell_data.items():

		prev_entries = cellDataSet[k].values

		if v not in prev_entries:

			allGood = False

			proc_results = process.extract(v, prev_entries, scorer=fuzz.partial_ratio)

			print(f'Entry of "{v}" is new to field "{k}" in cell data')
			print(f'Extant alternatives: {[pr[0] for pr in proc_results]}')
			print('\n')

	return allGood




# check fields for stim_params
def checkStimParams(stim_params, stimDataSet, splitDataSet = False):
	'''
	Checking that the fields or keys of new stim params are not unnecessarily unique
	'''

	if splitDataSet:
		_, stim_data_keys = separateDataKeys(stimDataSet)

		stimDataSet = stimDataSet[list(stim_data_keys)]


	allGood = True

	for k in stim_params:
		prev_fields = stimDataSet.columns
		if k not in prev_fields:
			
			allGood = False

			proc_results = process.extract(k, prev_fields, scorer=fuzz.partial_ratio)


			print(f'Key/Parameter of "{k}" is new to the stim parameters')
			print(f'Extant alternatives: {[pr[0] for pr in proc_results]}')
			print('\n')


	return allGood


def loadTrack(file_path):
	'''
	Unpickles a track object (hermes.Track)

	parameters
	----

	file_path : str
		file path of the pickled file	
	'''

	assert isinstance(file_path, str), 'File path must be a string'

	with open(file_path, 'rb') as f:
		return pickle.load(f)


class Track:

	'''
	Object for track reconstruction
	'''
	
	def __init__(self, experiment, track,
				lesion_depths, angle,
				shrinkage, px_scale, slice_width=50):
		'''
		
		Parameters
		----
		
		shrinkage : float
			percentage of linear distance LOST through
			processing.
			Ie, histology_dist = real_dist(1-shrinkage)
			real_dist = histology_dist / (1-s)
		
		pxl_scale : float
			How many pixels per micrometer (um) of primary 
			histology imagery.
			Eg, 0.3 means 0.3px per um or 3.3um per px
		'''

	#     print(inspect.getargvalues(inspect.currentframe()))
	#     print(locals())

		input_vars = {}
		input_arg_info = inspect.getargvalues(inspect.currentframe())

		# Take only input arguments prescribed by function definition
		input_vars.update(
			{a : input_arg_info.locals[a] 
			 for a in input_arg_info.args 
			 if a != 'self'}
		)
		
		self.__dict__.update(input_vars)

		# For saving any save paths
		self._absolute_paths = {}
		
	def setUnits(self, depths, unit_no):
		
		unit_depths = np.array(depths)
		
		# This sorting defines the sorting of all subsequent calculations
		sort_idx = np.argsort(unit_depths)
		
		self.unit_depths = unit_depths[sort_idx]
		self.unit_nos = np.array(unit_no)[sort_idx]
		
		
	def _hist_dist(self, real_dist):
		return real_dist * (1-self.shrinkage)
	
	def _real_dist(self, hist_dist):
		return hist_dist / (1-self.shrinkage)
	
	def _px_dist(self, hist_dist):
		return hist_dist * self.px_scale
	
	def _hist_from_pxl_dist(self, pxls):
		return pxls / self.px_scale
	
	def _cos(self, theta):
		return np.cos(np.radians(theta))
	
	def _sin(self, theta):
		return np.sin(np.radians(theta))
		
		
	def calcDepthValues(self):
		
		depth_values = {}
		# Take the deepest lesion as the reference point
		reference_depth = self.lesion_depths[-1]
		
		depth_values['rel_depth_real'] = self.unit_depths - reference_depth
		
		depth_values['rel_depth_hist'] = self._hist_dist(depth_values['rel_depth_real'])
		
		depth_values['rel_depth_hist_cut_plane'] = (
			depth_values['rel_depth_hist'] * self._cos(self.angle)
		)
		
		depth_values['rel_depth_cut_plane_pxl'] = (
			self._px_dist(depth_values['rel_depth_hist_cut_plane'])
		)
		
		depth_values['rel_sag_trans'] = np.abs(
			depth_values['rel_depth_hist'] * self._sin(self.angle)
		)
		
		depth_values['rel_sag_n_slice'] = np.abs(
			np.round(
			depth_values['rel_sag_trans'] / self.slice_width,
			decimals=2
		)
		)
		
		self.reference_depth = reference_depth
		self.depth_values = depth_values
		
		self.df_depth_values = pd.DataFrame.from_dict(self.depth_values)
		self.df_depth_values['unit'] = self.unit_nos
		self.df_depth_values = self.df_depth_values.set_index('unit')
		
		
	def setLayerDepths(self, borders, layers):
		'''
		Define layer borders, heading INTO cortex (from layer 1)
		
		Parameters
		---
		borders : list 
			define when the beginning of the relevant border is.
			values are presumed to be in pixels of the primary histology
			images, sharing the same px_scale as provided in self.__init__
			
			Relative to the last lesion location, and in cutting plane.
			
			MUST be in same order as layers - user discretion.
			AND, in numerical order (ie smallest depth first)
			
		layers : list[int]
			provide labels for the relevant borders, in ints
			MUST be in same order as borders.
			
			Will be converted to strings for labelling purposes
		'''
		
		# Check that borders are floats
		for b in borders:
			assert isinstance(b, (int,float)), 'borders must be floats'
			
		# check that borders are ordered
		for bi in range(len(borders)-1):
			assert borders[bi] < borders[bi+1], f'borders should be ordered, entries {bi} and {bi+1} are not'
			
		layer_labels = [
			str(l)
			for l in layers
		]
		
		self.layer_labels = layer_labels
		self.layer_borders = borders
			
			
		# Calculate layers each unit belongs in (where possible)
		
		unit_layer_labels = [None] * len(self.unit_depths)
		
		for ui, ud in enumerate( self.depth_values['rel_depth_cut_plane_pxl'] ):
			for li, lb in enumerate(self.layer_borders):
				
				if ud > lb:
					unit_layer_labels[ui] = self.layer_labels[li]
					
		self.unit_layer_labels = unit_layer_labels
		
		
		# Calculate relative position of unit within layer (where possible)
		
		unit_rel_layer_depth = [None] * len(self.unit_depths)
		
		for ui, ul in enumerate( self.unit_layer_labels):
			
			# unit has a layer label and that layer has known end
			if (ul != None) and ( ul != self.layer_labels[-1] ):
				
				lay_idx = self.layer_labels.index(ul)
				layer_beg, layer_end = self.layer_borders[lay_idx : lay_idx+2]
				
				# unit's location, relative to beginning of layer, normalised to length of layer
				location = round(
					(self.depth_values['rel_depth_cut_plane_pxl'][ui] - layer_beg) / 
					(layer_end - layer_beg),
					3
				)
				
				unit_rel_layer_depth[ui] = location
				
		self.unit_rel_layer_depth = unit_rel_layer_depth
		
		self.df_depth_values['layer'] = self.unit_layer_labels
		self.df_depth_values['rel_layer_depth'] = self.unit_rel_layer_depth


	def save(self):

		# Kinda given up on special path changing abilities and enforcement
		# This object gets saved in the current, path, that simple
		directory = pathl.Path('.')
		file_name = mk_track_file_name(track_obj=self)

		save_path = directory / file_name

		self._absolute_paths['save_path'] = save_path.absolute()
		self.save_path = save_path

		temp_path = save_path.with_suffix(save_path.suffix + '.tmp')

		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(save_path)

		print(f'Saved track object as {save_path}')




	# save
		# pickle, as usual, for whole object
		# BUT also, csv s of units and lesion
		# save in histology folder
