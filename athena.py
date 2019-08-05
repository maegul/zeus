import numpy as np
import pandas as pd

# from circ_stat import circ_stat as cs

import pickle
import inspect
# import glob
# import copy

import pathlib as pthl


# import collections

from . import hermes
from . import curveFuncs
from . import themis




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
		athena_obj = pickle.load(f)

		athena_obj._absolute_paths['previous_load_path'] = pthl.Path(file_path)

		return athena_obj
		
		

class Athena:

	def __repr__(self):

		summary = self.CellData.loc[:,['experiment', 'unit', 'cell', 'run']].groupby(by=['experiment', 'unit', 'cell']).count()

		# Number of unique indices from the multiindex of the grouping (experiment, unit, cell) ... ie, how many unique cells
		n_cells = summary.size

		n_animals = self.CellData.experiment.unique().size

		rep = f''

		rep += f"{self.PROJ_ID['name']}:\t{self.PROJ_ID['description']}\n"
		rep += f"{n_animals} Animals   |   {n_cells} Cells\n"
		rep += f"{summary}"

		return rep


	def __init__(self, name = None, description = None, path = '.'):

		'''
		path default is current directory.  Idea is all zeus files and athena files will be in same directory

		Decoupling of data and output path can be done at zeus level
		'''

		self.PROJ_ID = hermes.mk_project_ID(name = name, description = description)

		self.path = pthl.Path(path)

		# save absolute paths if necessary
		self._absolute_paths = {
			'path': self.path.absolute()
		}


	def save(self, file_name = None, 
				use_previous_load = True, use_absolute = False):

		'''
		Parameters
		----

		use_previous_load : boolean (True)
			on load, path used is saved to the athena obj before return
			under _absolute_paths['previous_load_path'].
			If True, use this path as the current save path.
			Helpful when athena_obj has been loaded into an environment
			with a different working directory from that of the object's
			source.

		use_absolute : boolean (False)
			If object has previously been saved (ie, hasattr('SavePath')),
			then use the absolute path attribute as a quick way to save
			to the original location of the saved object

			Particularly important when adding data to the athena object with
			functions like addTrack or addThemis, which then automatically save
			the athena object.
		'''

		if file_name is None:

			if not hasattr(self, 'SavePath'):
				file_name = self.path / ('Athena_' + self.PROJ_ID['name'] + '_' + self.PROJ_ID['ID'] + '.pkl')

				self.SavePath = file_name
				self._absolute_paths['SavePath'] = file_name.absolute()

			elif use_previous_load:
				file_name = self._absolute_paths['previous_load_path']

			elif use_absolute:
				file_name = self._absolute_paths['SavePath']

		else:
			file_name = pthl.Path(file_name)
			self.SavePath = file_name
			self._absolute_paths['SavePath'] = file_name.absolute()


		# Atomic saving (helpful?)
		temp_path = file_name.with_suffix(file_name.suffix + '.tmp')

		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(file_name)

		print(f'Saved pickle file to {str(file_name.parent)} as {str(file_name.name)}')


	def addThemis(self, themis_obj, check_redundancy=False):

		# Testing themis_obj is appropriate

		if isinstance(themis_obj, str):

			# Not sure how to use or store the path

			# themis_path = themis_obj
			with open(themis_obj, 'rb') as f:
				themis_obj = pickle.load(f)

		assert themis_obj.__class__.__name__ == 'Themis', 'Themis object is not a themis object?'


		hasCellID = hasattr(themis_obj, 'CELL_ID')
		hasStimParams = hasattr(themis_obj, 'STIM_PARAMS')

		assert hasCellID and hasStimParams, (
			f'Themis object does not have required data sets.\n'
			f'Cell ID: {hasCellID};  Stim Params: {hasStimParams}'
			)

		###
		# Cell Data init or append
		###

		# Init cell dataset as it does not exist yet
		if not hasattr(self, 'CellData'):
			self.CellData = hermes.initDataSet(
				themis_obj.CELL_KEY, themis_obj.RUN_KEY,
				themis_obj.CELL_ID, themis_obj.STIM_PARAMS,
				themis_obj.parameters['condition_type']
				)

		# Appened cell dataset
		else:
			assert themis_obj.RUN_KEY not in self.CellData.index[1], (
				f'Themis with RUN_KEY {themis_obj.RUN_KEY} already in athena.CellData'
				)
			
			self.CellData = hermes.appendDataSet(
				themis_obj.CELL_KEY, themis_obj.RUN_KEY,
				themis_obj.CELL_ID, themis_obj.STIM_PARAMS, 
				themis_obj.parameters['condition_type'],
				self.CellData,
				# hermes.appendDataSet has checking built in, if force is False
				force = (not check_redundancy)) 

		self.CellData.sort_index(inplace=True)

		###
		# Tuning Data init or append
		###

		# Take tuning data - init main Dataset
		if not hasattr(self, 'TunData'):
			self.TunData = themis_obj.cond_tuning_pd.copy()

		# Or, concat with existing
		else:
			assert themis_obj.RUN_KEY not in self.TunData.index[1], (
				f'Themis with RUN_KEY {themis_obj.RUN_KEY} already in athena.CellData'
				)

			self.TunData = pd.concat([self.TunData, themis_obj.cond_tuning_pd],
				sort = True)

		self.TunData.sort_index(inplace=True)


		# ensuring cond_type is first column
		cd_cols = self.CellData.columns.tolist()
		if cd_cols[0] != 'cond_type':
			cd_cols.insert(0, cd_cols.pop( cd_cols.index('cond_type')))
			self.CellData = self.CellData.reindex(cd_cols)

		# join condition and cell data
		# Relies on the cell key being the index for the join
		# Avoid overlapping columns by joining only those columns that are not present
		# in CellData

		cols_to_use = self.TunData.columns.difference(self.CellData.columns)
		self.Data = self.CellData.join(self.TunData[cols_to_use])


		self.save()
		
		themis_obj._save()



	def replaceThemis(self, themis_obj):

		# level 1 for run_key in the multi index
		self.TunData.drop(labels = themis_obj.RUN_KEY, level=1, inplace = True)
		self.CellData.drop(labels = themis_obj.RUN_KEY, level=1, inplace = True)

		self.addThemis(themis_obj)



	def addTrack(self, track):

		'''
		Add data from a track object

		track must be a track object (hermes.Track)
		'''
		
		# Check if track experiment matches current data in CellData
		if track.experiment not in self.CellData.experiment.unique():
			
			proc_results = hermes.checkEntries(track.experiment, self.CellData.experiment.unique())
			
			print(
				f'Track experiment {track.experiment} is not in project data \n'
				f'Close entries: {[p[0] for p in proc_results]} \n'
				'Track data not added'
				 )
			
			return

		extant_units = self.CellData.query('experiment == @track.experiment').unit.unique()

		for un in track.unit_nos:
			un = str(un)
			if un not in extant_units:
				proc_results = hermes.checkEntries(un, extant_units)
				
				print(
					f'Unit {un} no in project data for experiment {track.experiment} \n'
					f'Close entries: {[p[0] for p in proc_results]}'
					'track data not added'
				)
				
				return

		# Assign columns for layer and track information if not already there in Cell Data
		for col in ['layer', 'rel_layer_depth', 'track']:
			if col not in self.CellData.columns:
				self.CellData.assign(**{col: None})
		
		# Iterate through each unit, adding its information to the appropriate column
		for ui, un in enumerate(track.unit_nos):
			
			un = str(un)
			
			assignment_idx = (self.CellData.experiment == track.experiment) & (self.CellData.unit == un)
			
			self.CellData.loc[assignment_idx, 'layer'] = track.unit_layer_labels[ui]
			self.CellData.loc[assignment_idx, 'rel_layer_depth'] = track.unit_rel_layer_depth[ui]
			self.CellData.loc[assignment_idx, 'track'] = int(track.track)

		self.save(use_previous_load=True)




	def addAnalysisData(self, analData, key=None):

		'''
		Parameters
		----
		analData : dict
			Dictionary with results of analysis as items
			Either for a single run, or for multiple runs in a nested dict,
			see docs for key parameter
			
		key : str
			Default: None
			If None, analData presumed to be dict with first level keys
			providing run_keys.
			If key provided, analData is presumed to be a single row of analysis
			data, in a dictionary.	
		'''

		if not hasattr(self, 'AnalysisData'):
			self.AnalysisData = hermes.initAnalysisData(analData, key)

		else:
			self.AnalysisData = hermes.appendAnalysisData(analData, key, self.AnalysisData)

		self.AnalysisData.index.name = 'run_key'
		self.AnalysisData.sort_index()
		self.save()


	def replaceAnalysisData(self, analData, key=None):

		if key is not None:
			assert key in self.AnalysisData.index, \
				f'key {key} not in AnalysisData, cannot be replaced'
			self.AnalysisData.drop(index=key, inplace=True)
			self.AnalysisData = self.addAnalysisData(analData, key, self.AnalysisData)

		elif key is None:
			assert isinstance(analData, dict), \
				'analData must be a dictionary if no key is provided'

			keys = analData.keys()

			for k in keys:
				assert k in self.AnalysisData.index, \
					f'run_key {k} not in index'

			self.AnalysisData.drop(index=keys, inplace=True)
			self.AnalysisData = self.addAnalysisData(analData, key, self.AnalysisData)




	def getUniqueKeys(self, filtered_data, return_multi_index=False):

		'''
		Returns list of strings, each being a unique run key

		Parameters
		----
		filtered_data : dataframe
			A dataframe from an athena object filtered as desired.
			Must contain cell_id information (ie, CellData or Data)

		return_multi_index : boolean
			Switch for whether to return the MultiIndex object for each
			key.
			If so, two objects return, list of strings, and multiindex object
		'''

		idx = filtered_data.groupby(by=['experiment', 'unit', 'cell', 'run']).count().index

		keys = [
			hermes.mk_cell_key_from_iterable(i)
			for i in idx
		]

		if not return_multi_index:
			return keys
		elif return_multi_index:
			return keys, idx

		else:
			return keys


	def getCellIdFromRunKey(self, run_key):

		'''
		Retruns a dict with Cell Id data, as in themis_obj.CELL_ID
		'''

		cell_id = (self.CellData.xs(run_key, level=1)
						.loc[:, ['experiment','unit', 'cell', 'run']]
						.iloc[0] # dataframe to series by selecting row
						.to_dict()
					)

		return cell_id


	def getData(self, exp=None, u=None, c=None, r=None, key = None,
		columns = ['condition', 'max_resp']
		):
		'''
		For given cell metadata or key return the full data 
		(including cell,stim and tuning)

		Key must be run specific, ie, RUN_KEY

		Returns
		----
		Slice of self.Data as an array, 
		specific to cell metadata/key, and columns argument.

		Data is transposed into columnar organisation
		Ie, column names in first row with data in second row beneath
		'''
		
				
		if key is None:

			exp = str(exp)
			u = str(u)
			c = str(c)
			r = str(r)

			data = self.Data.query('(experiment == @exp) & (unit == @u) & (cell == @c) & (run == @r)')
			
		elif isinstance(key, str):	
			data = self.Data.loc[(slice(None), key),:]

		
		data = data.loc[:,columns].values.T

		return data

	# def addCurveFunc(self, name, curveFunc, overwrite=False):

	# 	'''
	# 	Curve is expected to be used for curve fitting (with scipy.optimize)
	# 	As such, it is presumed to take only positional arguments
	# 	And, it is expected that the first positional argument is 'x'


	# 	'''

	# 	if not hasattr(self, 'CurveFuncs' ):
	# 		self.CurveFuncs = dict()

	# 	assert isinstance(name, str), 'name must be a string'
	# 	assert callable(curveFunc), 'curveFunc must be a callable/function'

	# 	assert (name not in self.CurveFuncs) and (not overwrite), f'Name ({name}) already in use, overwrite? '


	# 	self.CurveFuncs[name] = curveFunc




	# def getCurveDoc(self, name):

	# 	assert hasattr(self, 'CurveFuncs'), 'CurveFuncs not initialised yet, use self.addCurveFunc'

	# 	print(self.CurveFuncs[name].__doc__)


	# def getCurveCode(self, name, return_string=False):

	# 	assert hasattr(self, 'CurveFuncs'), 'CurveFuncs not initialised yet, use self.addCurveFunc'

	# 	source = inspect.getsource(self.CurveFuncs[name])
	# 	if not return_string:
	# 		print( source  )
	# 		return

	# 	elif return_string:
	# 		return source 


	def addCurveFit(self, cell_data_key, name, popt_args, RSq_val, overwrite=False):

		if not hasattr(self, 'CurveFits'):
			self.CurveFits = dict()

		if cell_data_key in self.CurveFits:	
			assert name not in self.CurveFits[cell_data_key], f'function name ({name}) already in fits for {cell_data_key}'

		if cell_data_key not in self.CurveFits:
			self.CurveFits[cell_data_key] = dict()

		if name not in self.CurveFits[cell_data_key]:
			self.CurveFits[cell_data_key][name] = dict()

		func_args = inspect.getfullargspec(curveFuncs.__dict__[name])[0] # only need positional args
		func_args.remove('x') # don't need the x input arg (should come from tuning conditions)

		self.CurveFits[cell_data_key][name]['popt'] = dict(zip(func_args, popt_args))

		self.CurveFits[cell_data_key][name]['RSq'] = RSq_val



	def addCurveFitSuggestion(self, cell_data_key, name, popt_args, RSq_val, overwrite=False):

		if not hasattr(self, 'CurveFitsSugg'):
			self.CurveFitsSugg = dict()

		if cell_data_key in self.CurveFitsSugg:	
			assert name not in self.CurveFitsSugg[cell_data_key], f'function name ({name}) already in fits for {cell_data_key}'


		if cell_data_key not in self.CurveFitsSugg:
			self.CurveFitsSugg[cell_data_key] = dict()

		if name not in self.CurveFitsSugg[cell_data_key]:
			self.CurveFitsSugg[cell_data_key][name] = dict()

		func_args = inspect.getfullargspec(curveFuncs.__dict__[name])[0] # only need positional args
		func_args.remove('x') # don't need the x input arg (should come from tuning conditions)

		self.CurveFitsSugg[cell_data_key][name]['popt'] = dict(zip(func_args, popt_args))

		self.CurveFitsSugg[cell_data_key][name]['RSq'] = RSq_val


	def getCurve(self, key):
		curve_func_name = sorted(self.CurveFits[key].items(), key = lambda cv: cv[1]['RSq'])[-1][0]
		curve_func = getattr(curveFuncs, curve_func_name)
		
		return curve_func_name, curve_func



	def getCurveData(self, key, n_interp=None):
		'''
		Parameters
		----
		key : str
			Run Key for pulling data out for the specified run	

		n_interp : int
			Default: None
			If not None, determines number of data points to use in an interpolation
			of the conditions data (from min to max)
		
		'''


		curve_func_name, curve_func = self.getCurve(key)
		
		func_args = inspect.getfullargspec(curve_func)[0]
		func_args_param = self.CurveFits[key][curve_func_name]['popt']

		p0 = [
			func_args_param[arg]
			for arg in func_args[1:]
		]    
		
		data = self.getData(key = key)

		if n_interp is None:
			cond_data = data[0]
		else:
			assert isinstance(n_interp, int), 'n_interp must be an integer'
			cond_data = np.linspace(data[0].min(), data[0].max(), n_interp)

		y_curve = curve_func(cond_data, *p0)
		
		return np.vstack((
			cond_data, y_curve
		))


	def genThemisDirectory(self):

		directory = hermes.mk_themis_files_directory(self)

		self.ThemisDirectory = directory


	def getThemisObj(self, 
		run_key = None,
		experiment = None, unit = None, cell = None, run = None):


		'''
		Returns themis object using run_key or cell_id metadata
		'''

		assert hasattr(self, 'ThemisDirectory'), (
			'Project has not generated a ThemisDirectory yet, use self.genThemisDirectory'
			)

		if run_key is None:
			assert None not in (experiment, unit, cell, run), (
				'If not using a run_key, must provide all of experiment, unit, cell and run arguments'
				)

			run_key = hermes.mk_cell_key(*[experiment, unit, cell, run])

		assert self.ThemisDirectory[run_key] is not None, (
			f'Run key {run_key} does not have a path in the directory.\n'
			'Generate the directory again using self.genThemisDirectory to see if it can be found'
			)

		themis_obj = themis.load(self.ThemisDirectory[run_key])

		return themis_obj








def genRSq(xdata, ydata, curveFunc=None, opt_curveFuncArgs=None, curve_vals = None):

	'''
	
	Parameters
	_____

	curveFunc_name : str
		Reference to the dictionary key for a function stored in self.CurveFuncs

	opt_curveFuncArgs : dict, list
		Dict or list of optimised args to be passed to the function found at curveFunc_name
	'''

	if curve_vals is None:
		if isinstance(opt_curveFuncArgs, dict):
			curve = curveFunc(xdata, **opt_curveFuncArgs)
		elif isinstance(opt_curveFuncArgs, list) or hasattr(opt_curveFuncArgs, '__iter__'):
			curve = curveFunc(xdata, *opt_curveFuncArgs)

	else:
		assert hasattr(curve_vals, '__iter__'), 'curve vals must be an iterable'

		assert len(curve_vals) == len(ydata), 'ydata, xdata and curve_vals must have same size or length'

		curve = curve_vals

	
	res = ydata - curve
	
	ss_res = np.sum(res**2)
	
	ss_tot = np.sum(
		(ydata - np.mean(ydata))**2
	)
	
	return 1-(ss_res/ss_tot)




# def getCureveList(self, docs=True):


# 	for c,f in self.CurveFuncs.items():
# 		print(f'{c}{inspect.signature(f).__str__()}')
# 		print(f.__doc__)



