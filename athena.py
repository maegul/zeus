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
# from . import themis




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


	def save(self, file_name = None):

		if file_name is None:

			if not hasattr(self, 'SavePath'):
				file_name = self.path / ('Athena_' + self.PROJ_ID['name'] + '_' + self.PROJ_ID['ID'] + '.pkl')

				self.SavePath = file_name
				self._absolute_paths['SavePath'] = file_name.absolute()

			else:
				file_name = self.SavePath

		else:
			file_name = pthl.Path(file_name)
			self.SavePath = file_name
			self._absolute_paths['SavePath'] = file_name.absolute()


		# Atomic saving (helpful?)
		temp_path = file_name.with_suffix(file_name.suffix + '.tmp')

		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(file_name)

		print(f'Saved pickle file to {str(self.path)} as {str(file_name.name)}')


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

			self.CellData.sort_index()

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


		# ensuring cond_type is first column
		cd_cols = self.CellData.columns.tolist()
		if cd_cols[0] != 'cond_type':
			cd_cols.insert(0, cd_cols.pop( cd_cols.index('cond_type')))
			self.CellData = self.CellData.reindex(cd_cols)

		# join condition and cell data
		# Relies on the cell key being the index for the join
		self.Data = self.TunData.join(self.CellData)


		self.save()
		
		themis_obj._save()



	def replaceThemis(self, themis_obj):

		# level 1 for run_key in the multi index
		self.TunData.drop(labels = themis_obj.RUN_KEY, level=1, inplace = True)
		self.CellData.drop(labels = themis_obj.RUN_KEY, level=1, inplace = True)

		self.addThemis(themis_obj)



	def getUniqueKeys(self, filtered_data, return_multi_index=False):

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



