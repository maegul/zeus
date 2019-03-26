# import numpy as np
import pandas as pd

# from circ_stat import circ_stat as cs

import pickle
# import glob
# import copy

import pathlib as pthl


# import collections

from . import hermes
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


	def save(self):

		file_name = self.path / ('Athena_' + self.PROJ_ID['name'] + '_' + self.PROJ_ID['ID'] + '.pkl')

		self.SavePath = file_name
		self._absolute_paths['SavePath'] = file_name.absolute()

		# Atomic saving (helpful?)
		temp_path = file_name.with_suffix(file_name.suffix + '.tmp')

		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(file_name)

		print(f'Saved pickle file to {str(self.path)} as {str(file_name.name)}')


	def addThemis(self, themis_obj, check_redundancy=False):

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

		# Init cell dataset as it does not exist yet
		if not hasattr(self, 'CellData'):
			self.CellData = hermes.initDataSet(
				themis_obj.CELL_KEY, themis_obj.CELL_ID, themis_obj.STIM_PARAMS
				)

		# Appened cell dataset
		else:
			assert themis_obj.CELL_KEY not in self.CellData.index, (
				f'Themis with CELL_KEY {themis_obj.CELL_KEY} already in athena.CellData'
				)
			
			self.CellData = hermes.appendDataSet(
				themis_obj.CELL_KEY,
				themis_obj.CELL_ID, themis_obj.STIM_PARAMS, self.CellData,
				# hermes.appendDataSet has checking built in, if force is False
				force = (not check_redundancy)) 

		# Take tuning data - init main Dataset
		if not hasattr(self, 'TunData'):
			self.TunData = themis_obj.cond_tuning_pd.copy()

		# Or, concat with existing
		else:
			assert themis_obj.CELL_KEY not in self.TunData.index, (
				f'Themis with CELL_KEY {themis_obj.CELL_KEY} already in athena.CellData'
				)

			self.TunData = pd.concat([self.TunData, themis_obj.cond_tuning_pd],
				sort = True)

		# join condition and cell data
		# Relies on the cell key being the index for the join
		self.Data = self.TunData.join(self.CellData)

		self.save()
		
		themis_obj._save()



	def replaceThemis(self, themis_obj):

		self.TunData.drop(labels = themis_obj.CELL_KEY, inplace = True)
		self.CellData.drop(labels = themis_obj.CELL_KEY, inplace = True)

		self.addThemis(themis_obj)
