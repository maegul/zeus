import numpy as np
import pandas as pd

from time import time
from fuzzywuzzy import fuzz, process

# from collections import OrderedDict as OD
import inspect
import pathlib as pthl
import pickle

import os
import shutil
# import glob
import fnmatch

import re

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



# > Project ID

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


# > File Names

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

def mk_nb_file_name(nb_type = None, experiment=None, unit=None, run=None):
	'''
	Return file name for a notebook appropriate for the given arguments

	nb_type : str ('temp' | 'anal')
		temp -> templating notebok	
		anal -> analysis notebook
	'''

	if nb_type == 'temp':
		return f'{experiment}_u{unit}_r{run}.ipynb'.lower()
	elif nb_type == 'anal':
		return f'{experiment}_u{unit}_analysis.ipynb'.lower()



def mk_track_file_name(track_obj = None, experiment = None, track = None):

	if track_obj is not None:
		assert hasattr(track_obj, 'experiment'), 'Track object has no experiment attribute'
		assert hasattr(track_obj, 'track'), 'Track object has no track attribute'

		file_name = f'Track_{track_obj.experiment}_t{track_obj.track}.pkl'

	else:
		assert None not in (experiment, track), 'Must pass experiment and track variables if track_obj is None'

		file_name = f'Track_{experiment}_t{track}.pkl'

	return file_name



# > Make New NBs

def find_project_nb_templates(proj):
	'''
	Searches the path of proj for a folder that contains notebook templates

	Such folder expected to be something like nb_template or notebook_template

	regex: ^(nb|notebook).*(template), IGNORECASE

	returns relative path (relative to proj path) if match
	else, -1
	'''

	nb_check = re.compile(r'^(nb|notebook).*(template)', flags=re.I)

	for root, dirs, files in os.walk(proj._absolute_paths['path']):
		for d in dirs:
			if nb_check.match(d):
				return d

	return -1



# adjust below to have templating and analysis type argument so that one function for both

def mk_new_nb(proj, nb_type = None, experiment=None, unit=None, run=None):
	'''
	Makes a new notebook in the root path of the project provided

	Notebook is copied from a template.  Location of templates must be
	recorded in proj object.

	nb_type determines which template is copied

	Parameters
	----
	nb_type : str (template | analysis)
		What type of notebook to make (and therefore what type of template to copy)	
		temp -> for templating, and for a particular run 
		anal -> for analysing all the runs of a particular unit (run arg igored)
	'''

	template_names = dict(
		temp = pthl.Path('template.ipynb'),
		anal = pthl.Path('analysis.ipynb')
		)

	temp_types_need_run = ['temp']

	assert nb_type in template_names.keys(), (
		f'nb_type ({nb_type}) must be one of {template_names.keys()}'
		)

	assert hasattr(proj, '_NBTemplateDir'), (
						'Project has no nb_templates directory ' 
						'use proj.assignNBTemplateDir() to find, and '
						'if it does not exist, make it')

	assert None not in (experiment, unit), (
		'Assign values to key word arguments'
		)

	assert run is not None and nb_type not in temp_types_need_run, (
		f'If nb_type is in {temp_types_need_run}, must provide run arg'
		)

	assert isinstance(experiment, str), 'experiment must be string'
	experiment = experiment.lower()

	try:
		unit = int(unit)
	except ValueError:
		print(f'argument unit ({unit}) cannot be an integer; must be an integer')
		return -1

	try:
		run = int(run)
	except ValueError:
		print(f'argument run ({run}) cannot be an integer; must be integer')
		return -1


	# template_name = pthl.Path('template.ipynb')
	template_name = template_names[nb_type]
	template_file = proj._NBTemplateDir['abs_path'] / template_name

	new_file_name = mk_nb_file_name(nb_type=nb_type, experiment=experiment, unit=unit, run=run)
	new_nb_file_path = proj._absolute_paths['path'] / new_file_name

	assert template_file.exists(), (
		f'Notebook template of name {template_name} cannot be found '
		f'in NB template dir {proj._NBTemplateDir}'
		f'Create or rename notebook file to {template_name}'
		)

	assert not new_nb_file_path.exists(), (
		f'Templating notebook already exists with file name {new_file_name} '
		f'at location {proj._absolute_paths["path"]}'
		)

	print(f'Copying \n{template_file} to \n{new_nb_file_path} ...')

	shutil.copy2(template_file, new_nb_file_path)

	print('Done')

# > Directories

def mk_nb_files_directory(proj, nb_type = None, return_strays=False):
	'''
	nb_type : str (temp | anal)
		Specifies what type of notebook is to be searched for	
	'''


	nb_type_opts = ['temp', 'anal']
	assert nb_type in nb_type_opts, (
		f'nb_type must be either of {nb_type_opts}, not {nb_type}'
		)

	# define patterns for finding candidates for both nb types
	nb_file_patterns = {
		nbt : mk_nb_file_name(nb_type=nbt).replace('none', '*')
		for nbt in nb_type_opts
	}


	# define functions for generating no_cell keys from run keys
	# as appropriate for both nb types

	# matches anything that is (<exp>u<unit>)(c<cell>)(r<run>)
	no_cell_pat = re.compile(r'(.*u.+)(c.+)(r.+)')

	if nb_type == 'temp':
		no_cell = lambda run_key : no_cell_pat.sub(r'\1\3', run_key)

	elif nb_type == 'anal':
		no_cell = lambda run_key : no_cell_pat.sub(r'\1', run_key)

	nb_files, nb_file_paths, nb_relative_file_paths = [], [], []	

	runs = mk_cell_id_from_many_run_keys('all', proj)

	###
	# Find candidates
	###

	# generate matching pattern by calling file name function, with
	# default args and replacing all the lowercase 'none' with '*' to allow for
	# finding the variations in the directory

	# nb_file_pattern = mk_nb_file_name(nb_type=nb_type).replace('none', '*')
	# nb_file_pattern = mk_nb_file_name(nb_type='temp').replace('none', '*')
	nb_file_pattern = nb_file_patterns[nb_type]
	for root, dirs, files in os.walk(proj._absolute_paths['path']):

		# prune ipynb checkpoints from directories that will be walked
		if '.ipynb_checkpoints' in dirs:
			dirs.remove('.ipynb_checkpoints')

		matches = fnmatch.filter(files, nb_file_pattern)

		nb_files.extend(matches)
		nb_file_paths.extend(
			[
				os.path.join(root, m)
				for m in matches
			]
		)
		nb_relative_file_paths.extend(
			[
				os.path.relpath( os.path.join(root, m), proj._absolute_paths['path'])
				for m in matches
			]
		)

	run_file_paths = {
		k: None
		for k in runs.keys()
	}


	###
	# Try to match candidates to run keys
	###

	# Keys are run keys, without the cell
	# So that runs that differ only by the cell can have their notebooks cached and retrieved
	# as templating notebooks are specific to exp, unit, run (not cells)
	paths_cache = {}


	# go through all run keys
	for r, cell_id in runs.items():

		# convert each run key to a no_cell version (cell info removed)
		nc_run = no_cell(r)

		# check if key in paths_cache
		if nc_run in paths_cache:

			# If so, then take paths object from cache and use this
			run_file_paths[r] = paths_cache[nc_run]

		# If not, then look for the paths in the list and pop them if found
		else:

			# Generate putative file name, to be searched for in the list of candidates
			file_name = mk_nb_file_name(nb_type=nb_type,
				**{arg : cell_id[arg] for arg in ['experiment', 'unit', 'run']}
			)

			try:
				file_idx = nb_files.index(file_name)

				paths_object = dict(
						abs_path = pthl.Path(nb_file_paths.pop(file_idx)),
						rel_path = pthl.Path(nb_relative_file_paths.pop(file_idx))
					)

				run_file_paths[r] = paths_object
				paths_cache[nc_run] = paths_object # store paths object in cache

				nb_files.pop(file_idx) # pop file from file list

			except ValueError:
				pass

	if return_strays:
		return nb_file_paths, nb_relative_file_paths

	else:
		return run_file_paths








def mk_themis_files_directory(proj, return_strays=False):
	'''
	Generates directory of Themis files for a given project

	proj : Athena Obj

	Returns a dict of dicts {run_key : {abs_path : str, rel_path : str}}
	paths are str. abs_path paths are absolute, 
	relative paths are relative to the proj absolute_paths path


	Parameters
	----	
	return_strays : bool
		If true, don't return directory.
		Instead, return list of (absolute) file paths and relative file paths
		of those files that matched the regex search but not a particular run key
		stored in project.CellData
	'''

	# for the file names of the themis files
	themis_files = []
	# for the absolute paths of the themis files
	themis_file_paths = []
	# for the paths relative to proj file of the themis files
	themis_relative_file_paths = []

	# get the run keys from the proj dataframe
	runs = proj.CellData.index.get_level_values(level=1)

	# Finding the prefix used for themis files
	# So that changes to the relevant hermes function flow to this
	# function automatically
	themis_file_prefix = (
		re.match(
			r'([^_]*)_', # prefix = string without '_' before separator '_'
			mk_themis_file_name( # generate file name from first run key in proj.CellData
				**mk_cell_id_from_run_key(
						runs[0], proj
					)
				)
			)
			.group(1) # take first group
	)

	# walk through proj path
	# Look for all files that match the themis_file prefix and are ".pkl" files
	# add matched file names to themis_files
	# create absolute path and add to themis_file_paths
	for root, dirs, files in os.walk(proj._absolute_paths['path']):

		# all files that are themis files in this current directory
		matches =  fnmatch.filter(files, f'{themis_file_prefix}_*.pkl')

		themis_files.extend(matches)
		themis_file_paths.extend(
			[
				os.path.join(root, m)
				for m in matches
			]
		)
		themis_relative_file_paths.extend(
			[
				os.path.relpath( os.path.join(root, m), proj._absolute_paths['path'] )
				for m in matches
			]
			)


	# directory is dict: run_key -> absolute path | None
	# Initialise directory, where path remains None if not found in matched files
	run_file_paths = {
		k : None
		for k in runs
	}



	# Iterate through run_keys
	# generate appropriate file name for run key from hermes function
	# search for said file name in themis_files
	# wrap search in try,except and use index of successful search to retrieve absolute path from 
	# themis_file_paths (which is in the same order as themis_files)
	for r in runs:
		file_name = mk_themis_file_name(
			**mk_cell_id_from_run_key(r, proj)
		)

		try:
			file_idx = themis_files.index(file_name)

			run_file_paths[r] = dict(
				abs_path = pthl.Path(themis_file_paths.pop(file_idx)),
				rel_path = pthl.Path(themis_relative_file_paths.pop(file_idx))
				)

			# Not entirely happy about this popping business
			# For instance, how do duplictes work now ?!
			# Need to tidy up this code
			themis_files.pop(file_idx)

		# ValueError is for when list.index(obj) can't find object
		except ValueError:
			# leave value of None in directory dict
			pass

	if return_strays:
		return themis_file_paths, themis_relative_file_paths

	else:
		return run_file_paths


def mk_track_files_directory(proj, return_strays=False):
	'''
	Creates directory of track object files for the provided project object

	Parameters
	----
	return_strays : bool
		If true, don't return directory.
		Instead, return list of (absolute) file paths and relative file paths
		of those files that matched the regex search but not a particular run key
		stored in project.CellData
	'''

	# for the file names of the themis files
	track_files = []
	# for the absolute paths of the track files
	track_file_paths = []
	# for the paths relative to proj file of the track files
	track_relative_file_paths = []

	# get the run keys from the proj dataframe
	runs = proj.CellData.index.get_level_values(level=1)

	# Finding the prefix used for track files
	# So that changes to the relevant hermes function flow to this
	# function automatically

	# Pull out relevant data from cell that has some track data
	runs_with_track_data = proj.CellData.query('not track.isnull()')
	# check if any track data in project
	assert runs_with_track_data.index.size > 0, ('No track data assigned in project',
												 'no point in building directory')
	# pull out first experiment and track no
	track_exp, track_no = runs_with_track_data.iloc[0][['experiment', 'track']]

	track_file_prefix = (
		re.match(
			r'([^_]*)_', # prefix = string without '_' before separator '_'
			mk_track_file_name(experiment=track_exp,  track=track_no)
					)
			).group(1) # take first group

	# walk through directory
	# find and append matched files, paths and rel_paths
	for root, dirs, files in os.walk(proj._absolute_paths['path']):

		matches = fnmatch.filter(files, f'{track_file_prefix}_*.pkl')
		track_files.extend(matches)
		track_file_paths.extend(
			[
				os.path.join(root, m)
				for m in matches
			]
			)

		track_relative_file_paths.extend(
			[
				os.path.relpath( os.path.join(root, m), proj._absolute_paths['path'])
				for m in matches
			]
			)

	# initialise directory
	run_file_paths = {
		k: None
		for k in runs
	}

	# For storing paths for any given experiment x track combination so that all
	# runs of the same experiment and track number need only retrieve from the cache
	# after the first hit
	paths_cache = {}

	for r in runs:

		run_data = proj.CellData.xs(r, level=1)
		# For use in paths_cache as key
		# all track files are specific to an experiment and a track, thus this tuple as a key
		exp_track = tuple(run_data.loc[:, ['experiment', 'track']].values[0])


		# If no track data is assigned to the run in the project,
		# Then there is no way to create a putative file name for the track
		# Thus, without track data in the project, the file cannot appear in
		# the directory

		if pd.isnull(run_data.track.values[0]):
			continue

		else:

			# If in cache, take from cache
			if exp_track in paths_cache:
				run_file_paths[r] = paths_cache[exp_track]

			# If not in cache, generate file name and search candidate track files
			else:

				# create file names from run key
				file_name = mk_track_file_name(
					experiment=run_data.experiment.values[0], 
					track = run_data.track.values[0]
				)

				# Try to search for file_name in candidate track files
				try:
					file_idx = track_files.index(file_name)

					paths_object = dict(
						abs_path = pthl.Path(track_file_paths.pop(file_idx)),
						rel_path = pthl.Path(track_relative_file_paths.pop(file_idx))
						)

					track_files.pop(file_idx)

					# Append paths_object to list of paths and to cache using exp_track tuple key
					run_file_paths[r] = paths_object
					paths_cache[exp_track] = paths_object

				# if can't find file in list
				except ValueError:
					# leave value of None in directory dict
					pass

	if return_strays:
		return track_file_paths, track_relative_file_paths
	else:
		return run_file_paths


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


# > Cell ID and Key

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


def mk_cell_id_from_run_key(run_key, proj):

	'''
	Returns a dict similar to themis.CELL_ID that matches the run_key.

	This CELL_ID is derived from the CellData in proj (an athena obj)
	'''

	cell_id = (
			proj.CellData.xs(run_key, level=1)
					.loc[:, ['experiment', 'unit', 'cell', 'run']]
					.iloc[0] # to convert to series, by selecting row
					.to_dict()
		)

	return cell_id


def mk_cell_id_from_many_run_keys(run_keys, proj):
	'''
	Like mk_cell_id_from_run_key but for many run keys.
	Returns a dict of dicts : {run_key: {cell_id} }

	run_keys : iterable(strs) | 'all'
		if iterable, must be of strs/objects, each a run_key
		if 'all', then all runs are returned
	'''

	# slice(None) -> ':' for .loc[]
	run_keys = slice(None) if run_keys == 'all' else run_keys

	cell_idx = (
		proj.CellData.reset_index(level='cell_key', drop=True)
			.loc[run_keys, ['experiment', 'unit', 'cell', 'run']]
			.to_dict(orient='index')
		)

	return cell_idx



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
	object for track reconstruction
	'''
	
	# This can be made simpler
	# Take cut_plane distance, n_slides between lesions
	# Actual distance recorded in experiment and slide_width
	# automatically calculate shrinkage and theta (should be easier and more consistent)

	def __init__(self, experiment, track, lesion_depths_experimental,
				lesion_cut_plane_separation, lesion_n_slice_sag_separation, 
				px_scale, slice_width=50):
		'''
		initialise track object
		include information on track location and dimensions

		parameters
		----

		lesion_depths_experimental : iterable
			list or iterable of the depths at which lesions
			were made, as recorded at experiment time
			by a microdrive etc.

		lesion_cut_plane_separation : float (pxls)
			separation between the lesions, in the histology
			cutting plane, measured in pxls

		lesion_n_slice_sag_separation : float (slices)
			number of histological slices separating the two lesions
			in the plane orthogonal to the cutting plane

		pxl_scale : float
			how many pixels per micrometer (um) of primary 
			histology imagery.
			eg, 0.3 means 0.3px per um or 3.3um per px

		slice_width : float (default: 50)
			width of individual slices of the histology
		'''

	#     print(inspect.getargvalues(inspect.currentframe()))
	#     print(locals())

		input_vars = {}
		input_arg_info = inspect.getargvalues(inspect.currentframe())

		# take only input arguments prescribed by function definition
		input_vars.update(
			{a : input_arg_info.locals[a] 
			 for a in input_arg_info.args 
			 if a != 'self'}
		)
		
		# self.__dict__.update(input_vars)

		self.experiment = experiment
		self.track = track

		self.lesion_depths_experimental = sorted(lesion_depths_experimental)
		self.px_scale = px_scale
		self.slice_width = slice_width

		self.lesion_cut_plane_separation = lesion_cut_plane_separation / self.px_scale

		self.lesion_n_slice_sag_separation = lesion_n_slice_sag_separation
		self.lesion_sag_separation = lesion_n_slice_sag_separation * slice_width

		self.angle = self._atan(self.lesion_sag_separation , self.lesion_cut_plane_separation)

		self.shrinkage = 1 - (
				np.hypot(self.lesion_sag_separation, self.lesion_cut_plane_separation) /
				np.diff(self.lesion_depths_experimental)[0]
			)


		# for saving any save paths
		self._absolute_paths = {}
		

	def setUnits(self, depths, unit_no):

		'''
		add unit information to track object
		unit information will be immediately used to calculate 
		various figures regarding the estimated location of the units.
		this is done using self._calcDepthValues()

		parameters
		----

		depths : list/iterable (float)
			list of floats representing the depths of each
			unit of this track.
			
			as with self.__init__(), depths represent the depth
			as recorded at experiment time, by a microdrive etc.

		unit_no : list/iterable (int)
			list of ints representing the numbers assigned to
			the units of this track.
			must be in the same order as depths :
			[123, 647], [2, 3] -> unit 3 depth = 647
		'''
		
		unit_depths = np.array(depths)
		
		# this sorting defines the sorting of all subsequent calculations
		sort_idx = np.argsort(unit_depths)
		
		self.unit_depths = unit_depths[sort_idx]
		self.unit_nos = np.array(unit_no)[sort_idx]

		self._calcDepthValues()
		

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

	def _atan(self, opp, adj):
		return np.degrees(np.arctan(opp / adj))
		
	def _calcDepthValues(self):

		'''
		from track information and unit depths, calculate the
		location of the units relative to the cutting plane and
		last lesion

		creates following attributes:
		self.reference_depth
		self.depth_values
		self.df_depth_values
		'''

		depth_values = {}
		# take the deepest lesion as the reference point
		reference_depth = self.lesion_depths_experimental[-1]

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
		define layer borders, heading into cortex (from layer 1)
		
		parameters
		---
		borders : list 
			define when the beginning of the relevant border is.
			values are presumed to be in pixels of the primary histology
			images, sharing the same px_scale as provided in self.__init__
			
			relative to the last lesion location, and in cutting plane.
			
			must be in same order as layers - user discretion.
			and, in numerical order (ie smallest depth first)
			
		layers : list
			provide labels for the relevant borders, in ints
			must be in same order as borders.
			
			will be converted to strings for labelling purposes
		'''

		# check that borders are floats
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

				# calculate layers each unit belongs in (where possible)

		unit_layer_labels = [None] * len(self.unit_depths)

		for ui, ud in enumerate( self.depth_values['rel_depth_cut_plane_pxl'] ):
			for li, lb in enumerate(self.layer_borders):

				if ud > lb:
					unit_layer_labels[ui] = self.layer_labels[li]

					self.unit_layer_labels = unit_layer_labels

		# calculate relative position of unit within layer (where possible)

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

		# kinda given up on special path changing abilities and enforcement
		# this object gets saved in the current, path, that simple
		directory = pthl.Path('.')
		file_name = mk_track_file_name(track_obj=self)

		save_path = directory / file_name

		self._absolute_paths['save_path'] = save_path.absolute()
		self.save_path = save_path

		temp_path = save_path.with_suffix(save_path.suffix + '.tmp')

		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(save_path)

		print(f'saved track object as {save_path}')

