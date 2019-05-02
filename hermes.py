import pandas as pd

from time import time
from fuzzywuzzy import fuzz, process

# from collections import OrderedDict as OD
import inspect

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
def initDataSet(key, cell_data, stim_data):

	cell_data.update(stim_data)

	dataSet = pd.DataFrame(data=cell_data, index = [key])
	dataSet.index.name = 'key'

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



def appendDataSet(key, cell_data, stim_data, dataSet, force=False):


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

	cell_data = pd.Series(data = cell_data, name = key)

	# Relying on integrity check to prevent data being added twice
	dataSet = dataSet.append(cell_data, verify_integrity=True)

	
	return dataSet	



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


