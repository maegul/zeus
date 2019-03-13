
from . import hephaistos as heph, themis, athena

np = themis.np
plt = themis.plt


def quickHeph(path = None):
	'''
	Quickly produces a hephaistos unit 
	Much of the tdc work is still to be done by the user

	returns a hephaistos object
	'''

	unit = heph.Hephaistos(data_path = path)
	unit.readBasicSpkMat()
	unit.mkRawBin()
	unit.tdcInit()

	return unit


def quickThemis(data_path = None, cell_no = None, project = None):

	cell = themis.Themis(data_path = data_path, cell = cell_no)

	cell._assign_project(project)

	return cell


def copyCID(currentCell, prevCell, **updateKwArgs):

	cid = prevCell.CELL_ID.copy()

	cid.update(**updateKwArgs)

	currentCell._cell_ID(**cid)

	currentCell._cell_id_print()


def copySP(cell, prevCell, force=False, **updateKwArgs):

	sp = prevCell.STIM_PARAMS.copy()

	sp.update(**updateKwArgs)

	cell._stim_params(**sp, force=force)

	cell._stim_params_print()


def quickAthena(path = None):
	proj = athena.load(path)

	# Number of runs for each cell
	summary = proj.CellData.loc[:,['experiment', 'unit', 'cell', 'run']].groupby(by=['experiment', 'unit', 'cell']).count()

	# Number of unique indices from the multiindex of the grouping (experiment, unit, cell) ... ie, how many unique cells
	n_cells = summary.size

	n_animals = proj.CellData.experiment.unique().size

	print(f"{proj.PROJ_ID['name']}:\t{proj.PROJ_ID['description']}\n")
	print(f"{n_animals} Animals   |   {n_cells} Cells\n")
	print(summary)


	return proj





