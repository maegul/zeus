
from . import hephaistos as heph, themis

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


def copyCID(cell, prevCell, **updateKwArgs):

	cid = prevCell.CELL_ID.copy()

	cid.update(**updateKwArgs)

	cell._cell_ID(**cid)

	cell._cell_id_print()


def copySP(cell, prevCell, **updateKwArgs):

	sp = prevCell.STIM_PARAMS.copy()

	sp.update(**updateKwArgs)

	cell._stim_params(**sp)

	cell._stim_params_print()





