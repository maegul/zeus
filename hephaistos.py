'''
Library for generating spike times from raw data.

Intention is for spike times to be analysed by Zeus

Hephaistos -> god of craft and metalwork ... this is the hard work!
'''

import pathlib as pthl
from pprint import pprint

from pyth.zeus import spk2_mat_read as spkmr
import neo
import quantities as qnt

import numpy as np
import pandas as pd



# Decorator for inserting doc strings

def doc_string_add(inserted_func):
	'''
	Inserts doc string of inserted func into designated format `{0} locations
	If multiple functions provided, PRESUMES location is assigned in target doc string
	If not so provided, failur to insert will pass silently (as str.format() raises no errors)


	'''

	if not isinstance(inserted_func, list) and callable(inserted_func) :
		inserted_func = [inserted_func]
	for insfunc in inserted_func:
		assert callable(insfunc), f'Object {insfunc} is not a function, cannot insert doc string'

	def wrapper(func):

			func.__doc__ = func.__doc__.format(*[insfunc.__doc__ for insfunc in inserted_func])

			return func

	return wrapper



# Functions or classes?

# data mangemeent ... raw v output data?  ... need to keep in separate directories

class Hephaistos:

	'''
	'''

	def __init__(self, data_path='', output_path='.'):

		'''
			Class for managing raw data -> spike data for zeus

			Aim is to store key information (origin of raw data and actual spike data)

			NOTE - design of spk2_mat_read is to only call large data, not store internally

			Parameters
			__________

			data_path : str
				path to data file
			output_path : str
				path to location where
		'''

		self.Data_path = pthl.Path(data_path)
		self.Output_path = pthl.Path(output_path)

		assert self.Data_path.is_file(), 'data path not a file, please provide file with data'

		assert self.Output_path.is_dir(), 'output path is not a directory, need directory to store files in'

		self.RawBinChans = rawBinChans()

		# How manage paths?
			# Want to be able to easily pull out containing folders

	@doc_string_add(spkmr.spk2MatRead.__init__)
	def readSpkMat(self, **kwargs):

		'''

		Spike 2 Mat Read function args
		_______

		{0}

		'''

		self.Dat = spkmr.spk2MatRead(self.Data_path, **kwargs )


	def mkRawBin(self, channel = 'Ch1', filename = None, backend = 'neo', overwrite = False,
					dtype = 'float64', nb_channel = 1, sampling_rate = 2e4,
					comment = ''):
		'''
		Save raw binary file for TDC to data path

		One data trace at a time (TDC weird with multi channel)
		Default is Ch1 (open), but can specify

		'''



		backends = [
					'neo',
					'np'
					]

		# Make sure channel provided is available
		assert channel in self.Dat._traceChannelNames, 'channel name not in data, check available channels with self.Dat._traceChannelNames'


		assert backend in backends, f'not valid backed, select from {backends}'


		if not filename:
			# relying on pathlib attribute here (provides filename from end of path without ext)
			filename = self.Data_path.stem + '.raw'

		elif filename[-4:] != '.raw':
			# Add .raw extension to filename
			filename = filename + '.raw'

		# Path for writing binary file
		self.rawDataFilePath = self.Data_path.parent / filename


		args = locals()

		# Take out unwanted or boilerplate variables
		for bp_arg in ['self']:
			args.pop(bp_arg)

		args.update( dict(rawDataFilePath = self.rawDataFilePath))



		# Add channel information to catalogue of channels (also check if already been written)
		chan_check = self.RawBinChans.listChannels(chan_name=channel)
		assert not chan_check['processed'] or overwrite, f"{channel} already processed, written: {chan_check['written']}! overwrite: {overwrite}?"
		self.RawBinChans.addChannel(channel, args)

		# Should only need channel, backend and filename, rest from Dat
		# channel should be one of teh above channels!  Checked above, but just sayin
		signal = self.Dat.__getattribute__(channel).GetTraceData()
		signal = signal.astype(dtype)

		if backend == 'neo':

			neoBinDat = neo.io.RawBinarySignalIO(filename = self.rawDataFilePath, 
							dtype = dtype, nb_channel = nb_channel, 
							sampling_rate = sampling_rate*qnt.Hz)


			anasig = neo.core.AnalogSignal(signal = signal, units='V', 
											sampling_rate = sampling_rate*qnt.Hz, dtype = dtype)

			seg = neo.core.Segment()
			seg.analogsignals.append(anasig)

			neoBinDat.write_segment(seg)

			# Add parameter for whether binary file successfully written
			self.RawBinChans.__dict__[channel].written = True


		if backend == 'np':

			signal.tofile(self.rawDataFilePath)


			self.RawBinChans.__dict__[channel].written = True






class rawBinChans:
	'''
	To create objects that store information on channels converted to raw binary files
	'''


	class chan:
		def __init__(self, params):
			self.written = False
			self.params = params

			self.__dict__.update(params)					

		def pparams(self):
			pprint(self.params)


	def addChannel(self, chan_name, params):

		self.__dict__.update( { chan_name : self.chan(params) })

	def updateChannel(self, chan_name, params):

		assert self.listChannels(chan_name=chan_name), f'{chan_name} not available'

		assert isinstance(params, dict), 'params provided are not dictionary'

		getattr(self, chan_name).__dict__.update(params)


	def listChannels(self, chan_name=None):

		# Produce list of attributes of object that are of type self.chan (ie, are channels)
		existing_chans = [att for att in dir(self) if isinstance(getattr(self, att), self.chan)]

		if not chan_name:
			return existing_chans

		elif isinstance(chan_name, str):
			chan_processed = chan_name in existing_chans

			chan_written = self.__dict__[chan_name].written if chan_processed else False

			return dict(processed=chan_processed, written=chan_written)




