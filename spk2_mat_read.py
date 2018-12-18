'''
Objects for reading HDF5 (matlab >7.3) data files containing spike2 raw data

'''
# interval == sampling interval time
# length == number of samples
# start == ? ... perhaps a start time to synchronise other channels with different samples rates ...?
from pprint import pprint

import re
import h5py
import numpy as np


class basicMatRead:
	'''
		Read basic manually generated matlab file format.

		Generated using the SONVEIWER matlab libraries provided by CED



	'''

	def __init__(self, file_name, trace_channels=['Ch1'], marker_channel = 'Ch32',
					base = 'output_struct'):

		'''
		Parameters
		__________

		file_name : str
			file name or path to file containing HDF5 (.mat) data

		trace_channels : iterable
			Default - ['Ch1']
			
			Channel names that are to be extracted for easy access

			Names must match the keys used in the HDF5 file.
			PRESUMING data originates from spike2, CHannel names should simply 'Ch#' with
			"#" being the number of the channel in spike2.

			Ch32 is PRESUMED to be the marker channel

		baseRE : str regex
			RegEx for identifying the base key of the file. 
			This is orindarily derived from the name of the spike2 file.
			First group is used as the base

		markerRe : str regex
			As with baseRe, but for identifying the marker channel


		Attributes
		__________

		Channels : list
			Keys of all the channels in the data file

		channel attrs : HDF5 objects
			For each of the Trace Channels provided in the kwarg, there is
			an attribute of the same name of either the class Trace or spikeTemp
			This is for convenient access

		MarkTimes : array (1d)
			times of the markers from the marker channel

		MarkCodes : array (1d)
			Codes for each Marker Time.
			Should encode stimulus type ... needs to be derived from lab data

		TraceChannels : list
			All Trace channels encoded as Trace types in same order as 
			trace_channels (kwarg) (where possible) and _traceChannelKeys

		SpikeTempChannels : list
			As TraceChannels above


		Methods
		______

		listChanData
			Pretty Prints all channels and their attributes





		'''

		self.Dat = h5py.File(file_name)
		self._baseName = base
		self.DatBase = self.Dat[self._baseName]

		self.Channels = self.listChannels(prnt=False)


		self._markerChannelArg = marker_channel
		self._trace_chan_args = trace_channels

		self._marker_time_key = 'times'
		self._marker_code_key = 'codes'

		self._markChannel = self.DatBase[self._markerChannelArg]
		self.MarkTimes = self._markChannel[self._marker_time_key].value.flatten()
		self.MarkCodes = self._markChannel[self._marker_code_key].value.flatten()




		self.TraceChannels = []
		self.SpikeTempChannels = []

		self._traceChannelNames = []
		self._spikeTempChannelNames = []





		# Adding specified channels as attributes and filtering by type
		for tc in trace_channels:

			assert tc in self.Channels, f'Channel {tc} is not available'

			chan = self.DatBase[tc]



			self.__dict__.update({tc: trace(chan, tc)})

			self._traceChannelNames.append(tc)
			self.TraceChannels.append(getattr(self, tc))








	def listChannels(self, prnt=True):

		# With basic format, everything is under basename
		ks = list(self.Dat[self._baseName].keys())
		if prnt:
			pprint(ks)
		else:
			return ks

	def listChanData(self, prnt=True):
		self.ChanData = dict()

		for ch in self.Channels:
			self.ChanData.update({f'{ch}' : list(self.Dat[ch].keys())})

		if prnt:
			pprint(self.ChanData)
		else:

			return self.ChanData



class spk2MatRead:
	'''
	Read mat (HDF5) data formats exported from spike 2

	Object contains all relevant data, or has methods to accessing (larger) data
	from the HDF5 file.

	PRESUMES - data file has only trace and marker data.

	Currently, taking out templated spike data is limited.

	'''

	def __init__(self, file_name, trace_channels=['Ch1'], baseRe='(\w+)Ch\d+',
					markerRe = '_Ch32$'):

		'''
		Parameters
		__________

		file_name : str
			file name or path to file containing HDF5 (.mat) data

		trace_channels : iterable
			Default - ['Ch1']
			
			Channel names that are to be extracted for easy access

			Names must match the keys used in the HDF5 file.
			PRESUMING data originates from spike2, CHannel names should simply 'Ch#' with
			"#" being the number of the channel in spike2.

			Ch32 is PRESUMED to be the marker channel

		baseRE : str regex
			RegEx for identifying the base key of the file. 
			This is orindarily derived from the name of the spike2 file.
			First group is used as the base

		markerRe : str regex
			As with baseRe, but for identifying the marker channel


		Attributes
		__________

		Channels : list
			Keys of all the channels in the data file

		channel attrs : HDF5 objects
			For each of the Trace Channels provided in the kwarg, there is
			an attribute of the same name of either the class Trace or spikeTemp
			This is for convenient access

		MarkTimes : array (1d)
			times of the markers from the marker channel

		MarkCodes : array (1d)
			Codes for each Marker Time.
			Should encode stimulus type ... needs to be derived from lab data

		TraceChannels : list
			All Trace channels encoded as Trace types in same order as 
			trace_channels (kwarg) (where possible) and _traceChannelKeys

		SpikeTempChannels : list
			As TraceChannels above


		Methods
		______

		listChanData
			Pretty Prints all channels and their attributes





		'''

		self.Dat = h5py.File(file_name)

		self.baseRe = baseRe
		self._baseNameRe = re.compile(r''.join(baseRe))

		self.markerRe = markerRe
		self._markChannelRe = re.compile(r''.join(markerRe))

		self.Channels = self.listChannels(prnt=False)

		# using first channel as biasis (presumption)
		# using first group of regex as base ... error will be thrown if non existent
		self._baseName = self._baseNameRe.search(self.Channels[0]).group(1)

		# Inefficient / lazy use of regex here ... sighs
		for ch in self.Channels:
			if self._markChannelRe.search(ch):

				print(f'Channel {ch} treated as markers channel')

				self._markChannel = self.Dat[ch]

				# Take first element, as first deminsion is redundant
				self.MarkTimes = self._markChannel['times'][0]
				self.MarkCodes = self._markChannel['codes'][0]

				break

		else:
			print(f'No marker channel found.  Revise markerRe {markerRe}')


		self._trace_chan_args = trace_channels

		self.TraceChannels = []
		self.SpikeTempChannels = []

		self._traceChannelKeys = []
		self._spikeTempChannelKeys = []

		self._traceChannelNames = []
		self._spikeTempChannelNames = []

		# Keys EXCLUSIVE to types of channels (known through manual inspection)
		# Using sets for checking for inclusion
		trace_key_check = {'start'} 
		spk_key_check = {'times'}


		# Adding specified channels as attributes and filtering by type
		for tc in trace_channels:

			chan_key = self._baseName + tc
			assert chan_key in self.Channels, f'Channel "{self._baseName}" + "{tc}" is not available'

			chan = self.Dat[chan_key]


			if trace_key_check.issubset( set( chan.keys())):
				print(f'{tc} treated as trace')

				self.__dict__.update({tc: trace(chan, tc)})

				self._traceChannelKeys.append(chan_key) # for iterating, if it's useful
				self._traceChannelNames.append(tc)
				self.TraceChannels.append(getattr(self, tc))

			elif spk_key_check.issubset( set( chan.keys())):
				print(f'{tc} treated as spike template channel')

				self.__dict__.update({tc: spikeTemp(chan)})

				self._spikeTempChannelKeys.append(chan_key)
				self._spikeTempChannelNames.append(tc)
				self.SpikeTempChannels.append(getattr(self, tc)) # for iterating, if it's useful







	def listChannels(self, prnt=True):

		ks = list(self.Dat.keys())
		if prnt:
			pprint(ks)
		else:
			return list(self.Dat.keys())

	def listChanData(self, prnt=True):
		self.ChanData = dict()

		for ch in self.Channels:
			self.ChanData.update({f'{ch}' : list(self.Dat[ch].keys())})

		if prnt:
			pprint(self.ChanData)
		else:

			return self.ChanData



class trace:

	def __init__(self, chan, chan_name):

		self.Dat = chan
		self.Chan_Name = chan_name

		self.DatKeys = self.listDataSets(prnt=False)



		trace_param_keys = {'scale': 'scale',# not sure what this is (magnitude scaling?)
							'start': 'start',# time point of first data point for syncing?
							'length': 'size', # number of data points
							'interval': 'samp_time', # time between samples (samp freq ~ 1/interval)
							'offset': 'offset' # Don't know what this is
							}

		# Aim is to return the whole flattened array.  Often slow
		self.GetTraceData = lambda : chan['values'].value[0]

		# For slicing etc
		self.TraceData = chan['values']


		# Might be nice to in future turn all params to attributes ... need recursive dict / node / plane object update ... ?
		self.params = dict()

		for k, name in trace_param_keys.items():

			try:
				values = chan[k].value

				if values.size == 1:
					values = values[0]

				self.params.update({name: values })

			except:

				print(f'\n{k} is not an attribute of {self.Chan_Name} ... not the channel you were expecting???')


		# Create time array for each data point in trace data

		try:
			step = self.params['samp_time'][0]
			self.SampleFreq = 1 / step

			start = self.params['start']
			size = self.params['size']
			stop = start + (step * (size-1)) # first sample is at start, so one fewer steps to add, as endpoint is True and is counted

			self.Times = np.linspace(start=start, num=size, stop=stop, endpoint=True)

		except:
			# Can manage exceptions here more specifically
			print('Making time series failed ... appropriate params not available?')




	def listDataSets(self, prnt=True):

		ks = list(self.Dat.keys())
		if prnt:
			pprint(ks)
		else:
			return ks





class spikeTemp:

	def __init__(self, chan):

		self.Dat = chan

		self.DatKeys = self.listDataSets(prnt=False)

		for k in self.DatKeys:
			self.__dict__.update({ f'get_{k}': lambda: chan[k].value})



	def listDataSets(self, prnt=True):

		ks = list(self.Dat.keys())
		if prnt:
			pprint(ks)
		else:
			return list(self.Dat.keys())
