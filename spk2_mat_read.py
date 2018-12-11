# interval == sampling interval time
# length == number of samples
# start == ? ... perhaps a start time to synchronise other channels with different samples rates ...?
from pprint import pprint

import re
import h5py
import numpy as np




class spk2MatRead:
	'''
	Read mat (HDF5) data formats exported from spike 2

	PRESUMES - data file has only trace and marker data, perhaps two traces

	'''

	def __init__(self, file_name, trace_channels=['Ch1', 'Ch2'], baseRe='(\w+)Ch\d+'):

		self.Dat = h5py.File(file_name)

		self.baseRe = baseRe
		self._baseNameRe = re.compile(r''.join(baseRe))

		self._markChannelRe = re.compile('_Ch32$')

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

		self._trace_chan_args = trace_channels
		self.TraceChannels = []

		# Keys exclusive to trace channels (known through manual inspection)
		# Using sets for checking for inclusion
		trace_key_check = {'start'} 
		spk_key_check = {'times'}

		for tc in trace_channels:

			chan_name = self._baseName + tc
			assert chan_name in self.Channels, f'Channel "{self._baseName}" + "{tc}" is not available'

			chan = self.Dat[chan_name]
			self.TraceChannels.append(chan_name) # for iterating, if it's useful


			if trace_key_check.issubset( set( chan.keys())):
				print(f'{tc} treated as trace')
				self.__dict__.update({tc: trace(chan)})

			elif spk_key_check.issubset( set( chan.keys())):
				print(f'{tc} treated as spike template channel')
				self.__dict__.update({tc: spikeTemp(chan)})







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

	def __init__(self, chan):

		self.Dat = chan

		self.DatKeys = self.listDataSets(prnt=False)



		trace_param_keys = {'scale': 'scale',# not sure what this is (magnitude scaling?)
							'start': 'start',# time point of first data point for syncing?
							'length': 'size', # number of data points
							'interval': 'samp_time', # time between samples (samp freq ~ 1/interval)
							'offset': 'offset' # Don't know what this is
							}

		self.GetTraceData = lambda : chan['values'].value[0]


		# Might be nice to in future turn all params to attributes ... need recursive dict / node / plane object update ... ?
		self.params = dict()
		for k, name in trace_param_keys.items():

			try:
				values = chan[k].value

				if values.size == 1:
					values = values[0]

				self.params.update({name: values })

			except:

				print(f'\n{k} is not an attribute of {chan.name} ... not the channel you were expecting???')



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
			return list(self.Dat.keys())





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
