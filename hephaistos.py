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
from scipy.io.matlab import savemat
from scipy.stats import sem

from matplotlib import pyplot as plt

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

import tridesclous as tdc

from tridesclous import DataIO, CatalogueConstructor, Peeler

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
				path to location where spike data files will be placed
		'''

		self.Data_path = pthl.Path(data_path).absolute()
		self.Output_path = pthl.Path(output_path).absolute()

		assert self.Data_path.is_file(), 'data path not a file, please provide file with data'

		assert self.Output_path.is_dir(), 'output path is not a directory, need directory to store files in'

		# Add folder for all tdc data of the same name as the data file
		self._processing_data_path = self.Data_path.parent / self.Data_path.stem # stem removes suffix of file, useful to create new folder 
		if not self._processing_data_path.is_dir():
			self._processing_data_path.mkdir()

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


	@doc_string_add(spkmr.basicMatRead.__init__)
	def readBasicSpkMat(self, **kwargs):

		'''

		Basic Spike 2 Mat read function args
		_______

		{0}
		'''

		self.Dat = spkmr.basicMatRead(self.Data_path, **kwargs)


	def writeWaveClusMat(self, channel = 'Ch1'):

		mat_file_path = self._processing_data_path / ('wc_' + self.Data_path.stem)

		savemat(mat_file_path, dict(data = self.Dat.Ch1.GetTraceData(),
									sr = self.Dat.Ch1.SampleFreq),
				appendmat = True
				)




	def mkRawBin(self, channel = 'Ch1', filename = None, backend = 'neo', overwrite = False,
					dtype = None, nb_channel = 1, comment = ''):
		'''
		Save raw binary file for TDC to data path

		filename is inherited by default as file name of data path

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
			filename = f'{self.Data_path.stem}_{channel}.raw'

		elif filename[-4:] != '.raw':
			# Add .raw extension to filename
			filename = filename + '.raw'

		# Path for writing binary file
		self.rawDataFilePath = self._processing_data_path / filename

		sampling_rate = self.Dat.__dict__[channel].SampleFreq


		args = locals()

		# Take out unwanted or boilerplate variables
		for bp_arg in ['self']:
			args.pop(bp_arg)

		args.update( dict(rawDataFilePath = self.rawDataFilePath))



		# Check if channel has been processed before
		# Checking if processed and/or written before
		chan_check = self.RawBinChans.listChannels(chan_name=channel)
		# Either not processed before, or, processed before but overwrite is True
		assert not chan_check['processed'] or overwrite, f"{channel} already processed, written: {chan_check['written']}! overwrite: {overwrite}?"

		# Binary has been written before, but want to overwrite, so remove binary file
		# Just in case neo or numpy is unhappy with preexisting file
		if (self.rawDataFilePath.is_file()) and overwrite:
			# strangely, pathlib removes files with unlink
			self.rawDataFilePath.unlink()


		# Should only need channel, backend and filename, rest from Dat
		# channel should be one of teh above channels!  Checked above, but just sayin
		signal = self.Dat.__getattribute__(channel).GetTraceData()

		if dtype is not None:
			signal = signal.astype(dtype)
		else:
			dtype = signal.dtype.name

		args.update( dict(dtype = dtype))

		self.RawBinChans.addChannel(channel, args)

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
			self.RawBinChans.updateChannel(channel, {'written': True})


		if backend == 'np':

			signal.tofile(self.rawDataFilePath)

			self.RawBinChans.updateChannel(channel, {'written': True})


	def quickView(self, channel = 'Ch1', limit = 10000, viewLimit = 2):

		'''
		viewlimit is how many times greater the samples in full view are than limit

		limit is how many samples for downsampled view
		'''


		assert channel in self.Dat._traceChannelNames, 'channel name not in data, check available channels with self.Dat._traceChannelNames'

		# Bad to simply store the signal in memory ... could be replaced by an HDF5 read function
		signal = self.Dat.__getattribute__(channel).GetTraceData()
		# Sample rate
		sr = self.Dat.__getattribute__(channel).SampleFreq

		viewLimit *= limit

		# Number of samples skipped per down-sampled 'sample'
		DSIdx = np.floor(signal.size / limit).astype('int')
		testDataDS = signal[::DSIdx] # Downsampled dataset

		viewTime = np.arange(0, viewLimit/sr, 1/sr)
		dsTime = np.arange(testDataDS.size) * 1/sr * DSIdx

		win = pg.GraphicsWindow(title="Trace")
		win.resize(1000,600)
		win.setWindowTitle(f'{channel} Trace')

		# Enable antialiasing for prettier plots
		pg.setConfigOptions(antialias=True)

		p1 = win.addPlot(title='Whole Trace (downsampled)')
		p1.plot(y = testDataDS, x = dsTime)

		lr = pg.LinearRegionItem([0, viewLimit//DSIdx])
		lr.setZValue(-10)
		p1.addItem(lr)


		win.nextRow()
		p2 = win.addPlot(title="Trace (sliced)")

		p2plt = p2.plot(y=signal[:viewLimit], x = viewTime, pen=(200,90,170))


		def updatePlot():
			# getRegion returns in terms of x axis values
			# Convert real time to number of samples
		    idxs = [int(idx*sr) for idx in lr.getRegion()]
		    # 
		    viewTime = np.arange(0, (idxs[1]-idxs[0])/sr, 1/sr)

		    p2plt.setData(y = signal[idxs[0]:idxs[1]], x = viewTime)

		lr.sigRegionChanged.connect(updatePlot)

		# def updateRegion():
		    # lr.setRegion(p2.getViewBox().viewRange()[0])
		# p2.sigXRangeChanged.connect(updateRegion)

		updatePlot()

		# Do I need to return an object to cancel things?
		QtGui.QApplication.instance().exec_()
		



	def tdcInit(self, channel='Ch1'):
		'''
		DataIO directory and DataIO creation
		'''

		dirname = self._processing_data_path

		self.TDCDataIO = DataIO(dirname=dirname)

		# 'datasource_type' is a key created by initialising datasource in the info dict
		# If datasource_type is in info, than data source has been set already in the data_path for processing
		if 'datasource_type' not in self.TDCDataIO.info:
			assert channel in self.RawBinChans.listChannels(), 'Please provide appropriate channel'

			self.TDC_channel = channel

			dataSource = self.RawBinChans.__dict__[channel]

			assert dataSource.written, 'Channel not yet written, use self.mkRawBin()'


			dataSourceType = 'RawData'
			dataSourcePath = self.rawDataFilePath
			dataSourceDType = dataSource.dtype
			dataSourceSampleRate = dataSource.sampling_rate

			self.TDCDataIO.set_data_source(type=dataSourceType, filenames=[dataSourcePath],
				dtype=dataSourceDType, sample_rate=dataSourceSampleRate, total_channel=1)

			# Encode channel used in the TDC info.json for when using preprocessed data
			# ie, when 'datasource_type' is in info
			self.TDCDataIO.info['spkmr_channel'] = channel
			self.TDCDataIO.flush_info()

		else:
			previous_channel = self.TDCDataIO.info['spkmr_channel']

			assert previous_channel == channel, f'Channel arg ({channel}) is not same as previously encoded TDC channel ({previous_channel})'
			
			# recreating variable used above for future reference below
			self.TDC_channel = previous_channel

			print(f'DataSource already set. {previous_channel} used. \nProcess data again, or, Alter path to start again.\n')
			print(self.TDCDataIO)


	@doc_string_add(CatalogueConstructor.set_preprocessor_params)
	def tdcPreProcess(self, peak_sign='-', relative_threshold=5, highpass_freq=300, 
		lowpass_freq=5000, noiseEstimateDuration = 10, catalogue_duration=400,
		**kwargs):
		'''
		filter etc with signalprocessor
		Importantly, sets duration of signal to be used in the catalogue

		Relative threshold is in MAD (median deviations)

		noiseEstimateDuration is length of signal for estimating noise.
		Short durations should be fine, unless noise is not stable

		kwargs passed to set_preprocessor_params

		CatalogueConstructor.set_preprocessor_params Doc String
		-----
		{0}
		'''

		self.TDCCatConstructor = CatalogueConstructor(dataio = self.TDCDataIO)

		self.TDCCatConstructor.set_preprocessor_params(peak_sign=peak_sign, 
			relative_threshold=relative_threshold, highpass_freq=highpass_freq,
			lowpass_freq=lowpass_freq, **kwargs)

		self.TDCCatConstructor.estimate_signals_noise(duration=noiseEstimateDuration)

		self.TDCCatConstructor.run_signalprocessor(duration=catalogue_duration)

		pprint(self.TDCCatConstructor)
		pprint(self.TDCCatConstructor.info)



	def tdcSetUpCatalogue(self, n_left=-30, n_right=30, align_waveform=True, 
		alien_value_threshold = 100):
		'''
		Setup the catalogue and run both PCA and clustering

		n_left and n_right expressed in milliseconds, converted to samples in code

		alien_value_threshold is expressed MAD (I think)
		'''

		assert n_left < 0 and n_right > 0, 'n_left must be negative, n_right positive'

		self.TDC_nLeft = n_left
		self.TDC_nRight = n_right

		self.TDCCatConstructor.extract_some_waveforms(n_left = n_left, n_right = n_right,
			align_waveform=align_waveform)

		self.TDCCatConstructor.clean_waveforms(alien_value_threshold=alien_value_threshold)



	def tdcCatPCAClust(self, n_pca_components = 7, pca_method = 'global_pca',
		clust_method = 'kmeans', n_clusters=5):

		self.TDCCatConstructor.extract_some_features(method=pca_method, 
			n_components=n_pca_components)

		# So I am assuming that the pca method will be used, and none of the others!
		# Could add layer of abstraction for other methods ... later
		self.TDC_pca = self.TDCCatConstructor.projector.pca

		self.TDC_pca_expl_var = self.TDC_pca.explained_variance_
		self.TDC_pca_expl_var_rat = self.TDC_pca.explained_variance_ratio_
		self.TDC_pca_expl_var_sum = np.sum(self.TDC_pca_expl_var_rat)

		self.TDCCatConstructor.find_clusters(method = clust_method, n_clusters = n_clusters)

		self.TDC_clusters = self.TDCCatConstructor.clusters

		print(self.TDCCatConstructor)



	def tdcFindBetterBounds(self, mad_thresh = 1.01):
		'''
		attempts to find better bounds for waveforms

		very thin wrapper around find_good_limits()
		'''
		return self.TDCCatConstructor.find_good_limits(mad_threshold=mad_thresh)



	def tdcViewPCA(self, figsize = (10, 8), n_random_waveforms = 21):
		'''
		View, plot and analyse principal components coming out of the catalogue

		Currently using matplotlib for plotting
		'''

		plt.figure(figsize=figsize)
		plt.subplot(211)
		for i in range(self.TDC_pca.components_.shape[0]):
			factor = 2 * np.sqrt(self.TDC_pca_expl_var[i])

			plt.plot(np.abs(self.TDC_pca.components_[i,:])*factor,
				label = f'{i} ({np.round(self.TDC_pca_expl_var_rat[i], 3)})'
				)
		# plot center of spike	
		plt.axvline(x = -1*self.TDC_nLeft, dashes=[3,2], color='0.65')
		plt.legend()
		plt.title(f'Sum Var Explained: {self.TDC_pca_expl_var_sum}')

		# Plt some random waveforms
		plt.subplot(212)

		random_wf_idxs = np.random.randint(0, 
			self.TDCCatConstructor.some_waveforms.shape[0], n_random_waveforms
			)

		for idx in random_wf_idxs:
			plt.plot(self.TDCCatConstructor.some_waveforms[idx, ...].flatten(), color='steelBlue')
		plt.axvline(x = -1*self.TDC_nLeft, dashes=[3,2], color='0.65')






	def tdcOpenCatalogue(self):
		'''
		Open catalogue window

		Quit python app to close and return to interpreter environment
		'''

		app = pg.mkQApp()
		win = tdc.CatalogueWindow(self.TDCCatConstructor)
		win.show()
		app.exec_()

	def tdcPeel(self):
		'''
		Set up peeler, with catalogue, and run
		'''

		self.TDCCatConstructor.make_catalogue_for_peeler()

		self.TDCCatalogue = self.TDCDataIO.load_catalogue()

		self.TDCPeeler = Peeler(self.TDCDataIO)
		self.TDCPeeler.change_params(catalogue=self.TDCCatalogue)

		self.TDCPeeler.run()



	def tdcOpenPeeler(self):
		'''
		Open peeler window
		'''

		app = pg.mkQApp()
		win = tdc.PeelerWindow(dataio=self.TDCDataIO, catalogue=self.TDCCatalogue)
		win.show()
		app.exec_()


	def tdcExtractSpikes(self, merge_multi_unit = False, multi_unit_clusters = None,
		multi_unut_annot = 'mu', identify_collisions = True, collision_threshold = 10,
		filter_early_late = True):
		'''
		save spike data to file and to object

		collision_threshold : int
			number of samples where two spikes occurring within this many 
			samples or less will be idenified as having collided

			Relevant only if "identify_collisions" is True

		filter_early_late : Boolean
			(default True)
			Filter spikes that are too early or late in the trace to have a full waveform
			This is mainly to avoid errors with further operations with the waveforms.
			Spikes are removed only from the dataframe, and the index is not reset so that
			some sort of record is kept.
		'''

		spks = self.TDCDataIO.get_spikes()

		# Convert to dataframe for convenience, and take only index and cluster label (no jitter)
		self.Spikes = pd.DataFrame(spks[['index', 'cluster_label']])

		# Remove trash spikes (all cluster labels < 0)
		self.Spikes = self.Spikes.query('cluster_label >= 0')

		# Add time column based on index
		sample_rate = self.RawBinChans.__dict__[self.TDC_channel].sampling_rate
		self.Spikes['time'] = self.Spikes['index'] / sample_rate

		# get clus label for each spike
		clus_labels = self.Spikes.cluster_label
		# get annotation for each spike according to clus label
		annotations = self.TDC_clusters[clus_labels]['annotations']

		# add annotations to Spikes DF
		self.Spikes['annotation'] = annotations


		###
		# Merging Multi Units

		if merge_multi_unit:
			if multi_unit_clusters is None:

				self.Spikes.loc[self.Spikes['annotation']==multi_unut_annot, 'cluster_label'] = 111

			else:
				assert isinstance(multi_unit_clusters, list), 'multi unit clusters must be a list'
				multi_unit_spikes = self.Spikes.cluster_label.isin(multi_unit_clusters)
				# Assign cluster 111 to multi unit (using .loc as advised by pandas)
				self.Spikes.loc[multi_unit_spikes, 'cluster_label'] = 111



		clusterLabels = self.Spikes.cluster_label.unique()
		clusterLabels.sort()

		# Filter out multi unit aggregation
		nonMU_CL = clusterLabels[clusterLabels < 111]

		####
		# Identifying Collisions

		if identify_collisions:
			for cl in nonMU_CL: # collisions are expected with multi unit!
				indices = self.Spikes.loc[self.Spikes.cluster_label==cl, 'index']
				diff = np.diff(indices)
				# diff values start on second index, so concat with an initial False and
				# mask of diffs that are below the thresold (ie, are colliding)
				diff_mask = np.r_[False, diff <= collision_threshold]

				# Add collision column and insert diff mask
				self.Spikes.loc[self.Spikes.cluster_label==cl, 'collision'] = diff_mask



		####
		# Spike waveform avgs and std errors etc

		# Cluster averages done after merging if statement so that it flows into 
		# which clusters get averages

		# spans of waveforms for the catalgue used by the peeler
		catN_Left, catN_Right = self.TDCCatConstructor.catalogue['n_left'], self.TDCCatConstructor.catalogue['n_right']	

		# shape -> [n clusters (sorted), peak size (according to catConstructor.catalogue)]
		self.SpikeTemplates = np.squeeze(
			self.TDCCatConstructor.catalogue['centers0'][nonMU_CL, ...]
			)

		self.SpikeAvgs = np.zeros((nonMU_CL.size, catN_Right - catN_Left))
		self.SpikeStd = np.zeros_like(self.SpikeAvgs)
		self.SpikeSem = np.zeros_like(self.SpikeStd)

		# Sometimes a spike peak can be identified earlier than n_left, such that it cannot
		# have a whole waveform. Same with late spikes.
		# Filter out spikes with indexes that are too early or late.

		# presume only one segment, this zero as argument
		trace_length = self.TDCDataIO.datasource.get_segment_shape(0)[0] 

		early_late_idx = self.Spikes[(self.Spikes['index']<np.abs(catN_Left)) | (self.Spikes['index']>(trace_length-catN_Right))]
		self.Spikes.drop(labels=early_late_idx.index, axis=0, inplace=True)





		for i,e in enumerate(nonMU_CL):

			indices = self.Spikes.loc[self.Spikes.cluster_label==e, 'index']

			wfs = np.squeeze(
				self.TDCDataIO.get_some_waveforms(spike_indexes=indices, n_left=catN_Left, 
											n_right=catN_Right)
				)

			mean = np.mean(wfs, axis=0)
			std = np.std(wfs, axis=0)
			stdErr = sem(wfs, axis=0)

			self.SpikeAvgs[i,:] = mean
			self.SpikeStd[i, :] = std
			self.SpikeSem[i, :] = stdErr

		if merge_multi_unit:

			if multi_unit_clusters is None:
				multi_unit_clusters = self.Spikes.loc[self.Spikes.annotation=='mu', 'cluster_label'].unique()

			multi_unit_clusters.sort()

			self.MultiUnitTemps = np.squeeze(
				self.TDCCatConstructor.catalogue['centers0'][multi_unit_clusters, ...]
				)

			indicesMU = self.Spikes.loc[self.Spikes.cluster_label==111, 'index']
			wfsMU = np.squeeze(
				self.TDCDataIO.get_some_waveforms(spike_indexes=indicesMU, n_left=catN_Left,
											n_right=catN_Right)
				)

			self.MultiUnitAvg = np.mean(wfsMU, axis=0)
			self.MultiUnitStd = np.std(wfsMU, axis=0)
			self.MultiUnitSem = sem(wfsMU, axis=0)



		# make self.MultiUnitTemps, Avgs, STd etc
		# temps will be each of the clusters that were merged
		# avgs etc will the same ... single



class rawBinChans:
	'''
	To create objects that store information on channels converted to raw binary files
	'''




	def addChannel(self, chan_name, params):

		self.__dict__.update( { chan_name : chan(params) })


	def updateChannel(self, chan_name, params):

		assert self.listChannels(chan_name=chan_name), f'{chan_name} not available'

		assert isinstance(params, dict), 'params provided are not dictionary'

		getattr(self, chan_name).__dict__.update(params)


	def listChannels(self, chan_name=None):
		'''
		If chan_name == None (default)
			returns list of existing channels

		if chan_name == str (channel name)
			return {processed: True/False, written: True/False}

		'''

		# Produce list of attributes of object that are of type chan (ie, are channels)
		existing_chans = [att for att in dir(self) if isinstance(getattr(self, att), chan)]


		if chan_name is None:
			return existing_chans

		
		else:
			assert isinstance(chan_name, str), 'chan_name must be a string'

			chan_processed = chan_name in existing_chans
			chan_written = self.__dict__[chan_name].written if chan_processed else False

			return dict(processed=chan_processed, written=chan_written)


class chan:
	def __init__(self, params):
		self.written = False
		self.params = params

		self.__dict__.update(params)					

	def pparams(self):
		pprint(self.params)


