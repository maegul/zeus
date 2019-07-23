'''
Library for generating spike times from raw data.

Intention is for spike times to be analysed by Zeus

Hephaistos -> god of craft and metalwork ... this is the hard work!
'''

import pathlib as pthl
from pprint import pprint
import pickle

from . import spk2_mat_read as spkmr

from .utils import funcTracker, getMethods

import neo
import quantities as qnt

import numpy as np
import pandas as pd
from scipy.io.matlab import savemat
from scipy.stats import sem
from scipy import signal as sp_signal

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




def load(filename, tdc_refresh = False, print_state = True):
	with open(filename, 'rb') as f:
		loadedUnit = pickle.load(f)

	if tdc_refresh:
		loadedUnit._tdc_refresh()

	if print_state:
		print('\n*****\n')
		# for k,v in loadedUnit.TDCDataIO.info['HephTDCState'].items():
		# 	print(f'{k} ... {v}')
		print(loadedUnit)

		print('\n*****\n')

	# pprint(loadedUnit.TDCDataIO.info)

	return loadedUnit




class Hephaistos:

	'''
	'''

	def __repr__(self):

		methods = getMethods(Hephaistos)

		# List of methods to report, in order of preferred/necessary execution
		method_report = [
			'readSpkMat',
			'readBasicSpkMat',
			'mkRawBin',
			'tdcInit',
			'tdcPreProcess',
			'tdcSetUpCatalogue',
			'tdcCatPCAClust',
			'tdcPeel',
			'tdcExtractSpikes'			

		]

		conditional_methods = {
			'readBasicSpkMat' : ['writeWaveClusMat', 'quickView'],
			'tdcCatPCAClust' : [
				'tdcFindBetterBounds',
				'tdcViewPCA',
				'tdcOpenCatalogue'
				],
			'tdcPeel' : 'tdcOpenPeeler'
		}

		rep = ''

		for m in method_report:

			if m in methods:

				rep += '{0}\t{1}()\n'.format(
					u"\u2705" if m in self._function_tracking else u"\u274c",
					m
					)

			if m in conditional_methods:
				condMethods = conditional_methods[m]
				if isinstance(condMethods, list):
					for cm in condMethods:
						rep += '\t\t{0}()\n'.format(cm)
				else:
					rep += '\t\t{0}()\n'.format(condMethods)


		return rep



	def __init__(self, data_path=''):

		'''
			Class for managing raw data -> spike data for zeus

			Aim is to store key information (origin of raw data and actual spike data)

			NOTE - design of spk2_mat_read is to only call large data, not store internally

			NOTE - paths are used as relative paths, but also stored as absolute paths
			should it be useful later

			AS A RESULT - YOU MUST USE HEPHAISTOS IN THE DATA_PATH DIRECTORY


			Parameters
			__________

			data_path : str
				path to data file
		'''

		# tracking function executions
		# Useful for repr and user guidance
		self._function_tracking = dict()


		self.Data_path = pthl.Path(data_path).absolute()

		assert self.Data_path.is_file(), 'data path not a file, please provide file with data'


		# Add folder for all tdc data of the same name as the data file
		self._processing_data_path = self.Data_path.parent / self.Data_path.stem # stem removes suffix of file, useful to create new folder 
		if not self._processing_data_path.is_dir():
			self._processing_data_path.mkdir()


		# Make all paths relative

		paths = ('Data_path', '_processing_data_path')

		# Root for data storage is Data_path
		storage_root = self.Data_path.parents[0]

		# Check that currently in storage root (else relative paths won't work)
		assert storage_root == pthl.Path.cwd(), f'You are not in data path directory, relative paths cannot work.\n Move to {storage_root}'

		self._absolute_paths = {
				k: self.__dict__.get(k, None)
				for k in 
				paths
			}

		for p in paths:
			currentPath = self.__dict__[p]
			newPath = currentPath.relative_to(storage_root)
			assert newPath.exists(), f'New relative path ({newPath}) does not exist. \nGenerating relative path failed somehow?'

			self.__dict__[p] = newPath







		self.RawBinChans = rawBinChans()

		# How manage paths?
			# Want to be able to easily pull out containing folders

	@funcTracker
	@doc_string_add(spkmr.spk2MatRead.__init__)
	def readSpkMat(self, **kwargs):

		'''

		Spike 2 Mat Read function args
		_______

		{0}

		'''

		self.Dat = spkmr.spk2MatRead(self.Data_path, **kwargs )

		# Add marker data directly to self so that no need to reload matlab data
		# as self.Dat attribute when loading from pickled save file (as hdf5 not pickle-able)
		self.MarkTimes = self.Dat.MarkTimes
		self.MarkCodes = self.Dat.MarkCodes

	@funcTracker
	@doc_string_add(spkmr.basicMatRead.__init__)
	def readBasicSpkMat(self, **kwargs):

		'''

		Basic Spike 2 Mat read function args
		_______

		{0}
		'''

		self.Dat = spkmr.basicMatRead(self.Data_path, **kwargs)

		# Add marker data directly to self so that no need to reload matlab data
		# as self.Dat attribute when loading from pickled save file (as hdf5 not pickle-able)
		self.MarkTimes = self.Dat.MarkTimes
		self.MarkCodes = self.Dat.MarkCodes


	def writeWaveClusMat(self, channel = 'Ch1'):

		mat_file_path = self._processing_data_path / ('wc_' + self.Data_path.stem)

		savemat(mat_file_path, dict(data = self.Dat.Ch1.GetTraceData(),
									sr = self.Dat.Ch1.SampleFreq),
				appendmat = True
				)



	@funcTracker
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


	def quickView(self, channel = 'Ch1', limit = 10000, viewLimit = 2,
		filter=False, low_cut=300, high_cut=3000, order=5):

		'''
		viewlimit is how many times greater the samples in full view are than limit

		limit is how many samples for downsampled view
		'''


		assert channel in self.Dat._traceChannelNames, 'channel name not in data, check available channels with self.Dat._traceChannelNames'

		# Sample rate
		sr = self.Dat.__getattribute__(channel).SampleFreq
		# Bad to simply store the signal in memory ... could be replaced by an HDF5 read function

		if filter:

			def butter_bandpass(lowcut, highcut, fs, order=5):
				nyq = 0.5 * fs
				low = lowcut / nyq
				high = highcut / nyq
				b, a = sp_signal.butter(order, [low, high], btype='band')
				return b, a

			b,a = butter_bandpass(low_cut, high_cut, sr, order)

			signal = sp_signal.filtfilt(
				b, a,
				self.Dat.__getattribute__(channel).GetTraceData()
				)

		else:
			signal = self.Dat.__getattribute__(channel).GetTraceData()

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
		


	@funcTracker
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

			# Custom state recording
			self.TDCDataIO.info['HephTDCState'] = {}
			self.TDCDataIO.info['HephTDCState'].update(
					dict(DataIO = True)
				)
			self.TDCDataIO.flush_info()

		else:
			previous_channel = self.TDCDataIO.info['spkmr_channel']

			# assert previous_channel == channel, f'Channel arg ({channel}) is not same as previously encoded TDC channel ({previous_channel})'

			# recreating variable used above for future reference below
			self.TDC_channel = previous_channel

			print(f'DataSource already set. {previous_channel} used. \nProcess data again, or, Alter path to start again.\n')
			print(self.TDCDataIO)


	@funcTracker
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

		# for quick reference to measure of noise
		self.TDC_signals_mads, self.TDC_signals_medians = self.TDCCatConstructor.signals_mads, self.TDCCatConstructor.signals_medians

		# Custom state recording
		self.TDCDataIO.info['HephTDCState'].update(
				dict(PreProcess = True)
			)
		self.TDCDataIO.flush_info()


		pprint(self.TDCCatConstructor)
		pprint(self.TDCCatConstructor.info)


	@funcTracker
	def tdcSetUpCatalogue(self, n_left=-30, n_right=30, align_waveform=True, 
		alien_value_threshold = 100, catalogue_mode='rand'):
		'''
		Setup the catalogue and run both PCA and clustering

		n_left and n_right expressed in milliseconds, converted to samples in code

		alien_value_threshold is expressed MAD (I think)


		catalogue_mode : str
			mode for extract_some_waveforms
			'rand' - random
			'all' - all peaks
		'''

		assert n_left < 0 and n_right > 0, 'n_left must be negative, n_right positive'

		self.TDC_nLeft = n_left
		self.TDC_nRight = n_right


		self.TDCCatConstructor.extract_some_waveforms(n_left = n_left, n_right = n_right,
			align_waveform=align_waveform,
			mode = catalogue_mode)

		self.TDCCatConstructor.extract_some_noise()

		self.TDCCatConstructor.clean_waveforms(alien_value_threshold=alien_value_threshold)

		# Custom state recording
		self.TDCDataIO.info['HephTDCState'].update(
				dict(SetUpCatalogue = True)
			)
		self.TDCDataIO.flush_info()



	@funcTracker
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

		# Custom state recording
		self.TDCDataIO.info['HephTDCState'].update(
				dict(PCAClust = True)
			)
		self.TDCDataIO.flush_info()


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


	@funcTracker
	def tdcPeel(self, external=None, insert_actual_MAD_values = True):
		'''
		Set up peeler, with catalogue, and run

		external : str | Hephaistos object
			Pass an external unit, as a path to be loaded or directly
			The catalogue of this unit will be used for the purposes of peeling

		insert_actual_MAD_values : boolean
			Requires self.tdcPreProcess() to have been run (to estimate noise)
			Use noise estimages (medians and MAD values) of the CURRENT data in 
			the catalogue of the external data.
			This should adjust for any changes in amplitude that have occurred
			over time.
			The values are dynamically inserted into the self.TDCCatalogue 
			object, which is derived from the external unit.

		'''

		if not external:
			self.TDCCatConstructor.make_catalogue_for_peeler()

			# Useful to store here as I don't think tdc stores it persistently with the peeler
			self.TDCCatalogue = self.TDCDataIO.load_catalogue()
			self.TDC_clusters = self.TDCCatalogue['clusters']

		else:
			if isinstance(external, str):

				print(f'loading unit from path: {external}\n')
				external_unit = load(external, tdc_refresh=True)

			elif external.__class__.__name__ == 'Hephaistos':
				external_unit = external

			self.TDCCatalogue = external_unit.TDCDataIO.load_catalogue()
			self.TDC_clusters = self.TDCCatalogue['clusters']

			self.TDC_external_catalogue = external_unit.Data_path	

			# Insert noise measures of actual data into imported/external template
			if insert_actual_MAD_values:
				assert hasattr(self, 'TDC_signals_mads') and hasattr(self, 'TDC_signals_medians'), 'Need to have run self.tdcPreProcess to get estimate of noise and MADS'

				self.TDCCatalogue['signals_mads'] = self.TDC_signals_mads
				self.TDCCatalogue['signals_medians'] = self.TDC_signals_medians

		self.TDCPeeler = Peeler(self.TDCDataIO)
		self.TDCPeeler.change_params(catalogue=self.TDCCatalogue)

		self.TDCPeeler.run()

		# Custom state recording
		self.TDCDataIO.info['HephTDCState'].update(
				dict(Peel = True)
			)
		self.TDCDataIO.flush_info()
	



	def tdcOpenPeeler(self):
		'''
		Open peeler window
		'''

		app = pg.mkQApp()
		win = tdc.PeelerWindow(dataio=self.TDCDataIO, catalogue=self.TDCCatalogue)
		win.show()
		app.exec_()


	@funcTracker
	def tdcExtractSpikes(self, merge_multi_unit = False, multi_unit_clusters = None,
		multi_unut_annot = 'mu', identify_collisions = True, collision_threshold = 5,
		filter_early_late = True, 
		drop_duplicate_spikes=True, drop_collisions=True):
		'''
		save spike data to file and to object

		collision_threshold : int
			number of samples where two spikes occurring within this many 
			samples or less will be idenified as having collided

			Relevant only if "identify_collisions" is True

		drop_collisions : Boolean
			If collisions identified, remove them from self.Spikes

		drop_duplicate_spikes : Boolean
			If two spikes are of the same cluster

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
		# Get indices corresponding to each clus label
		clus_label_idxs = clus_labels.replace(self.TDCCatalogue['label_to_index'])

		# get annotation for each spike according to clus label 
		annotations = self.TDC_clusters[clus_label_idxs]['annotations']

		# add annotations to Spikes DF
		self.Spikes['annotation'] = annotations


		clusterLabels = self.Spikes.cluster_label.unique()
		clusterLabels.sort()


		###
		# Merging Multi Units

		# Generate list of multi unit clusters from annotations
		if multi_unit_clusters is None:
			multi_unit_index = self.Spikes['annotation']==multi_unut_annot
			multi_unit_clusters = list(
				self.Spikes.loc[multi_unit_index, 'cluster_label'].unique()
				)

		# Checking that all cluster labels provided are actually in spikes cluster_labels
		else:
			assert isinstance(multi_unit_clusters, list), 'multi_unit_clusters must be a list'

			# clusters = self.Spikes.cluster_label.unique()

			for muc in multi_unit_clusters:
				assert muc in clusterLabels, f'"{muc}" not in cluster labels ({clusterLabels})'



		if merge_multi_unit:

			assert isinstance(multi_unit_clusters, list), 'multi unit clusters must be a list'

			multi_unit_spikes = self.Spikes.cluster_label.isin(multi_unit_clusters)
			# Assign cluster 111 to multi unit (using .loc as advised by pandas)
			self.Spikes.loc[multi_unit_spikes, 'cluster_label'] = 111




		# Filter out multi unit aggregation
		nonMU_CL = [
			cl 
			for cl in clusterLabels
			if cl not in multi_unit_clusters+[111] # adding 111 if any merging has occurred
		]


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


		# Drop spikes where they are the same cluster at exactly the same time
		# Such duplications are likely the result of a cluster with variable amplitude
		# Such clusters could be a problem, but are possible
		if drop_duplicate_spikes:
			# Looking for duplication on time index and cluster label
			# Then, to allow duplication for multi unit, filtering out either multi unit cluster labels or the merged 111 label
			duplicates_boolean = ( 
				(self.Spikes.duplicated(subset=['index', 'cluster_label'])) 
				& (~self.Spikes.cluster_label.isin( multi_unit_clusters + [111] ))
				)

			duplicate_idxs = duplicates_boolean[duplicates_boolean].index #odd syntax, but works cuz series is boolean
			self.Spikes.drop(labels = duplicate_idxs, axis=0, inplace=True)

		if drop_collisions:
			assert identify_collisions, 'You need to identify collisions to drop them'

			collision_idxs = self.Spikes.query('collision == True').index
			self.Spikes.drop(labels = collision_idxs, axis=0, inplace=True)



		####
		# Spike waveform avgs and std errors etc

		# Cluster averages done after merging if statement so that it flows into 
		# which clusters get averages

		# spans of waveforms for the catalgue used by the peeler
		# catN_Left, catN_Right = self.TDCCatConstructor.catalogue['n_left'], self.TDCCatConstructor.catalogue['n_right']	
		catN_Left, catN_Right = self.TDCCatalogue['n_left'], self.TDCCatalogue['n_right']	

		df_col_sample_idx = np.arange(catN_Right - catN_Left)

		# Instantiate empty dataframes, with appropriate amount of columns (one for each sample)
		self.SpikeTemplates, self.SpikeAvgs, self.SpikeStd, self.SpikeSem = (
			pd.DataFrame(columns = df_col_sample_idx).rename_axis('cluster_label')
			for i in range(4)
			)
		# shape -> [n clusters (sorted), peak size (according to catConstructor.catalogue)]
		# self.SpikeTemplates = np.squeeze(
		# 	self.TDCCatalogue['centers0'][nonMU_CL, ...]
		# 	)

		# self.SpikeAvgs = np.zeros((nonMU_CL.size, catN_Right - catN_Left))
		# self.SpikeStd = np.zeros_like(self.SpikeAvgs)
		# self.SpikeSem = np.zeros_like(self.SpikeStd)



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

			# order of labels in tdc may be not in order of number
			# label_to_index records (I think) what the correspondence is, as key-value dict pairs
			label_to_idx = self.TDCCatalogue['label_to_index'][e]

			template = np.squeeze(
				self.TDCCatalogue['centers0'][label_to_idx,...]
				)
			# Add rows with .loc[index/cluster_label]
			self.SpikeTemplates.loc[e] = template
			self.SpikeAvgs.loc[e] = mean
			self.SpikeStd.loc[e] = std
			self.SpikeSem.loc[e] = stdErr


		# Recording templates and averages for the multi units

		if merge_multi_unit or (len(multi_unit_clusters) > 0):
			# multi_unit_clusters should be filled, if no arg passed, to be those marked with 'mu' in annotations


			multi_unit_clusters.sort()

			self.MultiUnitTemps = np.squeeze(
				self.TDCCatalogue['centers0'][multi_unit_clusters, ...]
				)

			# Allows for not merging multi units in Spikes table, but having a merge for mean spike shape
			# Perhaps in future, I want average for each mu cluster, but I suppose the template can fulfill that role
			mu_filt = (self.Spikes.annotation == 'mu') | (self.Spikes.cluster_label == 111)

			indicesMU = self.Spikes.loc[mu_filt, 'index']
			wfsMU = np.squeeze(
				self.TDCDataIO.get_some_waveforms(spike_indexes=indicesMU, n_left=catN_Left,
											n_right=catN_Right)
				)

			self.MultiUnitAvg = np.mean(wfsMU, axis=0)
			self.MultiUnitStd = np.std(wfsMU, axis=0)
			self.MultiUnitSem = sem(wfsMU, axis=0)


		# Custom state recording
		self.TDCDataIO.info['HephTDCState'].update(
				dict(ExtractSpikes = True)
			)
		self.TDCDataIO.flush_info()


		# make self.MultiUnitTemps, Avgs, STd etc
		# temps will be each of the clusters that were merged
		# avgs etc will the same ... single



	def _dropSpikesGenMask(self, cell_base=None, cell_collide = None, threshold = None):
		
		
		mask = (cell_collide > 0)
		
		for v in cell_base:
		
			newMask = (cell_collide > (v + threshold)) | ((v - threshold) > cell_collide )
			mask = mask & newMask
			
			
		return mask


	def tdcDropSpikesCell2Cell(self, cell_base=None, cell_collide = None, threshold = 0.0001):
		'''
		Drops spikes from one cell that collide with spikes from another
		Spikes are dropped from self.Spikes

		Parameters
		_____
		cell_base : int or str
			Cell label with which spikes may collide.  
			Spikes from this cell are not dropped 

		cell_collide : int or str
			Cell label for the spikes that may collide with cell_base spikes and will be dropped.

		threshold : float
			Time window within which spikes will be considered to have collided.
			Time window is made as (spike_time - threshold, spike_time + threshold)
		'''
		
		cell0 = self.Spikes.query('cluster_label == @cell_base').time
		cell0 = cell0.values
		
		cell1 = self.Spikes.query('cluster_label == @cell_collide').time
		cell1_idx = cell1.index
		cell1 = cell1.values
		
		mask = self._dropSpikesGenMask(cell_base=cell0, cell_collide=cell1, threshold=threshold)
		
		index_mask = cell1_idx[~mask]

		self.Spikes.drop(labels = index_mask, axis=0, inplace=True)
		


	def plotSpikeShape(self, use_template=False):
			
		spikeShapes = self.SpikeTemplates if use_template else self.SpikeAvgs
		
		for cell in self.SpikeAvgs.index:
			
			plt.plot(spikeShapes.loc[cell,:].values, label=f'Cell {cell}')
			
		plt.legend()


	def plotSpikeAvgsTemps(self, n_cols = 2, figsize = (12,6)):
		
		spikeShapes = self.SpikeAvgs
		spikeTemps = self.SpikeTemplates
		
		n_cells = self.SpikeAvgs.index.size
		
		n_rows = n_cells // n_cols + n_cells%n_cols
		
		plt.figure(figsize=figsize)
		
		for i, cell in enumerate(self.SpikeAvgs.index):
			
			plt.subplot(n_rows, n_cols, (i + 1))
			
			plt.plot(spikeShapes.loc[cell,:].values, label=f'Avg')
			plt.plot(spikeTemps.loc[cell,:].values, label='template')
			plt.title(f'Cell {cell}')
			
			if i == 0:
				plt.legend()
			
			
		plt.tight_layout()



	def save(self, file_name = None):
		'''
		pickle the unit
		'''

		if file_name is None:
			file_name = self.Data_path.parent / ('Heph_' + self.Data_path.stem + '.pkl')

		elif pthl.Path(file_name).suffix != '.pkl':
			file_name = self.Data_path.parent / pthl.Path(file_name).with_suffix('.pkl')

		self.SavePath = file_name.absolute()

		temp_path = self.SavePath.with_suffix(file_name.suffix + '.tmp')
		

		
		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(self.SavePath)

		print(f'Saved pickle file to {str(file_name.parent)} as {str(file_name.name)}')


	def __getstate__(self):
		state = self.__dict__.copy()

		# HDF5 Dat attribute prevents pickling, and is also redundant 
		# (it's) a read of data files

		try:
			del state['Dat']
		except KeyError:
			pass

		# Get rid of the TDC attributes, as they occupy much redundant space
		# Essentially the whole contents of the processing path is encoded in the
		# pickle (through all the memmap objects probably)
		# They will be re instantiated in setstate
		for tdcAttr in ['TDCDataIO', 'TDCCatConstructor', 'TDCPeeler']:
			try:	
				del state[tdcAttr]
			except KeyError:
				pass

		return state


	def __setstate__(self, state):
		# re-instantiate TDC objects
		# Just do all three, irrespective of whether they were already done
		# Except for dataio, only if tdcInit has been run (maybe just run tdcInit, then make cat and peeler, as they rely only on dataio object)

		self.__dict__.update(state)



	def _tdc_refresh(self):

		self.tdcInit()

		self.TDCCatConstructor = CatalogueConstructor(dataio = self.TDCDataIO)
		self.TDCPeeler = Peeler(self.TDCDataIO)




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



# 
# Utility functions
# 

def mkCorrelogram(spkTimesA, spkTimesB=None):

	a = spkTimesA
	b = spkTimesB if spkTimesB is not None else spkTimesA

	xs, ys = np.meshgrid(a,b)

	# subtract a times from b times to get how long after a spikes b spikes occur (7 - 3 = 4)
	diffs = ys - xs

	# Ie, this is an autocorrelation
	if spkTimesB is None:
		mask = np.ones_like(xs, dtype='bool' )
		np.fill_diagonal(mask, False)

		diffs = diffs[mask]

	return diffs.flatten()


def plotCorrelogram(unit, cells = [], histtype='stepfilled', width=15, bin_width=0.1, figsize=None):
	
	assert len(cells) <3, 'cannot yet do more than 2 cells for a grid plot'
	
	spikeTimes = [
		unit.Spikes.query('cluster_label==@c').time.values
		for c in cells
	]
	
	diffs = mkCorrelogram(*spikeTimes)
	# Work in milliseconds
	diffs *= 1000

	strt = width/2
	end = width/2

	bins = np.arange(-(strt+(bin_width/2)), end+bin_width+(bin_width/2), bin_width)
	
	if figsize is not None:
		plt.figure(figsize=figsize)
	plt.hist(diffs, bins=bins, histtype=histtype)



def plotCompareSpikesRuns(unit_a, unit_b, cells=[], run_labels=[], figsize=[12,6]):
	
	if isinstance(unit_a, str):
		unit_a = load(unit_a, print_state=False)
		
	if isinstance(unit_b, str):
		unit_b = load(unit_b, print_state=False)
		
	n_cells = len(cells)
	
	
	
	if len(run_labels) == 0:
		run_labels = ['Run a', 'Run b']
	
	
	plt.figure(figsize=figsize)
	
	for i,c in enumerate(cells):
		
		plt.subplot(n_cells, 3, (i*3+1))
		
		plt.plot(unit_a.SpikeTemplates.loc[c,:].values, label=run_labels[0], color='C0')
		plt.plot(unit_b.SpikeTemplates.loc[c,:].values, label=run_labels[1], color='C2')
		
		plt.title(f'Cell {c} - Temp')
		
		plt.legend()
		
		
		plt.subplot(n_cells, 3, (i*3+2))
		
		plt.plot(unit_a.SpikeAvgs.loc[c,:].values, label=run_labels[0], color='C1')
		plt.plot(unit_b.SpikeAvgs.loc[c,:].values, label=run_labels[1], color='C3')
		plt.title('Avg')
		
		plt.legend()
		
		
		plt.subplot(n_cells, 3, (i*3+3))

		plt.plot( unit_b.SpikeTemplates.loc[c,:].values - unit_a.SpikeTemplates.loc[c,:].values, 
				 label=f'Templ', color='C2')
		
		plt.plot( unit_b.SpikeAvgs.loc[c,:].values - unit_a.SpikeAvgs.loc[c,:].values, 
				 label=f'Avg', color='C3')
		
		plt.axhline(y=0, linestyle=':', color='0.65')
		
		plt.title(f'Difference ({run_labels[1]} - {run_labels[0]})')
		plt.legend()
		
	plt.tight_layout()
	