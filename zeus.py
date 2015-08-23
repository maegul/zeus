
import math
import numpy as np
import scipy as sp

import scipy.stats as stats
from scipy.ndimage import gaussian_filter1d
import numpy.fft as fft

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.style.use('ggplot')

import pandas as pd

import os
import glob

import cPickle

def load_raw_data(dir):
    
    """
    Loads files from the directory provided.
    
    Loads only txt files.
    
    Sorts and assigns expecting a number of either markerss, spikess times and spikess shape.
    
    Returns a dictionary of the imported data types
    """

    current_dir = os.getcwd()    
    
    os.chdir(dir)
    
    files = glob.glob('*.txt')
    file_names = []
    data = {}
    
    for f in files:
        f_name = f.lower()
        
        if f_name.find('mark') > -1:
            data['markers'] = np.loadtxt(f_name, skiprows=1)
            file_names.append(f)
            
        elif f_name.find('spike') > -1:
            data['spikes'] = np.loadtxt(f_name, skiprows=1)
            file_names.append(f)
            
        elif f_name.find('shape') > -1:
            data['shape'] = np.loadtxt(f_name, skiprows=1)
            file_names.append(f)
            
            
    os.chdir(current_dir)

            
    if len(data.keys()) != len(files):
        mesg = 'Not all of your file names are recognised; they may not have been imported appropriately'
        mesg2 = 'File names must contain the key words "mark", "spike" and/or "shape"'
        print mesg
        print mesg2
        print '\nFollowing files loaded successfully:\n'
        for i in file_names: print(i)
        return data

    
    elif len(data.keys()) == len(files):
        print 'All files imported and assigned'
        print '\nFollowing files loaded successfully:\n'
        for i in file_names: print(i)
        return data
        
def save(dataset):
    """
    Pickles a class instantiation of `Zeus`
    """
    
    directory = dataset.parameters['directory']
    data_label = dataset.parameters['data_label']
    
    
    with open('%s%s_object.pkl'%(directory, data_label), 'wb') as f:
        cPickle.dump(dataset, f)
        

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
        return cPickle.load(f)
        
        


#def plotform(ax, tickint=False):
#    
#
#
#    # remove top and right axes
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('left')
#    
#    # make ticks outward
#    ax.tick_params('both', direction='out')
#    
#    
##    fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
##    'weight' : 'normal', 'size' : 12}
#    
##    if tickint:
##        p.set_xticklabels(p.get_xticks().astype('int'), fontProperties)
##        p.set_yticklabels(p.get_yticks().astype('int'), fontProperties)
##    else:
##        p.set_xticklabels(p.get_xticks(), fontProperties)
##        p.set_yticklabels(p.get_yticks(), fontProperties)
#        
#        
#    #set the style of the major and minor grid lines, filled blocks
#    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
#    ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
#    ax.patch.set_facecolor('0.85')
#    ax.set_axisbelow(True)
#    
#    #set minor tick spacing to 1/2 of the major ticks
#    ax.xaxis.set_minor_locator(MultipleLocator( (plt.xticks()[0][1]-plt.xticks()[0][0]) / 2.0 ))
#    ax.yaxis.set_minor_locator(MultipleLocator( (plt.yticks()[0][1]-plt.yticks()[0][0]) / 2.0 ))
#    
#    #remove axis border
##    for child in ax.get_children():
##        if isinstance(child, matplotlib.spines.Spine):
##            child.set_alpha(0)
#       
#    #restyle the tick lines
#    for line in ax.get_xticklines() + ax.get_yticklines():
#        line.set_markersize(5)
#        line.set_color("gray")
#        line.set_markeredgewidth(1.4)
#    
#    #remove the minor tick lines    
#    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
#        line.set_markersize(0)
#
#    
#    
#    if ax.legend_ <> None:
#        lg = ax.legend_
#        lg.get_frame().set_linewidth(0)
#        lg.get_frame().set_alpha(0.5)
#        
#    plt.tight_layout()


def steppify(arr, axis='x'):
    """
    expands datasets to create de facto step (step-post) plots
    """
    
    if axis == 'x':
        newarr = np.r_[arr[0], np.dstack((arr[1:], arr[1:])).flatten()]
    
    elif axis == 'y':
        newarr = np.r_[np.dstack((arr[:-1], arr[:-1])).flatten(), arr[-1]]
    
    else:
        print 'your axes in steppify are improperly identified'

    return newarr


#==============================================================================
# Generalised Analysis Functions
#==============================================================================

def fourier(data, temp_freq, axis, output = 'amplitude'):
    """
    Generalised fourier analysis for drifting grating stimuli.
    
    
    Parameters
    __________
    
    axis : int
        Axis along which to perform fourier analysis
    
    """
        
    
    # take largest possible multiple of F1 from PSTH.
    # Generate freq and fft
    # generate amplitude
    # return amplitude, F0, F1 and F2 values
    
    

def bootstrap(data, alpha=0.05, n_bootstrap = 2000, func=None, **func_args):
    
    """Generalised bootstrapping function for data structured as con x trial x bins(time).
    
    Parameters
    __________
    data : array (3-D)
        Data to be bootstrapped.  Must be three dimensional with 
        trials along the second dimension (conditions on the first, time on the third)
    
    alpha : float
        Determines the confidence interval, where the probability of a type I error
        is `alpha`.  confidence intervals will be returned at the (100*alpha/2) and
        (100*(1-alpha/2)) percentiles.  Default is 0.05 (ie, 5%).
        
    n_bootstrap : int
        Number of iterations for the bootstrap distribution.  Default is 2000.
        Advisable that between 1000 and 10,000.
        
    func : function
        A function to operate on the bins or time dimension of the data before
        confidence intervals are calculated (philosophy being that the most
        non-normal data ought to be bootstrapped the most).
        
        Must accept 3-D data set, operate on the third dimension (axis=2),
        and return a 3-D data set.  If a pre-defined function is sought, wrapping
        in a custom function to meet these needs is advised.
        
    
    Returns
    _______
    out : tuple of arrays
        Tuple of 2 arrays, first the positive confidence interval, second, negative
        condidence interval, both at percentiles corresponding to a `alpha` chance
        of a type I error.  Each array is 2-D, with conditions on first axis, and bins
        (or time) on the second axis.
    """
    
    assert data.ndim == 3, 'Data is not 3-dimensional.  Function only works for 3-D data.'    
    
    # Trials form the second dimension
    n_trials = data.shape[1]
    
    # generate randomised bootstrap resamples as random indices
    bootstrap_index = np.random.randint(0, n_trials, 
                                        (n_trials, n_bootstrap) )
    
    # For each bin in the histogram, randomly samples from the results
    # of each trial and repeats, effectively, n_bootstrap times                
    trials_bootstrap = data[:, bootstrap_index, :]
    
    # dimension one is the trials, zero is the conditions; this averaging 
    # goes across the trials creating a PSTH for each condition, and,
    # importantly, for each bootstrap resample
    avg_bootstrap = trials_bootstrap.mean(axis=1)
    
    if func:
        avg_bootstrap = func(avg_bootstrap, **func_args)
        
    # find percentile values for each bin along the bootstrap resamples,
    # which are on axis 1                                              
    CI_pos = np.percentile(avg_bootstrap, 100*(1 - (alpha/2.)), 
                                axis=1)
    CI_neg = np.percentile(avg_bootstrap, 100*(alpha/2.), 
                                axis=1)


    return CI_pos, CI_neg
    


## Class Definitions

class Zeus:
    
    
    
    def __init__(self, data_label, dir):
        
        """
        Instatiates the output of `zeus.load` as self.data.
        
        `zeus.load` argument, directory, is passed from class instantiation.
        
        For markers, spikes and shape, an attribute is directly .
        
        self.parameters instantiated, and directory argument appended.
        
        """
        
        # Load files and assign to attributes
        self.data = load_raw_data(dir)
        
        try:
            self.markers = self.data['markers']
        except:
            pass
        
        try:
            self.spikes = self.data['spikes']
        except:
            pass

        try:
            self.shape = self.data['shape']
        except:
            pass


        curr_dir = os.getcwd()
        self.parameters = {}
        
        self.parameters['directory'] = dir
        self.parameters['ops_directory'] = curr_dir
        
        self.parameters['data_label'] = str(data_label)
        
        if (self.parameters['directory'][:3] == '../') or (self.parameters['directory'][:2] == './'):
            print '\n\nFull directories should be employed for more accurate data records \n\n'
        
        print ('\nRun the following methods:\n self.sort()\n self.conditions()\n'
                ' self._analyse() \n\n'
                'Saving and plotting functions may then be used.')
        
        
    def _sort(self, conditions = 9, trials = 10, stim_len = None,
              bin_width=0.02, auto_spont=False, sigma=3, alpha=0.05, 
              n_bootstrap=2000):
        
        """Sorts spikes into conditions and trials on the basis of the markers
        loaded by `zeus.__init__` and the expected amount of conditions and
        trials defined by the passed arguments.  Then, generates histograms, 
        or PSTHS, for each of the conditions, averaged across the trials.
        
        Parameters
        __________
        conditions : int
            Amount of stimulus conditions in the experiment.
            
        trials : int
            Amount of trials of each condition in the experiment.
            
        stim_len : float, optional
            Default is `None`.
            Defines the length of the stimulus for a single trial.
            If not defined, it is calculated as the maximum difference in time
            between markers.
            
        bin_width : float, default = 0.02
            Determines the width of the bins in histograms or PSTHs
            
        auto_spont : Boolean
            Not yet complete.  Determines whether to calculate a spontaneous rate automatically.
            Currently, only sure method for this is to count spikes before the first stimulus.
            More methods are likelty to be used later.
            
        sigma : float
            Defines the standard deviation, in bins, of the gaussian kernel applied
            to time-series data to smoothen the data.  To obtain the standard 
            deviation in units of absolute time -> `self.bin_width * sigma`.
            
        alpha : float
            To be passed to global method `bootstrap()`
            Determines the confidence interval, where the probability of a type I error
            is `alpha`.  confidence intervals will be returned at the (100*alpha/2) and
            (100*(1-alpha/2)) percentiles.  Default is 0.05 (ie, 5%).
        
        n_bootstrap : int
            To be passed to global method `bootstrap()`
            Number of iterations for the bootstrap distribution.  Default is 2000.
            Advisable that between 1000 and 10,000.

        
        Returns
        _______
        
        self.conditions_trials : dict of lists
            Three dimensional array containing the REAL spike times of 
            each trial of each condition.
            
        self.conditions_trials_zeroed : dict of lists
            Three dimensional array containing the spike times, relative to the
            commencement of the relevant trial, of each trial of each 
            condition.
            
        self.parameters : dict, pre-instantiated
            'conditions' (Number of conditions)
            'trials' (Number)

        self.bin_width : float
            Width of the bins used in generating the histograms.  Defines in seconds.
            
        self.bins : array (1-D)
            The bind boundaries used in generating the histograms
            
        self.conditions_trials_hist : array (3-D)
            PSTH/histograms for each trial for each condition.  
            Provides spike counts for each trial
            First dimension (axis=0) - conditions
            Second dimension (axis=1) - trials
            Third dimension (axis=2) - bins (or time)
            
        self.conditions_hist_mean : array (2-D)
            As `self.conditions_trials_hist` but averaged across trials, 
            provides the mean spike count.
            Second dimension (axis=1) - bins (or time),
            
        self.conditions_hist_stderr : array (2-D)
            As `self.conditions_hist_mean` but with standard error of the mean rather than the mean
            itself provided.
            
        self.conditions_trials_sdf : array (3-D)
            As `self.conditions_trials_hist` but with gaussian smoothing applied
            to each trial. 
            This gaussian is 5-sigmas wide on each side (10 total) and has 
            `sigma` standard deviation. It is the same gaussian applied for 
            the spike density function.
            
        self.spike_density_function : array (2-D)
            As `self.conitions_hist_mean` but with gaussian smoothing applied to
            the mean values.
            
        self.CI_pos : array (2-D)
            As `self.conditions_hist_stderr` but with bootstrap derived confidence
            intervals, for the positive side.  
            Accompanies `self.spike_density_function`, as the intervals are derived
            from a population sampled from smoothened values.  
            Generated by global method `bootstrap()`.

        """
        # Loading parameters to parameters dictionary
        self.parameters['conditions'] = conditions
        self.parameters['trials'] = trials

        self.parameters['sdf_alpha'] = alpha
        self.parameters['sdf_number_bootstrap'] = n_bootstrap     
        self.parameters['sdf_sigma'] = sigma       

        self.parameters['bin_width'] = bin_width
        
#==============================================================================
# Post-marker buffer & Spontaneous calculation needed.
#==============================================================================
        
        assert self.markers.size == (conditions * trials), (
            'The number of recorded markers' 
            'does not match your dataset specifications.'
            '\n#Conditions * #trials != #markers.\n'
            'Markers will probably have to be inserted.  Run self._marker_diag() to diagnose')
        
        
        # determine stimulus length
        if stim_len:
            self.parameters['stimulus_length'] = float(stim_len)
        else:
            marker_diff = np.subtract( self.markers[1:], self.markers[:-1])
            stim_len = np.max( marker_diff)
            
            self.parameters['stimulus_length'] = stim_len
            
            print('\nstimulus length has been calculated from the set of markers; ' 
                  'taking the maximum distance between all adjacent pairs of markers')
            
        
        # Initialising 
        self.conditions_trials = {}
        self.conditions_trials_zeroed = {}
        
        ## Sorting into trials and conditions 
        # one condition at a time
        for c in range(conditions):
            
            ct = [] # temp for condition_trials
            ctz = [] # temp for condition_trials_zeroed
            
            # one trial at a time
            for t in self.markers[c::conditions]: # Goes thorugh every trial of each condition
                
                # extract spike times for each trial
                trial = np.extract( (self.spikes >= t) & (self.spikes < (t+stim_len)), 
                                   self.spikes)
                
                # produces list of all spike times for the condition, each trial listed
                # consecutively
                ct.append(trial)
                # zero to trial onset
                ctz.append( (trial - t))
                
            # list of spike times listed in dictionary and keyed by condition number    
            self.conditions_trials[c] = ct
            self.conditions_trials_zeroed[c] = ctz
            
#==============================================================================
# Generating PSTHs - Previously in analyse()
#==============================================================================

        
        self.bin_width = bin_width
        
        
        # bin boundaries, in time.  Fundamental to basic PSTH analysis and so made an attribute
        self.bins = np.arange(0, stim_len, bin_width)
        
        n_con = conditions
        n_trial = trials
        
        # Initialising.  Will contain, for each condition, a PSTH for each trial
        self.conditions_trials_hist = np.zeros((n_con, n_trial, self.bins.size - 1))
        
        # One condition, and trial, at a time
        for c in range(n_con):
            
            trials_hist = np.zeros((n_trial, self.bins.size - 1))
            
            for t in range(n_trial):
                
                trials_hist[t] = np.histogram(self.conditions_trials_zeroed[c][t], self.bins)[0]
                
            # PSTH for each trial added for the condition    
            self.conditions_trials_hist[c] = trials_hist
            
        # Average and Standard Error, for each bin, calculated across the trials (second axis)    
        self.conditions_hist_mean = np.mean(self.conditions_trials_hist, axis=1)
        self.conditions_hist_stderr = stats.sem(self.conditions_trials_hist, axis= 1, ddof=0)


#==============================================================================
# Create more robust spontaneous function
#==============================================================================
        # Calculate spontaneous rate from unstimulated portion of the record
        if auto_spont:
            index = self.spikes < self.markers[0]
            self.spont_rate = self.spikes[index].size / self.markers[0]



#==============================================================================
# Spike Density Function generation - previously _sdf()            
#==============================================================================

        
        def gauss_smooth(data, sigma):
            """
            Wrapper for guassian smooting, to be passed to bootstrap function
            """
                        
            
            # make the kernel 5 sigmas wide in each direction
            kernel = stats.norm.pdf(np.arange(-5*sigma, (5*sigma)+1), scale=sigma)
            
            return sp.ndimage.convolve1d(data, kernel, axis=2)

        kernel = stats.norm.pdf(np.arange(-5*sigma, (5*sigma)+1), scale=sigma)

                                              
        self.CI_pos, self.CI_neg = bootstrap(self.conditions_trials_hist, alpha=alpha,
                                             func = gauss_smooth, **{'sigma':sigma})
                               
        # Gaussian smoothing for each trial
        self.conditions_trials_sdf = sp.ndimage.convolve1d(self.conditions_trials_hist,
                                                            kernel, axis=2)
        # Gaussian smoothing over condition averages, for a spike density function
        self.spike_dens_func = sp.ndimage.convolve1d(self.conditions_hist_mean, 
                                                     kernel, axis=1)
            
            
    def _marker_diag(self, auto_replace = False):
        """Attempts to pinpoint where markers have been lost or missed in the marker stream.
        
        Parameters
        __________
        auto_replace : Boolean
            Not yet active.  Intention is to allow auto-replacenemtn of missing markers
            
        Returns
        _______
        table : dict
            Dictionary containing table of where markers are missing and key numbers regarding
            how many markers are missing.
        """
        
        num_missing_markers = (self.parameters['conditions']*self.parameters['trials']) - self.markers.size        
        
        if self.parameters['conditions']*self.parameters['trials'] == (self.markers.size):
            print(
                'Stimulus parameters appear to be appropriate.  Running this function is '
                'unlikely to be necessary.  If anomalies are found, timing errors during the '
                'experimen are likely to be the culprit, and perhaps irretrivable.')
        
        output = {}

        marker_diff = np.r_[0.0, ( np.subtract( self.markers[1:], self.markers[:-1]))]
        
        multiple = np.around((marker_diff / float(np.median(marker_diff))) - 1, 1)
    
        anomaly = np.greater(marker_diff, np.ones_like(marker_diff) * 1.2 * np.median(marker_diff))
        
        table = pd.DataFrame({'markers': self.markers, 'multiple': multiple, 
                              'difference': marker_diff, 'anomalous': anomaly})

        
        # Anomalous markers, require insertion
        
        anomalies = table[table['anomalous'] == True]
        output['bad_marks'] = anomalies
        output['bad_marks_index'] = anomalies.index
        output['num_mark_found'] = np.around(anomalies['multiple'].sum(), 1)
        output['num_missing_mark'] = num_missing_markers
        output['num_missing_mark_at_beg/end'] = max(output['num_missing_mark']-output['num_mark_found'], 0)
        output['maximum_marer_diff'] = marker_diff.max()
        output['median_marker_diff'] = np.median(marker_diff)
            
#==============================================================================
# Add:
#            
#            Auto Replace function.  Create a parameter reflecting, and do not touch the original
#            file.
#            Capacity to tell whether marker missing from front or end of experiment.
#            Use average spike rates to test, or require feedback from teh original data file.        
#==============================================================================
            
        return output
            

            
            
             
             

    def _conditions(self, beg=-90, intvl=20, con_type='orientation', stim='bar', 
                    biphasic=True, unit='deg', con_list=[], temp_freq = 2):
                        
        """Generates condition and stimulus descriptions both as numbers and strings.
        
        Parameters
        __________
        
        beg : float
            First condition parameter/description.  Presumes a linear series for
            list of conditions.  If not true, use `list` parameter to provide 
            conditions list manually
            
        intvl : float
            Increment by which condition series increases from `beg`.  Increment 
            will be iterated as many times as there are conditions (derived from 
            `self.parameters`)
            
        con_type : string
            Condition type.  Stimulus parameter varying in experiment.  
            Select either - `orientation` - `spat_freq` - `temporal freq`.  
            Default is `orientation`.
            
        stim : string
            Stimulus type.  
            Select either - `bar` - `grating`.
            Default is `bar`.
            
        biphasic : Boolean
            Whether bar is biphasic, that is, whether it reverses direction in the
            one trial.
            
        unit : string
            Units of the condition parameters.  Default is `deg`.
            
        con_list : list (of strings)
            List of condition descriptions.  Must be strings convertible to floats,
            and provide one for each condition.
            
        temp_freq : float
            temporal frequency of drifting grating.  Required if `stim = 'grating'`.
            
        
        Returns
        _______
        
        self.conditions : array (1-D)
            Condition descriptions in numbers
            
        self.conditions2 : array (1-D)
            For when `biphasic == True`.  Lists the secondary condition descriptions for the
            returning bar.
            
        self.cond_label : list (of strings)
            String descriptions of conditions, chiefly for plotting purposes.
        """
        
        
        con_types = ['orientation', 'spat_freq', 'temporal_freq']
        stims = ['bar', 'grating']
        
        assert con_type.lower() in con_types, 'con_type invalid.  select from %s'%con_types
        assert stim.lower() in stims, 'stim invalid.  select from %s' %stims
        
        n_con = self.parameters['conditions']
        
        self.parameters['condition_type'] = con_type.lower()
        self.parameters['condition_unit'] = unit.capitalize()
        self.parameters['stimulus'] = stim.lower()
        
        if stim.lower() == stims[1]:
            # Gratings are GENERALLY not biphasic
            self.parameters['biphasic'] = 'N/A'
        else:
            self.parameters['biphasic'] = biphasic
        
        if stim == 'grating':
            self.parameters['temp_freq'] = float(temp_freq)
            
            # Sample rate must be a multiple of F1/temp_freq for it to be a frequency measured
            # in the FFT.
            samp_rate = 1/float(self.bin_width)
            
            
            assert samp_rate % temp_freq == 0., ('Bin_width (%s) is incompatible wih obtaining' 
                                                 'an FFT containing the specified temp_freq (%s). '
                                                 'The sampling frequency (1/bin_width) must be a'
                                                 'multiple of the temp_freq. \n\n Try as a' 
                                                 'bin_width %s and rerun self._sort().'
                                                 % (self.bin_width, temp_freq, 
                                                    1/(np.ceil(samp_rate/float(temp_freq))*temp_freq)))
        
        self.cond_label = []

        
        def circ(ori):
            """Ensures all orientation values are between 0 and 360 degrees.
            """
            ori[ori<-360] += 720
            ori[ori<0] += 360
            ori[ori>360] -= 360
            ori[ori>720] -= 720
            return ori

        # if list of conditions provided directly
        if len(con_list) > 0:
            
            assert len(con_list) == n_con, ('the number of labels provided '
                                        'manually (%s) does not match the '
                                        'number of conditions (%s).' % 
                                        (len(con_list), n_con))
                                        
            assert all(isinstance(l, str) for l in con_list), ('not all the '
                                                           'labels provided '
                                                           'are strings')
                                                                                                         
            
            self.cond_label = con_list
            
            # Relying on numpy conversion error should list be unable to convert to float.
            self.conditions = np.array(con_list).astype('f32')
                
        
        # if condition tpye is orientation
        elif con_type.lower() == con_types[0]:
            
            self.conditions = circ(np.arange(beg, beg+(n_con*intvl), intvl))
            
            assert len(self.conditions) == n_con, ('The amount of condition labels (%s) '
                                            'and conditions (%s) do not match; '
                                            'check your condition parameters' % 
                                            (self.cond_label.size, n_con))
            
            if biphasic:
                
                # self.conditions has been defined as an np.ndarray
                self.conditions2 = circ(self.conditions + 180) 

                for c in range(n_con):
                    label = '%s / %s %s' %(self.conditions[c], self.conditions2[c],
                                           self.parameters['condition_unit'])
                    self.cond_label.append(label)
                    
            else:
                
                for c in range(n_con):
                    label = '%s %s' %(self.conditions[c],
                                      self.parameters['condition_unit'])
                    self.cond_label.append(label)
                    
                    
        elif con_type.lower() == con_types[1]:
            self.conditions = np.arange(beg, beg + (n_con*intvl), intvl)
            
            assert len(self.conditions) == n_con, ('The amount of condition labels (%s) '
                                            'and conditions (%s) do not match; '
                                            'check your condition parameters' % 
                                            (self.cond_label.size, n_con))

            for c in range(n_con):
                label = '%s %s' %(self.conditions[c], self.parameters['condition_unit'])
                self.cond_label.append(label)

                    
        # if condition type is not predefined in this method            
        elif not con_type.lower() in con_types:
            
            self.conditions = np.arange(beg, beg+(n_con*intvl), intvl)
            
            for c in range(n_con):
                
                label = '%s %s' %(self.conditions[c],
                                  self.parameters['condition_unit'])
                self.cond_label.append(label)




    def _analyse(self, source='sdf', alpha = 0.05, n_bootstrap = 2000, 
                frequency=True, figsize = (15, 8)):
                    
        
        """
        
        Parameters
        __________
        
        alpha : float
            For bootstrap confidence intervals generated for the fourier analysis.
        
        n_bootstrap : int
            For bootstrap confidence intervals generated for the fourier analysis.
        
        """
        ## Add parameters to parameters dictionary
        
        # Organising source selection - raw and mov_avg not develoepd fully yet.
        sources = {'sdf': (self.spike_dens_func, self.CI_pos, self.CI_neg), 
                   'mov_avg': 'doesnt exist yet, call it self.spike_mov_avg', 
                   'raw': (self.conditions_hist_mean, 
                           self.conditions_hist_mean + 2*self.conditions_hist_stderr, 
                           self.conditions_hist_mean - 2*self.conditions_hist_stderr)}
                   
        assert source.lower() in sources.keys(), ('Tuning source data "%s" is invalid '
                                                    'select one of %s' %(source, sources.keys()))      
       
## Need to expand this functionality to the mean and CI_pos and CI_neg.  Doing so for
       # raw and moving average is not a priority, using sdf and bootstrap is pretty good.
       # overall aim is to clean this function up to accomadte a number of tuning functinos
       # in a clear and easy to use fasion.
       
        # Assign sources.
        
        
        
        
        n_con = self.parameters['conditions']
        
        
        # Get response values
        
        # values for transient bar responses
        if self.parameters['stimulus'] == 'bar':
            
            resp, CI_pos, CI_neg = sources[source.lower()]
            
            
            if self.parameters['biphasic']:
        
                half = self.bins.size/2

                max_val_arg = (resp[:, :half].argmax(axis=1),
                               resp[:, half:].argmax(axis=1)+half)
                                    
                max_val = (resp[:, :half].max(axis=1),
                           resp[:, half:].max(axis=1))
                           
                               
                max_val_CI_neg = (CI_neg[np.arange(n_con), max_val_arg[0]],
                                  CI_neg[np.arange(n_con), max_val_arg[1]])
                                  
                max_val_CI_pos = (CI_pos[np.arange(n_con), max_val_arg[0]],
                                  CI_pos[np.arange(n_con), max_val_arg[1]])
                                  
                self.cond_tuning = np.vstack((np.hstack((self.conditions, 
                                                        self.conditions2)),
                                             np.hstack(max_val),
                                             np.hstack(max_val_CI_neg),
                                             np.hstack(max_val_CI_pos)))
                                             
                # Convert to Hertz - design choice is to keep all PSTH datasets as raw average spike
                # counts, with easy option of seeing frequency in the plotting, but converting to 
                # Hertz for all condition tuning data.
                self.cond_tuning[1,:] *= (1/self.bin_width)
                                             
                              
                # Column labels for pd.dataframe of tuning data
                ci_perc = (100 * (1 - self.parameters['sdf_alpha']))
                idx = ['Conditions (%s)'% self.parameters['condition_unit'], 
                       'Peak Resp', 'Pos CI (%s%%)'%ci_perc, 'Neg CI (%s%%)'%ci_perc]
                # transpose of tuning array to data frame object       
                self.cond_tuning_pd = pd.DataFrame(self.cond_tuning.transpose(), columns=idx)
        
        
            if not self.parameters['biphasic']:
            #non biphasic version of above

                max_val_arg = resp[:, :].argmax(axis=1)
                                    
                max_val = resp[:, :half].max(axis=1)
                           
                               
                max_val_CI_neg = CI_neg[np.arange(n_con), max_val_arg]
                                  
                max_val_CI_pos = CI_pos[np.arange(n_con), max_val_arg]
                                  
                self.cond_tuning = np.vstack((self.conditions,
                                             max_val,
                                             max_val_CI_neg,
                                             max_val_CI_pos))
                                             
                # Convert to Hertz - design choice is to keep all PSTH datasets as raw average spike
                # counts, with easy option of seeing frequency in the plotting, but converting to 
                # Hertz for all condition tuning data.
                self.cond_tuning[1,:] *= (1/self.bin_width)
                
                
                # Column labels for pd.dataframe of tuning data
                ci_perc = (100 * (1 - self.parameters['sdf_alpha']))
                idx = ['Conditions (%s)'% self.parameters['condition_unit'], 
                       'Peak Resp', 'Pos CI (%s%%)'%ci_perc, 'Neg CI (%s%%)'%ci_perc]
                       
                # transpose of tuning array to data frame object       
                self.cond_tuning_pd = pd.DataFrame(self.cond_tuning.transpose(), columns=idx)
        
        
        # values for sinusoids/gratings
        ## Note issue of temporal frequency tuning - need variable tf.
        if self.parameters['stimulus'] == 'grating':
            
            self.parameters['fft_alpha'] = alpha
            self.parameters['fft_number_bootstrap'] = n_bootstrap
            
            if source == 'sdf':
                print ('WARNING, using a smoothed/filtered dataset will artificially increase'
                       'the amplitude of the DC component and decrease that of the F1') 
            
            sources = {'sdf': self.conditions_trials_sdf,
                       'mov_avg': "doesn't exist yet (?)",
                       'raw': self.conditions_trials_hist}
            
            resp = sources[source]
            
            temp_freq = self.parameters['temp_freq']
            stim_len = self.parameters['stimulus_length']
            
            # ensuring that the temp_freq is measured in the FFT whilst taking the maximum time.
            # on the basis of delt-f = 1 / n*del-t; stim_len*F1=factor; 1/(bin_width*F1)=min bins
            # number times greater than minimum can fit in stim_length            
            factor = np.floor(stim_len * temp_freq).astype('int')
            
            # number of bins to take - the window size necessary for temp_freq to be measured
            bins_take = np.floor(factor / (self.bin_width * temp_freq)).astype('int')

            # Frequency axis generation
            self.freq = fft.rfftfreq(bins_take, self.bin_width)
            
            #Checkign whether the temp_freq is in the FFT.
            assert self.freq[factor] == temp_freq, ('The calculated FFT F1 frequency (%s)'
                                                       'does not equal the Stimulus temp_freq (%s)'
                                                       %(self.freq[bins_take], temp_freq))

            # Fourier Transform
            self.conditions_trials_fourier = fft.rfft(resp[:,:,:bins_take], axis=2)
            
            # Amplitude (peak-to-peak)
            self.conditions_trials_ampl = np.abs(self.conditions_trials_fourier)
            
            # normalising to dataset size, except the DC.
            self.conditions_trials_ampl[:,:,0] *= 1 / float(bins_take)
            self.conditions_trials_ampl[:,:,1:] *= 2 / float(bins_take)
            
            
            # Mean amplitudes and bootstrapped CI_intervals            
            self.conditions_ampl_mean = np.mean(self.conditions_trials_ampl, axis=1)
            
            CI_pos, CI_neg = bootstrap(self.conditions_trials_ampl, alpha=alpha, 
                                       n_bootstrap=n_bootstrap)
            self.conditions_ampl_CI_pos, self.conditions_ampl_CI_neg = CI_pos, CI_neg
            
            # isolating F0, F1, and F2 responses and compiling into a single table.
            conditions_f0 = self.conditions_ampl_mean[:,0]
            conditions_f1 = self.conditions_ampl_mean[:,factor]
            conditions_f2 = self.conditions_ampl_mean[:,2*factor]
            
            # Condition Tuning array
            self.cond_tuning = np.vstack((self.conditions,
                                         conditions_f0, CI_pos[:,0], CI_neg[:,0],
                                         conditions_f1, CI_pos[:,factor], CI_neg[:,factor],
                                         conditions_f2, CI_pos[:,2*factor], CI_neg[:,2*factor],
                                         conditions_f1/conditions_f0))
            
            # Convert to Hertz - design choice is to keep all PSTH datasets as raw average spike
            # counts, with easy option of seeing frequency in the plotting, but converting to 
            # Hertz for all condition tuning data.
            
            self.cond_tuning[1:-1,:] *= (1/self.bin_width)
            
            # Column labels for pd.dataframe of tuning data
            ci_perc = (100 * (1 - self.parameters['fft_alpha']))
            idx = ['Conditions (%s)'% self.parameters['condition_unit'], 
                   'F0', 'F0 Pos CI (%s%%)'%ci_perc, 'F0 Neg CI (%s%%)'%ci_perc, 
                   'F1', 'F1 Pos CI (%s%%)'%ci_perc, 'F1 Neg CI (%s%%)'%ci_perc,
                   'F2', 'F2 Pos CI (%s%%)'%ci_perc, 'F2 Neg CI (%s%%)'%ci_perc,
                   'F1/F0 Ratio']
            # transpose of tuning array to data frame object       
            self.cond_tuning_pd = pd.DataFrame(self.cond_tuning.transpose(), columns=idx)
        
        
        
############        
        # for orientation data, the orientation angles can get scrambled due to the circ() function
        # rotating the angles around.  This orders them numerically in the final cond_tuning
        
        if self.parameters['condition_type'] == 'orientation':
            self.cond_tuning = self.cond_tuning[:,self.cond_tuning[0].argsort()]
            self.cond_tuning_pd.sort(self.cond_tuning_pd.columns[0], inplace=True)

        
        
    def _out(self):
        directory = self.parameters['directory']
        data_label = self.parameters['data_label']
        con_type = self.parameters['condition_type']
        
        param = pd.Series(self.parameters)
        with open('%s%s_parameters.csv'%(directory, data_label), 'w') as f:
            param.to_csv(f)
        
        
        with open('%s%s_%s_tuning.csv'%(directory, data_label, con_type), 'w') as f:
            self.cond_tuning_pd.to_csv(f)
            
        ## Other files to save etc
            # PSTH trial and means and SDF - pd with multiindexing and excel write?
            # make convenience functions for reading particular files?
            
    def _plot_psth(self, plot_type = 'hist', frequency = True, 
                   smooth = False, smooth_type = 'gauss', sigma = 3, 
                   mov_avg_window = 3, figsize=(15, 8), format=True):

        """Plots the PSTHs (histograms) for each of the conditions
        
        Parameters
        __________
        plot_type : string
            Type of plot for hte PSTHs.  Either 'hist', 'line', 'sdf'.
            
        frequency : Boolean
            Convert spike counts to spike frequency in Hz
            
        smooth : Boolean
            Plot, on top, a spike density curve derived from gaussian kernel convolution of
            the PSTH bins.
            Applies to both 'hist' and 'line' plot types.
            
        smooth_type : string
            type of smoothing to employ over the PSTH data.
            'gauss' employs a gaussian kernel of `sigma` standard deviation (in bins)
            'mov_avg' employs a wimple square kernal of `sigma` width (in bins)
            
        sigma : float
            If `smooth_type == "gauss"`: standard deviation of the gaussian kernel 
            used to produce the smoothed curve.
            If `smooth_type == "mov_avg"`: width of the square kernel.
            This parameter is dimesionless.  It refers to the number of bins.  The standard
            deviation or window width in units of time can be calculated from 
            `sigma * self.bin_width`.
            
        mov_avg : Boolean
            Not developed yet.  Alternative to density curve is a simple moving average.
            
        
        Returns
        _______
        Plots : function
            Passes an `ax.plot()` function call.
        
        
        """
        plot_types = ['hist', 'line', 'sdf']
        assert plot_type.lower() in plot_types, ('plot_type invalid,\n select from %s' %plot_types)
        
        assert smooth_type.lower() in ['gauss', 'mov_avg'], ('Smoothing type invalid, select '
                                                            'either "gauss" or "mov_avg"')
        # Generate smoothed curves                                        
        if smooth:
            if smooth_type.lower() == 'gauss':
                # Gaussian smoothing 
                spd = gaussian_filter1d(self.conditions_hist_mean, sigma, axis=1)
                
            if smooth_type.lower() == 'mov_avg':
                # Moving average smoothing
                spd = sp.ndimage.convolve1d(self.conditions_hist_mean,
                                            np.ones((sigma))/float(sigma), axis=1) 

            

        n_con = self.parameters['conditions']
        bin_width = self.bin_width
        
        cols = 3
        rows = math.trunc(n_con/3.) + (n_con%3.)        
        
                
        fig = plt.figure(figsize = figsize)

                 
        if plot_type.lower() == plot_types[0]:
            
            # stepping the dataset for a fill_between plot by stacking and 
            # flattening to double up the values and bin boundaries
            # this is substantially faster than the bar plot for dataset this
            # large
            
            step_bins = np.column_stack((self.bins, self.bins)).flatten()[1:-1]
            step_hist_mean = np.column_stack((self.conditions_hist_mean.flatten(),
                                              self.conditions_hist_mean.flatten()))
            step_hist_mean = step_hist_mean.flatten()
            step_hist_mean = step_hist_mean.reshape((n_con, (self.bins.size-1)*2))
            
            for c in range(n_con):
                ax = fig.add_subplot(rows, cols, c+1)
                
                ax.set_ylim(0, self.conditions_hist_mean.max() * 1.14)       
                ax.set_xlim(0, self.bins[-1])
                
                ax.fill_between(step_bins, step_hist_mean[c], lw=0, 
                                facecolor='0.3', zorder=2)
                                
                ax.set_title(self.cond_label[c])
                
                

                if smooth:
                    ax.plot(self.bins[:-1] + 0.5*bin_width, spd[c], 
                            linestyle='-', color='FireBrick', linewidth=4, alpha=0.8, 
                            zorder=3)
                    
            # convert ylabl to frequency units
                    
                if frequency:
                    freq_label = np.round(ax.get_yticks() * (1 / bin_width), 
                                          decimals=1)
                    ax.set_yticklabels( freq_label)
                
                    # ylabel for all left most subplots
                for sub_plt in np.arange(1, rows*cols, cols):
                    if sub_plt == (c+1):
                        if frequency:
                            ax.set_ylabel('Frequency')
                        else:
                            ax.set_ylabel('Average count')
                
#                plotform(ax)



        elif plot_type.lower() == plot_types[1]:
            
            for c in range(n_con):
                
                pos_stderr = self.conditions_hist_mean + 2*self.conditions_hist_stderr
                neg_stderr = self.conditions_hist_mean - 2*self.conditions_hist_stderr
    
                ax = fig.add_subplot(rows, cols, c+1)
                
                ax.set_ylim(0, pos_stderr.max() * 1.14)       
                ax.set_xlim(0, self.bins[-1])            
                
                
                ax.plot(self.bins[:-1] + 0.5*bin_width, self.conditions_hist_mean[c],
                        color='0.28', linewidth=1)
    
                ax.fill_between(self.bins[:-1] + 0.5*bin_width, 
                                pos_stderr[c], 
                                neg_stderr[c],
                                color='0.6', alpha=0.6)
                                
                ax.set_title(self.cond_label[c])
                
    
    
                if smooth:
                    ax.plot(self.bins[:-1] + 0.5*bin_width, spd[c], 
                            linestyle='-', color='FireBrick', linewidth=3, alpha=0.8, 
                            zorder=3)
    
    
                    
                # convert ylabl to frequency units
                        
                if frequency:
                    freq_label = np.round(ax.get_yticks() * (1 / bin_width), 
                                          decimals=1)
                    ax.set_yticklabels( freq_label)
                
                    # ylabel for all left most subplots
                for sub_plt in np.arange(1, rows*cols, cols):
                    if sub_plt == (c+1):
                        if frequency:
                            ax.set_ylabel('Frequency')
                        else:
                            ax.set_ylabel('Average count')
                

#                plotform(ax)
                
                
                
        elif plot_type.lower() == plot_types[2]:
            
            cols = 3
            rows = math.trunc(n_con/3.) + (n_con%3.)        
            
            
            for c in range(n_con):
                ax = fig.add_subplot(rows, cols, c+1)
                
                ax.set_ylim(0, self.CI_pos.max() * 1.14)       
                ax.set_xlim(0, self.bins[-1])            
                
                
                ax.plot(self.bins[:-1] + 0.5*bin_width, self.spike_dens_func[c],
                        lw=2, color='#036eb6', zorder=1)
    
                ax.fill_between(self.bins[:-1] + 0.5*bin_width, 
                                self.CI_neg[c], self.CI_pos[c],
                                color='#036eb6', alpha=0.3)
                                
                ax.set_title(self.cond_label[c])
                                
                if frequency:
                    freq_label = np.round(ax.get_yticks() * (1 / bin_width), 
                                              decimals=1)
                    ax.set_yticklabels( freq_label)
                    
                for sub_plt in np.arange(1, rows*cols, cols):
                    if sub_plt == (c+1):
                        if frequency:
                            ax.set_ylabel('Frequency')
                        else:
                            ax.set_ylabel('Average count')
    
          # bug with this and teh macosx backend      
#        plt.tight_layout()
        plt.subplots_adjust(hspace=0.45)

    def _plot_psth_flat(self, sigma=5, figsize = (15, 8)):
        """ Do all conditions side by side with a small filter
        """
    
        gaus_filt = sp.ndimage.gaussian_filter1d
        all_resp = gaus_filt(self.conditions_hist_mean.flatten(), sigma)
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
        
        ax.plot(all_resp, linestyle='-', color='0.28')
        
        n_con = self.parameters['conditions']
        con_mark = np.arange(0, (self.bins.size -1) * n_con, self.bins.size -1)
                
        ax.xaxis.set_ticks(con_mark)
        ax.xaxis.set_ticklabels(self.cond_label)
        
        freq_label = np.round(ax.get_yticks() * (1/self.bin_width),
                              decimals = 1)
        ax.set_yticklabels(freq_label)
        ax.set_ylabel('Frequency')
        
        for label in ax.xaxis.get_majorticklabels():
            label.set_horizontalalignment('left')
            
        ax.set_xlim(0, (self.bins.size -1) * n_con)
        
        # bug with macosx backend
#        plt.tight_layout()
        plt.subplots_adjust(hspace=0.45)
        
            
    
    def _plot_tuning(self, plot_type = 'linear', frequency=True, modulation='both',
                     xaxis_scale = 'linear', xaxis_base = '2',
                     figsize = (10, 6)):
        
        """
        
        Parameters
        __________
        
        modulation : string
            For grating stimuli, which modulation to plot.
            Either 'F1', 'F0' or 'both'.  Default is 'both'.
            
        
        """
        
        def wrap(array):
            """ 
            returns wrapped version of array, with the first element
            in each rwo appended to the end of the same rows.
            
            For polar plotting or angles/orientations
            """
            
            assert array.ndim == 2, ('works on 2D arrays only; expecting '
                                     'self.cond_tuning, which should be 2D')
                                            
            return np.column_stack((array, array[..., 0]))
        

        assert plot_type.lower() in ['polar', 'linear'], ('plot type unrecognised.  Use either '
                                                        ' "polar" or "linear"')
                                                        
        assert xaxis_scale.lower() in ['linear', 'log'], ('xaxis_scale must be either "linear" or "log"')                                                

        if plot_type.lower() == 'polar':
            assert self.parameters['condition_type'] == 'orientation', ('Polar plotting makes '
                                                                        'sense only for orientation '
                                                                        'tuning data.')
                                                                        
            assert xaxis_scale.lower() == 'linear', ('Log scale for a polar plot ...?')

                                                            
        
        xaxis_scale_kw = {'basex':float(xaxis_base)}
        
        if self.parameters['stimulus'] == 'grating':
            
            if modulation.lower() == 'f0':
                data = self.cond_tuning[(0, 1, 2, 3),:]
            
            if modulation.lower() == 'f1':
                data = self.cond_tuning[(0, 4, 5, 6),:]

            if modulation.lower() == 'both':
                data = self.cond_tuning[:7,:]
                
        if self.parameters['stimulus'] == 'bar':
            data = self.cond_tuning


                
        if self.parameters['condition_type'] == 'orientation':
        
            if plot_type.lower() == 'polar':
                
                fig = plt.figure(figsize=figsize)                
                ax=fig.add_subplot(111, polar=True)
                
                data = wrap(data)
                
                ax.plot(np.radians(data[0]), data[1], color='#507df5', lw=2)
                
                ax.fill_between(np.radians(data[0]), data[2], data[3], 
                                           color='#507df5', alpha=0.5, lw=0)
                                                              
            
            elif plot_type.lower() == 'linear':
                
                fig = plt.figure(figsize=figsize)
                ax = fig.add_subplot(111)
                                
                
                ax.plot(data[0], data[1], color='#507df5', lw=2)
                ax.fill_between(data[0], data[2], data[3], color='#507df5',
                                alpha=0.5, lw=0)
                                
                ax.set_xlabel(self.parameters['condition_type'].capitalize()+' '+self.parameters['condition_unit'])
                    

                
        else:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
                        
            
            if data.shape[0] == 7:
                ax.plot(data[0], data[1], color='#507df5', lw=2, label='F0')
                ax.fill_between(data[0], data[2], data[3], color='#507df5',
                                alpha=0.5, lw=0)
                                
                ax.plot(data[0], data[4], color='#bc5692', lw=2, label='F1')
                ax.fill_between(data[0], data[5], data[6], color='#bc5692',
                                alpha=0.5, lw=0)
                ax.legend()
            else:
                ax.plot(data[0], data[1], color='#507df5', lw=2)
                ax.fill_between(data[0], data[2], data[3], color='#507df5',
                                alpha=0.5, lw=0)
                            
            ax.set_xlabel(self.parameters['condition_type'].capitalize()+' '+self.parameters['condition_unit'])
            ax.set_xscale(xaxis_scale, **xaxis_scale_kw)
            
            
        
        ax.set_ylabel('Response (Hz)')
            
            
                
## spontaneous rate
                # plotting - subtract or line?
                # handling a blank condition - presume always last?
                # Three ways - blank condition, or blank portion of stimulus, or, pre/post-conditions

## Post-Marker Buffer
    # Unsure of implementation.  Perhaps at analysis, just take only a certain portion of PSTH.
    # Good to keep all data that we can.                
                
## Multiple sources for tuning or just sdf?  If multiple, need cleaner means of management.
#                Options include, raw mean, moving average.  All should be dealt with in the same way.
                
## Saving
                

## in relation to Pandas integration ... general I/O.

## convenience function for listing all attributes

## Conditions / Trials load (for when marker issue doesn't exist)
## Pandas Integration -> question of what to provide a pandas form for
# key advantage will be outputting to excel or csv easily
# possibly provide functions for viewing and excel purposes that provide
# returns and not attributes

## ATHENA!!! (analysis!!)
# Orientation Tuning
# Spike Shape (simple), maybe a zeus function

np.convolve()