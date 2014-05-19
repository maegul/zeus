# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 16:06:22 2014

@author: errollloyd
"""
import math
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.ndimage import gaussian_filter1d



import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

import os
import glob

def load(dir):
    
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

    
    elif len(data.keys()) == len(files):
        print 'All files imported and assigned'
        print '\nFollowing files loaded successfully:\n'
        for i in file_names: print(i)
        return data
        
    
        
        


def plotform(ax, tickint=False):
    


    # remove top and right axes
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    # make ticks outward
    ax.tick_params('both', direction='out')
    
    
#    fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
#    'weight' : 'normal', 'size' : 12}
    
#    if tickint:
#        p.set_xticklabels(p.get_xticks().astype('int'), fontProperties)
#        p.set_yticklabels(p.get_yticks().astype('int'), fontProperties)
#    else:
#        p.set_xticklabels(p.get_xticks(), fontProperties)
#        p.set_yticklabels(p.get_yticks(), fontProperties)
        
        
    #set the style of the major and minor grid lines, filled blocks
    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
    ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
    ax.patch.set_facecolor('0.85')
    ax.set_axisbelow(True)
    
    #set minor tick spacing to 1/2 of the major ticks
    ax.xaxis.set_minor_locator(MultipleLocator( (plt.xticks()[0][1]-plt.xticks()[0][0]) / 2.0 ))
    ax.yaxis.set_minor_locator(MultipleLocator( (plt.yticks()[0][1]-plt.yticks()[0][0]) / 2.0 ))
    
    #remove axis border
#    for child in ax.get_children():
#        if isinstance(child, matplotlib.spines.Spine):
#            child.set_alpha(0)
       
    #restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)
    
    #remove the minor tick lines    
    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)

    
    
    if ax.legend_ <> None:
        lg = ax.legend_
        lg.get_frame().set_linewidth(0)
        lg.get_frame().set_alpha(0.5)
        
    plt.tight_layout()


def steppify(arr, axis='x'):
    """
    expands datasets to create de factor step (step-post) plots
    """
    
    if axis == 'x':
        newarr = np.r_[arr[0], np.dstack((arr[1:], arr[1:])).flatten()]
    
    elif axis == 'y':
        newarr = np.r_[np.dstack((arr[:-1], arr[:-1])).flatten(), arr[-1]]
    
    else:
        print 'your axes in steppify are improperly identified'

    return newarr



## Class Definitions

class Zeus:
    
    
    
    def __init__(self, dir):
        
        # Load files and assign to attributes
        self.data = load(dir)
        
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

        self.parameters = {}
        self.parameters['directory'] = dir
        
        
    def _sort(self, conditions = 9, trials = 10):
        
        
        if self.markers.size != (conditions * trials):
            print 'The number of recorded markers does not match your data specifications'
            exit()
        
        self.parameters['conditions'] = conditions
        self.parameters['trials'] = trials        
        
        # this should be a parameter while bin_width should be an attribute
        # as it is arbitrary to the analysis
        stim_len = np.max( np.subtract( self.markers[1:], self.markers[:-1]))
        self.parameters['stimulus_length'] = stim_len
        
        
        self.conditions_trials = {}
        self.conditions_trials_zeroed = {}
        
        
        for c in range(conditions):
            
            ct = [] # temp for condition_trials
            ctz = [] # temp for condition_trials_zeroed
            
            for t in self.markers[c::conditions]:
                
                trial = np.extract( (self.spikes >= t) & (self.spikes < (t+stim_len)), 
                                   self.spikes)
                
                ct.append(trial)
                ctz.append( (trial - t))
                
                
            self.conditions_trials[c] = ct
            self.conditions_trials_zeroed[c] = ctz
        
    
    
    def _analyse(self, bin_width=0.02):
        
        self.bin_width = bin_width
        
        stim_len = self.parameters['stimulus_length']
        
        self.bins = np.arange(0, stim_len, bin_width)
        
        n_con = self.parameters['conditions']
        n_trial = self.parameters['trials']
        
        self.conditions_trials_hist = np.zeros((n_con, n_trial, self.bins.size - 1))
        
        
        for c in range(n_con):
            
            trials_hist = np.zeros((n_trial, self.bins.size - 1))
            
            for t in range(n_trial):
                
                trials_hist[t] = np.histogram(self.conditions_trials_zeroed[c][t], self.bins)[0]
                
            self.conditions_trials_hist[c] = trials_hist
            
            
        self.conditions_hist_mean = np.mean(self.conditions_trials_hist, axis=1)
        self.conditions_hist_stderr = stats.sem(self.conditions_trials_hist, axis= 1, ddof=0)
     


    def _conditions(self, beg=-90, intvl=20, type='orientation', stim='bar', 
                    biphasic=True, unit='deg', input='bounds', list=[]):
                        
        con_types = ['orientation', 'spat_freq', 'temporal_freq']
        
        n_con = self.parameters['conditions']
        
        self.parameters['condition_type'] = type
        self.parameters['condition_unit'] = unit.capitalize()
        self.parameters['stimulus'] = stim
        self.parameters['biphasic'] = biphasic
        
        self.con_label = []

        
        def circ(ori):
            ori[ori<0] += 360
            ori[ori>360] -= 360
            ori[ori>720] -= 720
            return ori

        
        if type.lower() == con_types[0]:
            
            if biphasic:
                
                self.con = circ(np.arange(beg, beg+(n_con*intvl), intvl))
                self.con2 = circ(self.con + 180)

                for c in range(n_con):
                    label = '%s / %s %s' %(self.con[c], self.con2[c],
                                           self.parameters['condition_unit'])
                    self.con_label.append(label)
                    
            else:
                self.con = circ(np.arange(beg, beg+(n_con*intvl), intvl))
                
                for c in range(n_con):
                    label = '%s %s' %(self.con[c],
                                      self.parameters['condition_unit'])
                    self.con_label.append(label)
                    
                    
        if not type.lower() in con_types:
            
            self.con = np.arange(beg, beg+(n_con*intvl), intvl)
            
            for c in range(n_con):
                
                label = '%s %s' %(self.con[c],
                                  self.parameters['condition_unit'])
                self.con_label.append(label)

        
        
        

    def _psth(self, figsize=(15, 8), format=True, plot_type = 'hist', 
             frequency = True, density = False, sigma = 3, mov_avg = False,
             mov_avg_window = 3):


        n_con = self.parameters['conditions']
        bin_width = self.bin_width
        
        cols = 3
        rows = math.trunc(n_con/3.) + (n_con%3.)        
        
                
        fig = plt.figure(figsize = figsize)

                 
        if plot_type.lower() == 'hist':
            
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
                                
                ax.set_title(self.con_label[c])
                
                

                if density:
                    
                    ### add simply moving average (which can have std err too)
                    
                    # Density function ... should be done across the whole
                    # data set in .Analyse()
                           
                    spd = gaussian_filter1d(self.conditions_hist_mean[c], sigma)
                    
                    
                    ax.plot(self.bins[:-1] + 0.5*bin_width, spd, 
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
                
                plotform(ax)



        elif plot_type.lower() == 'line':
            
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
                                
                ax.set_title(self.con_label[c])
                
    
    
                if density:
                    
                    ### add simply moving average (which can have std err too)
                    
                    # Density function ... should be done across the whole
                    # data set in .Analyse()
                           
                    spd = gaussian_filter1d(self.conditions_hist_mean[c], sigma)
                    
                    
                    ax.plot(self.bins[:-1] + 0.5*bin_width, spd, 
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
                

                plotform(ax)

    def _psth_flat(self, sigma=5, figsize = (15, 8)):
        # Do all conditions side by side with a small filter
    
        gaus_filt = sp.ndimage.gaussian_filter1d
        all_resp = gaus_filt(self.conditions_hist_mean.flatten(), sigma)
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
        
        ax.plot(all_resp, linestyle='-', color='0.28')
        
        n_con = self.parameters['conditions']
        con_mark = np.arange(0, (self.bins.size -1) * n_con, self.bins.size -1)
                
        ax.xaxis.set_ticks(con_mark)
        ax.xaxis.set_ticklabels(self.con_label)
        
        for label in ax.xaxis.get_majorticklabels():
            label.set_horizontalalignment('left')
            
        ax.set_xlim(0, (self.bins.size -1) * n_con)
        
        plotform(ax)
        
        
    def _sdf(self, frequency = True, conf_int = True, conf_int_type = 'boostrap', 
            sigma=3, alpha=0.05, n_bootstrap=2000, plot=True, figsize=(15, 8)):
        """
        Generate a Spike Density Function (SDF) for the data set using 
        a smoothened form of the PSTH and confidence intervals derived from
        bootstrapping
        """
        
        n_trials = self.parameters['trials']
        n_con = self.parameters['conditions']
        bin_width = self.bin_width
        
        # make the kernel 5 sigmas wide in each direction
        kernel = stats.norm.pdf(np.arange(-5*sigma, (5*sigma)+1), scale=sigma)
        
        bootstrap_index = np.random.randint(0, n_trials, 
                                            (n_trials, n_bootstrap) )
        
        # For each bin in the histogram, randomly samples from the results
        # of each trial and repeats, effectively, n_bootstrap times                
        trials_bootstrap = self.conditions_trials_hist[:, bootstrap_index, :]
        
        # dimension one is the trials, zero is the conditions; this averaging 
        # goes across the trials creating a PSTH for each condition, and,
        # importantly, for each bootstrap resample
        PSTH_bootstrap = trials_bootstrap.mean(axis=1)
        
        # smoothing along the bins (ie, a SDF), where dimension zero is the
        # conditions, one the bootstrap resamples and two the bins
        trials_bootstrap_filt = sp.ndimage.convolve1d(PSTH_bootstrap, kernel, 
                                                      axis = 2)
        
        # find percentile values for each bin along the bootstrap resamples,
        # which are on axis 1                                              
        self.CI_pos = np.percentile(trials_bootstrap_filt, 100*(1 - (alpha/2.)), 
                               axis=1)
        self.CI_neg = np.percentile(trials_bootstrap_filt, 100*(alpha/2.), 
                               axis=1)
                               
        # the spike density function (smoothed PSTH)                       
        self.spike_dens_func = sp.ndimage.convolve1d(self.conditions_hist_mean, 
                                                     kernel, axis=1)
                               
        # Plotting
        if plot:
            cols = 3
            rows = math.trunc(n_con/3.) + (n_con%3.)        
            
            fig = plt.figure(figsize = figsize)
            
            for c in range(n_con):
                ax = fig.add_subplot(rows, cols, c+1)
                
                ax.set_ylim(0, self.CI_pos.max() * 1.14)       
                ax.set_xlim(0, self.bins[-1])            
                
                
                ax.plot(self.bins[:-1] + 0.5*bin_width, self.spike_dens_func[c],
                        lw=2, color='#036eb6', zorder=1)
    
                ax.fill_between(self.bins[:-1] + 0.5*bin_width, 
                                self.CI_neg[c], self.CI_pos[c],
                                color='#036eb6', alpha=0.3)
                                
                ax.set_title(self.con_label[c])
                                
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
    
                
                plotform(ax)
            
    
    def _tuning(self, ori_linear=False, source='sdf'):
        
        n_con = self.parameters['conditions']
        
        if source.lower() == 'sdf':
            if self.parameters['condition_type'] == 'orientation':
                if self.parameters['biphasic']:
                    
                    half = self.bins.size/2

                    max_val_arg = (self.spike_dens_func[:, :half].argmax(axis=1),
                                   self.spike_dens_func[:, half:].argmax(axis=1)+half)
                                        
                    max_val = (self.spike_dens_func[:, :half].max(axis=1),
                               self.spike_dens_func[:, half:].max(axis=1))
                               
                                   
                    max_val_CI_neg = (self.CI_neg[np.arange(n_con), max_val_arg[0]],
                                      self.CI_neg[np.arange(n_con), max_val_arg[1]])
                                      
                    max_val_CI_pos = (self.CI_pos[np.arange(n_con), max_val_arg[0]],
                                      self.CI_pos[np.arange(n_con), max_val_arg[1]])
                                      
                    self.con_tuning = np.vstack((np.hstack((self.con, self.con2)),
                                                 np.hstack(max_val),
                                                 np.hstack(max_val_CI_neg),
                                                 np.hstack(max_val_CI_pos)))
                                  
                    self.con_tuning = self.con_tuning[:,self.con_tuning[0].argsort()]
            
                
                
        
        