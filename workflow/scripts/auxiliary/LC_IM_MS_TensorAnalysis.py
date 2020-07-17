"""
This is the core module of the "Hydrogen-Deuterium-eXchange Liquid-chromatography Ion-mobility-separation Mass-spectrometry Interactive Tensor analysis" (HDX-LIMIT) pipeline.

HDX-LIMIT defines classes which automate processing of a timeseries of 3D LC-IMS-MS data extracted from .mzML.gz files.


The module contains 3 Data-Classes and 2 Processing-Classes:
    
    Data-Classes: 
        - DataTensor: An object created from LC-IM-MS data extracted from a .mzML.gz file, with methods to perform a deconvoluting Tensor-Factorization.
        - Factor: Objects created and stored within the DataTensor through a Tensor-Factorization - the n-dimensional analog to multaplicative-factorization.
        - IsotopicCluster: Objects within a factor identified as 'looking sufficiently like' an MS protein signal.
    
    Processing-Classes:
        - TensorGenerator: Main class for creating tensors from incoming data, handles all timeseries and charge state data. Stores all DataTensor outputs internally. 
        - PathOptimizer: Accepts all IsotopicClusters from TensorGenerator, generates time-series of cluseters by score-based optimization from bootstrap of 'series-space'.



The TensorGenerator and PathOptimizer classes are instantiated once per 'protein retention-time group' (RT-group):

    In the HDX-LIMIT approach, a 'protein RT-group' is a collection of isotopic-clusters identified as close to the expected mass of a library protein-of-interest, for all 
    non-redundant charge states that fall within a parameterized range of LC retention-time. Depending on the width of a protein's elution window, multiple RT-groups
    can be identified for a single sample protein - and are handled separately to increase the stability of the tensor factorization.

A TensorGenerator creates a DataTensor for all charges of an RT-group, for all HDX times. The TensorGenerator iterates over the timeseries, instantiates all DataTensors,
factorizes them, and saves their resulting Factors and IsotopicClusters internally as attributes, discarding the DataTensor objects at each timepoint.

The IsotopicClusters are passed to an instance of PathOptimizer, which creates a collection of 'sample-paths' - set of bootstrapped m/z centroid time-series created by 
analytically generating a series of trajectories with ~200 different unfolding rates. These sample-paths are populated with the IsotopicClusters from the TensorGenerator 
which are closest to the expected COM of that path at each timepoint. Each path is optimized to a cost-function by finding the single-best-substitution (the cluster at 
any timepoint which will best-minimze the cost function), and this is repeated until no substitutions yeild benefits. The series with the lowest-score is kept as that 
most-expected to contain the correct IsotopeClusters at each timepoint, and representing the unfolding dynamics of the protein. These series are then fed downstream to 
scripts which determine the rate of unfolding for each residue and the overall DeltaG_unfolding of the protein. 
"""



##########################################################################################################################################################################################################################################
### IMPORTS ##############################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################



### LIBRARIES ###
import re
import os
import sys
import glob
import math
import zlib
import copy
import time
import ipdb
import pymzml
import pickle
import tensorly
import peakutils
import statistics
import importlib.util

### ALIASES ###
import scipy as sp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import _pickle as cpickle
import matplotlib.pyplot as plt

### FUNCTIONS ###
from Bio.PDB import *
from scipy.stats import gmean
from scipy.fftpack import fft
from scipy.signal import argrelmax
from Bio.SeqUtils import ProtParam
from IPython.display import display
from mpl_toolkits.mplot3d import Axes3D
from tensorly.decomposition import non_negative_parafac
from scipy.spatial.distance import pdist, squareform
from scipy.ndimage.filters import gaussian_filter

### BOKEH ###
from bokeh.plotting import figure
from bokeh.palettes import Spectral6
from bokeh.transform import linear_cmap
from bokeh.layouts import gridplot, column
from bokeh.models import HoverTool, ColorBar, Text, Div, Whisker
from bokeh.models.glyphs import MultiLine, Line
from bokeh.io import save, output_file
from bokeh.models.callbacks import CustomJS
from bokeh.models.sources import ColumnDataSource, CDSView
from bokeh.models.filters import Filter, GroupFilter, IndexFilter

### SPECLOAD ###
hxtools_spec = importlib.util.spec_from_file_location("hxtools", "scripts/auxiliary/hxtools.py")
hxtools = importlib.util.module_from_spec(hxtools_spec)
hxtools_spec.loader.exec_module(hxtools)



##########################################################################################################################################################################################################################################
### Classes ##############################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################



###
### Class - DataTensor:
### Holds the data for a single identified cluster, generated from subspace_selector.ipynb
###

class DataTensor:
    
    def __init__(
        self, 
        source_file,
        tensor_idx, 
        timepoint_idx, 
        name, 
        total_mass_window,
        n_concatenated,
        charge_states,
        **kwargs
        ):

        
        ###Set Common Attributes###

        #Positional args
        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.name = name
        self.total_mass_window = total_mass_window
        self.n_concatenated = n_concatenated
        self.charge_states = charge_states

        #Keyword Args
        if kwargs is not None:
            kws = list(kwargs.keys())
            if 'rts' in kws:
                self.rts = np.array(kwargs['rts'])
            if 'dts' in kws: 
                self.dts = np.array(kwargs['dts'])
            if 'seq_out' in kws:
                self.seq_out = np.array(kwargs['seq_out'])
            if 'int_seq_out' in kws and kwargs['int_seq_out'] is not None: 
                self.int_seq_out = np.array(kwargs['int_seq_out'])
                self.int_seq_out_float = self.int_seq_out.astype('float64')
                self.int_grid_out = np.reshape(self.int_seq_out_float, (len(self.rts), len(self.dts), 50))
                self.int_gauss_grids = self.gauss(self.int_grid_out)
            if 'concat_dt_idxs' in kws:
                self.concat_dt_idxs = kwargs['concat_dt_idxs']
            if 'concatenated_grid' in kws:
                self.concatenated_grid = kwargs['concatenated_grid']
            if 'lows' in kws:
                self.lows = kwargs['lows']
            if 'highs' in kws:
                self.highs = kwargs['highs']
            if 'abs_mz_low' in kws:
                self.mz_bin_low = kwargs['abs_mz_low']

        ###Compute Instance Values###
        
        #Handle normal case: new DataTensor from output of isolate_tensors.py
        if self.n_concatenated == 1:

            self.grid_out = np.reshape(self.seq_out, (len(self.rts), len(self.dts)))

            #For creating full_grid_out, a 3d array with all mz dimensions having the same length
            #Find min and max mz range values:
            #check head and tail of each grid_out index ndarray and compare them
            self.mz_bin_high = 0
            self.mz_bin_low = 1E9
            for j in range(len(self.grid_out)):
                for k in range(len(self.grid_out[j])):
                    if np.shape(self.grid_out[j][k])[0] != 0:
                        if self.grid_out[j][k][0].item(0) < self.mz_bin_low:
                            self.mz_bin_low = self.grid_out[j][k][0].item(0)
                        if self.grid_out[j][k][-1].item(0) > self.mz_bin_high:
                            self.mz_bin_high = self.grid_out[j][k][-1].item(0)

            #create zero array with range of bin indices
            self.mz_bins = np.arange(self.mz_bin_low, self.mz_bin_high, 0.002)
            self.mz_len = len(self.mz_bins)

            #create empty space with dimensions matching grid_out and m/z indices
            self.full_grid_out = np.zeros((np.shape(self.grid_out)[0], np.shape(self.grid_out)[1], self.mz_len))

            #determine and apply mapping of nonzero grid_out values to full_grid_out
            for l in range(len(self.grid_out)):
                for m in range(len(self.grid_out[l])):
                    if len(self.grid_out[l][m]) == 0:
                        pass

                    else:
                        if len(self.grid_out[l][m]) != 0:
                            low_idx = np.searchsorted(self.mz_bins, self.grid_out[l][m][0,0])
                            high_idx = np.searchsorted(self.mz_bins, self.grid_out[l][m][-1,0])
                            if high_idx - low_idx == len(self.grid_out[l][m][:,0]):
                                self.indices = np.arange(low_idx, high_idx, 1)
                            else:
                                self.indices = np.clip(np.searchsorted(self.mz_bins, self.grid_out[l][m][:,0]), 0, len(self.mz_bins)-1)
                        else:
                            self.indices = np.clip(np.searchsorted(self.mz_bins, self.grid_out[l][m][:,0]), 0, len(self.mz_bins)-1)
                        
                        if self.indices[-1] == len(self.mz_bins):
                            self.indices[-1] = self.indices[-1]-1
                        try:
                            self.full_grid_out[l,m][self.indices] = self.grid_out[l][m][:,1]
                        except:
                            ipdb.set_trace()
                            print("this stops the iterator")

            self.full_grid_out = self.full_grid_out

            self.full_gauss_grids = self.gauss(self.full_grid_out)
           
        #Handle concatenated tensor case, kwarg values cannot be computed internally and must be passed 
        else: 
            if not all('dts' in kwargs and \
                'rts' in kwargs and \
                'lows' in kwargs and \
                'highs' in kwargs and \
                'concatenated_grid' in kwargs and \
                'abs_mz_low' in kwargs and \
                'concat_dt_idxs' in kwargs):

                print("Concatenated Tensor Missing Required Values")
                sys.exit()
        
    #Takes tensor input and gaussian filter parameters, outputs filtered data
    def gauss(self, grid, rt_sig = 3, dt_sig = 1):

        gauss_grid = np.zeros(np.shape(grid))
        for i in range(np.shape(grid)[2]):
            gauss_grid[:,:,i] = gaussian_filter(grid[:,:,i],(rt_sig, dt_sig))
        return gauss_grid
    
    #Takes length of mz_bins to interpolated to, and optional gaussian filter parameters
    #Returns the interpolated tensor, length of interpolated axis, interpolated low_lims and high_lims
    def interpolate(self, grid_in, new_mz_len, gauss_params = None):
        
        if gauss_params != None:
            grid = self.gauss(grid_in, gauss_params[0], gauss_params[1])
        else:
            grid = grid_in
            
        test_points = []
        z_axis = np.clip(np.linspace(0, np.shape(grid)[2], new_mz_len), 0, np.shape(grid)[2]-1)
        for n in range(np.shape(grid)[0]):
            for o in range(np.shape(grid)[1]):
                for p in z_axis:
                    test_points.append((n, o, p))
        
        x, y, z = np.arange(np.shape(grid)[0]), np.arange(np.shape(grid)[1]), np.arange(np.shape(grid)[2])
        interpolation_function = sp.interpolate.RegularGridInterpolator(points = [x, y, z], values = grid)
        interpolated_out = interpolation_function(test_points)
        interpolated_out = np.reshape(interpolated_out, (np.shape(grid)[0], np.shape(grid)[1], new_mz_len))
        
        interpolated_bin_mzs = np.linspace(self.mz_bin_low, self.mz_bin_high, new_mz_len)
        interpolated_low_lims = np.searchsorted(interpolated_bin_mzs, self.mz_bins[self.lows])
        interpolated_high_lims = np.searchsorted(interpolated_bin_mzs, self.mz_bins[self.highs])
        
        return [interpolated_out, interpolated_low_lims, interpolated_high_lims]

    #Decomp series takes low n_factors to high n_factors and stores lists of factors in a list within the DataTensor
    def decomposition_series(
        self, 
        n_factors_low, 
        n_factors_high, 
        new_mz_len = None,
        gauss_params = None
        ):

        self.decomps = []
        self.decomp_times = []
        for i in range(n_factors_low, n_factors_high+1):
            t1 = time.time()
            self.decomps.append(self.factorize(i, new_mz_len, gauss_params))
            t2 = time.time()
            self.decomp_times.append(t2-t1)
            print(str(i-n_factors_low+1)+" of "+str(n_factors_high-n_factors_low+1)+" T+"+str(t2-t1))
        
    #Takes single specified data type from a DataTensor object and desired number of factors, returns list of factorization components. Input lows and highs if 
    #using interpolated data, use n_concatenated to pass the number of data tensors concatenated together to the component object
    #Gauss_params must be tuple of len=2
    def factorize(
        self, 
        n_factors, 
        new_mz_len = None, 
        gauss_params = None
        ):
        
        factors = []
        if self.n_concatenated != 1: 
            grid, lows, highs, concat_dt_idxs = self.concatenated_grid, self.lows, self.highs, self.concat_dt_idxs
        else:    
            if new_mz_len != None:
                if gauss_params != None:
                    grid, lows, highs, concat_dt_idxs = interpolate(self.full_grid_out, new_mz_len, gauss_params[0], gauss_params[1]), None
                else: 
                    grid, lows, highs, concat_dt_idxs = interpolate(self.full_grid_out, new_mz_len), None
            else: 
                lows, highs, concat_dt_idxs = self.lows, self.highs, None
                if gauss_params != None:
                    grid = self.gauss(self.full_grid_out, gauss_params[0], gauss_params[1])
                else:
                    grid = self.full_grid_out
                
        
        #Multiplies all values outside of integration box boundaries by 0, good for reducing noise because we only care about data within the bounds
        zero_mult = np.zeros((np.shape(grid)))
        for lo, hi in zip(lows, highs):
            zero_mult[:,:,lo:hi] = 1
        grid *= zero_mult


        nnp = non_negative_parafac(grid, n_factors, init = 'svd')# or 'random' #, n_iter_max = 50) can be used to limit time of calculations at cost of higher reconstruction error
        for i in range(n_factors):
            factors.append(
                Factor(
                    source_file = self.source_file,
                    tensor_idx = self.tensor_idx, 
                    timepoint_idx = self.timepoint_idx, 
                    name = self.name, 
                    charge_states = self.charge_states, 
                    rts = nnp[1][0].T[i], 
                    dts = nnp[1][1].T[i], 
                    mz_data = nnp[1][2].T[i], 
                    factor_idx = i, 
                    n_factors = n_factors, 
                    lows = lows, 
                    highs = highs, 
                    abs_mz_low = self.mz_bin_low, 
                    n_concatenated = self.n_concatenated, 
                    concat_dt_idxs = concat_dt_idxs,
                    total_mass_window = self.total_mass_window
                    )
                )

        return factors
    
###    
### Class - Factor:        
### Holds data from a single component of a non_negative_parafac factorization
### constructed as: (nnp[0].T[i], nnp[1].T[i], nnp[2].T[i], i, n, self.lows, self.highs, self.n_concatenated)
###

class Factor:
    
    def __init__(
            self, 
            source_file,
            tensor_idx, 
            timepoint_idx, 
            name, 
            charge_states, 
            rts, 
            dts, 
            mz_data, 
            factor_idx, 
            n_factors, 
            lows, 
            highs, 
            abs_mz_low, 
            n_concatenated, 
            concat_dt_idxs,
            total_mass_window
            ):

        ###Set Attributes###

        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.name = name
        self.charge_states = charge_states   
        self.rts = rts
        self.dts = dts
        self.mz_data = mz_data
        self.auc = sum(mz_data)
        self.factor_idx = factor_idx
        self.n_factors = n_factors
        self.lows = lows
        self.highs = highs
        self.abs_mz_low = abs_mz_low
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.total_mass_window = total_mass_window

        ###Compute Instance Values###
    
        #integrate within expected peak bounds and create boolean mask of expected peak bounds called grate
        self.integrated_mz_data = []
        self.grate = np.resize(False, (len(self.mz_data)))
        for i, j in zip(self.lows, self.highs):
            self.integrated_mz_data.append(sum(self.mz_data[i:j]))
            self.grate[i:j] = True       
        self.grate_sum = sum(self.mz_data[self.grate])
        
        self.max_rtdt = max(self.rts) * max(self.dts)
        self.outer_rtdt = sum(sum(np.outer(self.rts, self.dts)))
        
#This can be a shared function
        self.integration_box_centers = []
        for i, j in zip(self.lows, self.highs):
            self.integration_box_centers.append(i+((j-i)/2))
           
        self.box_intensities = self.mz_data[self.grate]
        self.max_peak_height = max(self.box_intensities)
        self.mz_peaks = sp.signal.find_peaks(self.mz_data, height = 0.01)[0] # TODO consider replacing with prominence, check in notebook
        #Unused at factor level
        #self.peak_error, self.peaks_chosen = self.peak_error(self.mz_data, self.mz_peaks, self.integration_box_centers, self.max_peak_height)
                
        #this is a poor implementation, at least use list comprehensions   TODO  
        self.box_dist_avg = 0
        for i in range(1,len(self.integration_box_centers)):
            self.box_dist_avg += self.integration_box_centers[i]-self.integration_box_centers[i-1]
        self.box_dist_avg = self.box_dist_avg/(len(self.integration_box_centers)-1)
                
        #Writes to self.isotope_clusters 
        self.find_isotope_clusters(5, height = 0.5) #heuristic height value, should be high-level param TODO - Will require passage through DataTensor class

    #Uses find_window function to identify portions of the integrated mz dimension that look 'isotopic-cluster-like', saves as Factor attribute
    def find_isotope_clusters(self, peak_width, **kwargs):
        self.isotope_clusters = []
        peaks = sp.signal.find_peaks(self.integrated_mz_data, **kwargs)[0]
        if len(peaks) == 0:
            return
        else:
            cluster_idx = 0
            for i in range(len(peaks)):
                integrated_indices = self.find_window(self.integrated_mz_data, peaks[i], peak_width)
                if integrated_indices != None:
                    self.isotope_clusters.append(
                        IsotopeCluster(
                            charge_states = self.charge_states, 
                            factor_mz_data = copy.deepcopy(self.mz_data), 
                            source_file = self.source_file,
                            tensor_idx = self.tensor_idx, 
                            timepoint_idx = self.timepoint_idx, 
                            n_factors = self.n_factors, 
                            factor_idx = self.factor_idx, 
                            cluster_idx = cluster_idx, 
                            low_idx = self.lows[integrated_indices[0]]-math.ceil(self.box_dist_avg/2), 
                            high_idx = self.highs[integrated_indices[1]]+math.ceil(self.box_dist_avg/2), 
                            lows = self.lows, 
                            highs = self.highs, 
                            grate = self.grate, 
                            rts = self.rts, 
                            dts = self.dts, 
                            max_rtdt = self.max_rtdt, 
                            outer_rtdt = self.outer_rtdt, 
                            box_dist_avg = self.box_dist_avg, 
                            abs_mz_low = self.abs_mz_low, 
                            n_concatenated = self.n_concatenated, 
                            concat_dt_idxs = self.concat_dt_idxs,
                            total_mass_window = self.total_mass_window
                            )
                        )
                    cluster_idx += 1
            return

    #heuristically identifies 'things that look like acceptable isotopic clusters' in integrated mz dimension, roughly gaussian allowing some inflection points from noise
    def find_window(self, array, peak_idx, width):
        rflag = True
        lflag = True
        if peak_idx == 0:
            win_low = 0
            lflag = False
        if peak_idx == len(array)-1:
            win_high = len(array)-1
            rflag = False
            
        idx = peak_idx+1 
        if idx < len(array)-1: #Check if idx is last idx
            if array[idx] < array[peak_idx]/5: #Peak is likely not an IC if peak > 5 x neighbors
                win_high = idx
                rflag = False

        while rflag:
            #make sure looking ahead won't throw error
            if idx+1 < len(array):
                #if idx+1 goes down, and its height is greater than 20% of the max peak
                if array[idx+1] < array[idx] and array[idx+1] > array[peak_idx]/5:
                    idx += 1
                #if above check fails, test conditions separately
                else: 
                    if array[idx+1] < array[peak_idx]/5:
                        win_high = idx
                        rflag = False
                    else:
                        #first check if upward point is more than 5x the height of the base peak
                        if array[idx+1] < array[peak_idx]*5:
                            #look one point past upward point, if its below the last point and the next point continues down, continue
                            if idx+2 < len(array):
                                if array[idx+2] < array[idx+1]:
                                    if idx+3 < len(array):
                                        if array[idx+3] < array[idx+2]:
                                            idx += 3
                                        else: #point 3 ahead goes up, do not keep sawtooth tail
                                            win_high = idx
                                            rflag = False   
                                    else: #point 2 past idx is end of array, set as high limit
                                        win_high = idx+2
                                        rflag = False
                                else: #points continue increasing two ahead of idx, stop at idx
                                    win_high = idx
                                    rflag = False
                            else: #upward point is end of array, do not keep
                                win_high = idx
                                rflag = False
                        else: #upward point is major spike / base peak is minor, end search
                            win_high = idx
                            rflag = False
            else: #idx is downward and end of array
                    win_high = idx
                    rflag = False        

        idx = peak_idx-1
        if idx >= 0:
            if array[idx] < array[peak_idx]/5:
                win_low=idx
                lflag=False

        while lflag:
            if idx-1 >= 0: #make sure looking ahead won't throw error
                if array[idx-1] < array[idx] and array[idx-1] > array[peak_idx]/5: #if idx-1 goes down, and its height is greater than 20% of the max peak
                    idx -= 1
                #if above check fails, test conditions separately
                else:
                    if array[idx-1] < array[peak_idx]/5:
                        win_low = idx
                        lflag = False
                    else:
                        #first check if upward point is more than 5x the height of the base peak
                        if array[idx-1] < array[peak_idx]*5:
                            #look one point past upward point, if its below the last point and the next point continues down, continue
                            if idx-2 >= 0:
                                if array[idx-2] < array[idx]:
                                    if idx-3 >= 0:
                                        if array[idx-3] < array[idx-2]:
                                            idx -= 3
                                        else:  #point 3 ahead goes up, do not keep sawtooth tail
                                            win_low = idx
                                            lflag = False
                                    else: #point 2 past idx is end of array, set as high limit
                                        win_low = idx-2
                                        lflag = False
                                else: #points continue increasing two ahead of idx, stop at idx
                                    win_low = idx
                                    lflag = False
                            else: #upward point is start of array, do not keep
                                win_low = idx
                                lflag = False
                        else: #upward point is major spike / base peak is minor, end search
                            win_low = idx
                            lflag = False
            else: #idx is start of array
                win_low = idx
                lflag = False    
                    
        if win_high-win_low < width:
            return None
        else:
            return [win_low, win_high]
    
###
### Class - IsotopeCluster:
### Contains mz data of mz region identified to have isotope-cluster-like characteristics, stores full data of parent factor
###

class IsotopeCluster:
    
    def __init__(
        self, 
        charge_states, 
        factor_mz_data,
        source_file, 
        tensor_idx, 
        timepoint_idx, 
        n_factors, 
        factor_idx, 
        cluster_idx, 
        low_idx, 
        high_idx, 
        total_mass_window,
        lows, 
        highs, 
        grate, 
        rts, 
        dts, 
        max_rtdt, 
        outer_rtdt, 
        box_dist_avg, 
        abs_mz_low, 
        n_concatenated, 
        concat_dt_idxs,
        ):
        

        ###Set Attributes###

        self.charge_states = charge_states
        self.factor_mz_data = factor_mz_data
        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.n_factors = n_factors
        self.factor_idx = factor_idx
        self.cluster_idx = cluster_idx
        self.low_idx = low_idx
        self.high_idx = high_idx
        self.lows = lows
        self.highs = highs
        self.grate = grate
        self.rts = rts
        self.dts = dts
        self.max_rtdt = max_rtdt
        self.outer_rtdt = outer_rtdt
        self.box_dist_avg = box_dist_avg
        self.abs_mz_low = abs_mz_low
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.total_mass_window = total_mass_window


        ###Calculate Scoring Requirements###

        #Create array of expected peak positions
        self.integration_box_centers = []
        for i, j in zip(self.lows, self.highs):
            self.integration_box_centers.append(i+((j-i)/2))
        
        #prune factor_mz to get window around cluster that is consistent between charge-states
        self.cluster_mz_data = self.factor_mz_data
        self.cluster_mz_data[0:self.low_idx] = 0
        self.cluster_mz_data[self.high_idx:] = 0
        
        #integrate area of IC
        self.auc = sum(self.cluster_mz_data)
        
        #Values of isotope cluster that fall within the grate of expected peak bounds
        self.box_intensities = self.cluster_mz_data[self.grate]
        self.grate_sum = sum(self.box_intensities)
        
        #identify peaks and find error from expected peak positions using raw mz
        self.mz_peaks = sp.signal.find_peaks(self.factor_mz_data, distance = self.box_dist_avg/10)[0]
        self.max_peak_height = max(self.box_intensities)
        self.peak_error, self.peaks_chosen = self.find_peak_error(self.cluster_mz_data, 
                                                        self.mz_peaks, 
                                                        self.integration_box_centers[np.searchsorted(self.lows, self.low_idx): np.searchsorted(self.highs, self.high_idx)], 
                                                        self.max_peak_height)
        
        #subtract baseline from IC mz values, recompute intrinsic values with new array
        self.baseline = peakutils.baseline(self.cluster_mz_data[self.low_idx:self.high_idx], 6) #6 degree curve seems to work well
        self.baseline_subtracted_mz = self.cluster_mz_data
        self.baseline_subtracted_mz[self.low_idx: self.high_idx] = self.cluster_mz_data[self.low_idx: self.high_idx]-self.baseline
        self.baseline_auc = sum(self.baseline_subtracted_mz)
        self.log_baseline_auc = np.log(self.baseline_auc)
        self.baseline_box_intensities = self.baseline_subtracted_mz[self.grate]
        self.baseline_grate_sum = sum(self.baseline_box_intensities)
        self.baseline_max_peak_height = max(self.baseline_box_intensities)
        self.baseline_peak_error, self.baseline_peaks_chosen = self.find_peak_error(
                                                                          self.baseline_subtracted_mz,
                                                                          self.mz_peaks,
                                                                          self.integration_box_centers[np.searchsorted(self.lows, self.low_idx): np.searchsorted(self.highs, self.high_idx)],
                                                                          self.baseline_max_peak_height
                                                                          )

        #create integrated mz array, indexed by integration box
        baseline_int_mz = []
        for lo, hi in zip(self.lows, self.highs):
            baseline_int_mz.append(sum(self.baseline_subtracted_mz[lo:hi]))
        self.baseline_integrated_mz = np.asarray(baseline_int_mz)

        #Cache int_mz and rt scoring values
        self.baseline_integrated_mz_norm = self.baseline_integrated_mz/np.linalg.norm(self.baseline_integrated_mz)
        self.baseline_integrated_mz_com = sp.ndimage.measurements.center_of_mass(self.baseline_integrated_mz)[0] #COM in IC integrated bin dimension
        self.baseline_integrated_mz_std = np.average((np.arange(len(self.baseline_integrated_mz)) - self.baseline_integrated_mz_com)**2, weights = self.baseline_integrated_mz)**0.5

        self.rt_norm = self.rts/np.linalg.norm(self.rts)
        self.rt_com = sp.ndimage.measurements.center_of_mass(self.rts)[0]
        
        #Cache DT values
        #If DT is concatenated, return list of coms and norms of single rts relative to bin numbers, a single_dt distribution starts at 0. If only one charge state, return list of len=1
        if self.concat_dt_idxs is not None:
            single_dts = []
            #generate list of single dts
            single_dts.append(self.dts[:self.concat_dt_idxs[0]])
            for i in range(len(self.charge_states)-1):
                single_dts.append(self.dts[self.concat_dt_idxs[i]:self.concat_dt_idxs[i+1]])
            
            self.single_dts = single_dts
            self.dt_coms = [sp.ndimage.measurements.center_of_mass(dt)[0] for dt in single_dts]
            self.dt_norms = [dt/np.linalg.norm(dt) for dt in single_dts]
        else:
            self.dt_coms = [sp.ndimage.measurements.center_of_mass(self.dts)[0]]
            self.dt_norms = [self.dts/np.linalg.norm(self.dts)] 
        

        if self.n_concatenated == 1:
            self.abs_mz_com = self.find_mz_com(self.total_mass_window) 
        else: 
            self.abs_mz_com = "Concatenated, N/A, see IC.baseline_integrated_mz_com"


        #format useful values to be read by pandas
        self.info_tuple = (
                        self.source_file, #Filename of data used to create parent DataTensor
                        self.tensor_idx, #Library master list row of parent-DataTensor
                        self.n_factors, #Number of factors in parent decomposition
                        self.factor_idx, #Index of IC parent-factor in DataTensor.decomps[]
                        self.cluster_idx, #Index of IC in parent-factor.isotope_clusters[]
                        self.charge_states, #List of charge states in IC
                        self.n_concatenated, #number of source tensors IC parent-DataTensor was made from
                        self.low_idx, #Low bin index corresponding to Factor-level bins
                        self.high_idx, #High bin index corresponding to Factor-level bins
                        self.baseline_auc, #Baseline-subtracted AUC (BAUC)
                        self.baseline_grate_sum, #Baseline-subtracted grate area sum (BGS)
                        self.baseline_peak_error, #Baseline-subtracted version of peak-error (BPE)
                        self.baseline_integrated_mz_com, #Center of mass in added mass units
                        self.abs_mz_com, #IsotopicCluster center-of-mass in absolute m/z dimension
                        self.rts, #Array of RT values
                        self.dts, #Array of DT values, if a tensor is concatenated,this is taken from the last tensor in the list, can be seen in tensor_idx
                        np.arange(0, len(self.baseline_integrated_mz), 1),
                        self.baseline_integrated_mz #Array of baseline-subtracted integrated mz intensity values
                        )

        #instantiate to make available for setting in PathOptimizer
        self.bokeh_tuple = None #bokeh plot info tuple
        self.single_sub_scores = None #differences between winning score and score if this IC were substituted, list of signed values
        self.undeut_ground_dot_products = None

    #uses internal 
    def find_mz_com(self, tensor_mass_range):
        factor_mz_range = tensor_mass_range/self.charge_states[0]
        factor_mz_bin_step = factor_mz_range/len(self.factor_mz_data)
        left_mz_dist = (self.highs[math.floor(self.baseline_integrated_mz_com)] - self.integration_box_centers[math.floor(self.baseline_integrated_mz_com)]) * factor_mz_bin_step #MZ dist from center to right bound, times bin_to_mz factor
        right_mz_dist = (self.integration_box_centers[math.floor(self.baseline_integrated_mz_com)+1] - self.lows[math.floor(self.baseline_integrated_mz_com)+1]) * factor_mz_bin_step #MZ dist from left bound to center, times bin_to_mz factor

        if (self.baseline_integrated_mz_com - math.floor(self.baseline_integrated_mz_com)) <= 0.5:
            major = self.abs_mz_low + (factor_mz_bin_step * self.integration_box_centers[math.floor(self.baseline_integrated_mz_com)])
            minor = left_mz_dist * (self.baseline_integrated_mz_com - math.floor(self.baseline_integrated_mz_com))
            abs_mz_com = major + minor
        else:
            major = self.abs_mz_low + (factor_mz_bin_step * self.lows[math.floor(self.baseline_integrated_mz_com)] + 1)
            minor = right_mz_dist * (self.baseline_integrated_mz_com - math.floor(self.baseline_integrated_mz_com) - 0.5)
            abs_mz_com = major + minor
        
        return abs_mz_com
    
    #calculates sum of distances between nearest prominent peaks and expected peak centers in IC
    def find_peak_error(
        self,
        source, 
        mz_peaks, 
        integration_box_centers,
        max_peak_height
        ):

        peak_error = 0
        peaks_chosen = []
        peaks_total_height = 0
        match_idx = np.searchsorted(mz_peaks, integration_box_centers)
        if len(mz_peaks) > 0:
            for i in range(len(match_idx)):
                #handle peaks list of length 1
                if len(mz_peaks) == 1:
                    peak_error += abs(integration_box_centers[i]-mz_peaks[0])*(source[mz_peaks[0]])
                    peaks_chosen.append(mz_peaks[0])
                    peaks_total_height += source[mz_peaks[0]]
                else:                  
                    #check if place to be inserted is leftmost of peaks
                    if match_idx[i] == 0:
                        peak_error += abs(integration_box_centers[i]-mz_peaks[match_idx[i]])*(source[mz_peaks[match_idx[i]]])
                        peaks_chosen.append(mz_peaks[match_idx[i]])
                        peaks_total_height += source[mz_peaks[match_idx[i]]]
                    else:
                        #check if insertion position is rightmost of peaks
                        if match_idx[i] == len(mz_peaks):
                            peak_error += abs(integration_box_centers[i]-mz_peaks[-1])*(source[mz_peaks[-1]])
                            peaks_chosen.append(mz_peaks[-1])
                            peaks_total_height += source[mz_peaks[-1]]
                        else:
                            #handle case where distances between peaks are the same, pick biggest peak
                            if abs(integration_box_centers[i]-mz_peaks[match_idx[i]]) == abs(integration_box_centers[i]-mz_peaks[match_idx[i]-1]):
                                peak_error += max([abs(integration_box_centers[i]-mz_peaks[match_idx[i]])*(source[mz_peaks[match_idx[i]]]),
                                            abs(integration_box_centers[i]-mz_peaks[match_idx[i]-1])*(source[mz_peaks[match_idx[i]-1]])])
                                if (abs(integration_box_centers[i]-mz_peaks[match_idx[i]])*(source[mz_peaks[match_idx[i]]]) >
                                    abs(integration_box_centers[i]-mz_peaks[match_idx[i]-1])*(source[mz_peaks[match_idx[i]-1]])):
                                    peaks_chosen.append(mz_peaks[match_idx[i]])
                                    peaks_total_height += source[mz_peaks[match_idx[i]]]
                                else:
                                    peaks_chosen.append(mz_peaks[match_idx[i]-1])
                                    peaks_total_height += source[mz_peaks[match_idx[i]-1]]
                            else:
                                #only need to check left hand side differences because of left-hand default of searchsorted algorithm
                                #now check which peak is closer, left or right. This poses problems as there may be very close peaks which are not
                                #actually significant in height but which pass filtering.
                                if abs(integration_box_centers[i]-mz_peaks[match_idx[i]]) < abs(integration_box_centers[i]-mz_peaks[match_idx[i]-1]):
                                    peak_error += abs(integration_box_centers[i]-mz_peaks[match_idx[i]])*(source[mz_peaks[match_idx[i]]])
                                    peaks_chosen.append(mz_peaks[match_idx[i]])
                                    peaks_total_height += source[mz_peaks[match_idx[i]]]
                                else:
                                    peak_error += abs(integration_box_centers[i]-mz_peaks[match_idx[i]-1])*(source[mz_peaks[match_idx[i]-1]])
                                    peaks_chosen.append(mz_peaks[match_idx[i]-1])
                                    peaks_total_height += source[mz_peaks[match_idx[i]-1]]
            
            box_dist_total = 0
            for i in range(1,len(integration_box_centers)):
                box_dist_total += integration_box_centers[i]-integration_box_centers[i-1]
            
            
            peak_error = (peak_error / peaks_total_height / (box_dist_total/(len(integration_box_centers)-1)))
            return peak_error, peaks_chosen
    
###
### Class - TensorGenerator:
### Collects and clusters all hdx timepoint files associated with a single protein name, creates and factorizes DataTensors, extracts isotopic clusters for PathOptimizer
###

class TensorGenerator:  
    """
    Accepts a protein name associated with .mzML files in a specified directory, collects and clusters those files before generating DataTensors and IsotopeClusters from them.
    
    name = <string> protein file name 
    path = <string> directory to look for .mzML files in
    library_info = <pd.DataFrame> contains all information about identified proteins 
    inputs = <list of strings> filenames of all pickle.zlib tensors for protein name
    n_factors = <int> of 
    """
    
    ###Class Attributes###
    hd_mass_diff = 1.006277
    c13_mass_diff = 1.00335

    def __init__(
        self, 
        name, 
        library_info,
        inputs, #passes a flat list of limit-readable tensors for all charge-states of all timepoints of a single protein RT-group
        timepoints,
        **kwargs
        ):

        ###Set Instance Attributes###

        self.name = name
        self.library_info = library_info
        self.inputs = inputs
        self.timepoints = timepoints
        

        if kwargs is not None: #TODO: Default control values should be placed in a check statement below (if hasattr(self, 'name of value'): Do) do this for other kwargs fxns
            for key in kwargs.keys():
                setattr(self, key, kwargs[key])

        if not hasattr(self, 'low_mass_margin'):
            self.low_mass_margin = 10
        if not hasattr(self,'high_mass_margin'):
            self.high_mass_margin = 17
        if not hasattr(self, 'ppm_radius'):
            self.ppm_radius = 30
        if not hasattr(self, 'n_factors_low'):
            self.n_factors_low = 1
        if not hasattr(self, 'n_factors_high'):
            self.n_factors_high = 3
        if not hasattr(self, 'gauss_params'):
            self.gauss_params = (3,1)
            
            
        self.max_peak_center = len(self.library_info.loc[self.library_info['name'] == self.name]['sequence'].values[0])
        self.total_isotopes = self.max_peak_center+self.high_mass_margin 
        self.total_mass_window = self.low_mass_margin+self.total_isotopes

        self.est_peak_gaps = [0] + list(np.linspace(self.c13_mass_diff, self.hd_mass_diff, 7)) + [self.hd_mass_diff for x in range(self.total_isotopes - 8)]
        self.cum_peak_gaps = np.cumsum(self.est_peak_gaps)

        self.mz_lows = []
        self.mz_highs = []
        self.mz_centers = []
        self.low_lims = []
        self.high_lims = []

        ###This block is redundant, limit to protein name indices? - make lists into dicts with keys corresponding to charge state indices TODO
        count = 0
        for i in range(len(self.library_info)):
            self.mz_centers.append(self.library_info['obs_mz'].values[i] + (self.cum_peak_gaps / self.library_info['charge'].values[i]))
            self.mz_lows.append(self.library_info['obs_mz'].values[i] - (self.low_mass_margin/self.library_info['charge'].values[i]))
            self.mz_highs.append(self.library_info['obs_mz'].values[i] + (self.total_isotopes/self.library_info['charge'].values[i]))

            self.low_lims.append(self.mz_centers[count] * ((1000000.0 - self.ppm_radius)/1000000.0))
            self.high_lims.append(self.mz_centers[count] * ((1000000.0 + self.ppm_radius)/1000000.0))
            count += 1
        ###End block

        self.mz_low = min(self.mz_lows)
        self.mz_high = max(self.mz_highs)
    
    #Instantiates all DataTensors identified as part of a protein RT-group (Collections of peaks identified as the POI within a ), factorizes Tensors,
    #identifies IsotopicClusters within Factors, scores ICs, concatenates Tensors if deemed beneficial/possible, saves all timepoint ICs and Factors internally. 
    #Must be called externally on an instance of TensorGenerator.   
    def generate_tensors(self):

        #will both have nested-list structure with outer-length = n_timepoints
        all_tp_clusters = []
        all_tp_factors = []

        #Iterate over timepoints in ascending order
        for tp in range(len(self.timepoints)):
            print("Starting "+str(self.timepoints[tp])+" seconds, "+str(tp+1)+" of "+str(len(self.timepoints)))
            tp_factors = []
            #iterate over each replicate of a given timepoint
            for fn in self.inputs[tp]:
                print("Sourcefile: "+fn+", "+str(self.inputs[tp].index(fn)+1)+" of "+str(len(self.inputs[tp])))

                #grab tensor datafiles from input based on their source MS-datafile
                charges = [tensor_path for tensor_path in self.inputs[tp] if fn in tensor_path]
                
                DataTensors = []
                
                #iterate over each charge-state of the rt-group
                for charge in charges:
                    print("File: "+charge)
                    
                    lib_idx = int(charge.split("/")[-1].split("_")[0]) # expects format: path/to/{LibraryIdx}_{protName}_{tp}.cpickle.zlib
                    output = limit_read(charge)

                    #Instantitate DataTensor
                    newDataTensor = DataTensor(
                                        source_file = charge,
                                        tensor_idx = lib_idx,
                                        timepoint_idx = tp, 
                                        name = self.name, 
                                        total_mass_window = self.total_mass_window,
                                        n_concatenated = 1,
                                        charge_states = [self.library_info["charge"].values[lib_idx]], 
                                        rts = output[0], 
                                        dts= output[1], 
                                        seq_out = output[2], 
                                        int_seq_out = None
                                        )

                    #translate low/high limits on expected peak locations from absolute mz to DataTensor instance's mz bin range TODO - Change to dict to remove redundancy
                    newDataTensor.lows = np.searchsorted(newDataTensor.mz_bins, self.low_lims[lib_idx])
                    newDataTensor.highs = np.searchsorted(newDataTensor.mz_bins, self.high_lims[lib_idx])
                    
                    newDataTensor.decomposition_series(n_factors_low = self.n_factors_low, n_factors_high = self.n_factors_high, gauss_params = self.gauss_params)
                    DataTensors.append(newDataTensor)

                fn_factors = []
                fn_clusters = []
                fn_cluster_info = []
                for tensor in DataTensors:
                    for decomp in tensor.decomps:
                        for factor in decomp:
                            fn_factors.append(factor)
                            for ic in factor.isotope_clusters:
                                fn_clusters.append(ic)
                                fn_cluster_info.append(ic.info_tuple) 

                try:
                    protein_clusters = pd.DataFrame(fn_cluster_info)
                    protein_clusters.columns = [
                        "source_file",
                        "tensor_idx", 
                        "n_factors", 
                        "factor_idx", 
                        "cluster_idx",
                        "charge_states", 
                        "n_concatenated",
                        "mz_bin_low", #from factor bins
                        "mz_bin_high", #from factor bins
                        "baseline_subtracted_area_under_curve", 
                        "baseline_subtracted_grate_sum",  
                        "baseline_subtracted_peak_error",
                        "baseline_integrated_mz_com", 
                        "abs_mz_com",
                        "rts", 
                        "dts", 
                        "int_mz_x",
                        "int_mz_y"]
                except:
                    ipdb.set_trace()

                top_fn_clusters = protein_clusters.sort_values(by = ['baseline_subtracted_peak_error'], ascending = False)[:math.ceil(len(protein_clusters)/5)]

                #cheks that each tensor gives 85% of the even split contribution to the top 20% of clusters, if less the tensor is not included in the concatenated data
                #85% is an arbitrary choice of cutoff value, empirically determine point at which adding DT to combined harms quality
                concat_ins = []
                concat_charges = []
                concat_dt_lens = []
                if len(output)>0:    
                    for i in range(len(DataTensors)):
                        contribution = len(top_fn_clusters.loc[top_fn_clusters['tensor_idx'] == i])
                        if contribution > (math.floor(math.ceil(len(protein_clusters)/5)/len(DataTensors))*0.85):
                            concat_ins.append(DataTensors[i].interpolate(DataTensors[i].full_grid_out, 5000, (3,1))[0])
                            concat_charges.append(DataTensors[i].charge_states[0])
                            concat_dt_lens.append(len(DataTensors[i].dts))
                    if len(concat_ins)>1:
                        print("CONCATENATING")
                        try:
                            interpolated_lows, interpolated_highs = DataTensors[0].interpolate(DataTensors[0].full_grid_out, new_mz_len = 5000, gauss_params = self.gauss_params)[1:3]
                            concatenated = np.concatenate(concat_ins, axis = 1)
                            concat_dt_idxs = np.cumsum(concat_dt_lens)
                            mz_bin_low = min([x.mz_bin_low for x in DataTensors])
                           
                            #RT is taken from last tensor read in, consider averaging?
                            DataTensors.append(
                                DataTensor(
                                    source_file = "Concatenated, see IC.charge_states",
                                    tensor_idx = "Concatenated, see IC.charge_states",
                                    timepoint_idx = tp,
                                    name = self.name, 
                                    total_mass_window = self.total_mass_window,
                                    charge_states = concat_charges, 
                                    rts = output[0], 
                                    dts = concatenated, 
                                    n_concatenated = len(concat_ins), 
                                    concat_dt_idxs = concat_dt_idxs, 
                                    concatenated_grid = concatenated, 
                                    lows = interpolated_lows, 
                                    highs = interpolated_highs, 
                                    abs_mz_low = mz_bin_low
                                    )
                                )
                            DataTensors[-1].decomposition_series(n_factors_low, n_factors_high)   

                            #Add concatenated factors to fn_factors, ICs emptied later
                            for decomp in DataTensors[-1].decomps:
                                for factor in decomp:
                                    fn_factors.append(factor)
                        
                        except ValueError:
                            #handle length mismatches in future, now just skip it and show the error TODO
                            #often a +- 1 error
                            print("Length mistmatch in the tensors to be combined")

                #append list of fn_factor lists to tp_factors
                tp_factors.append(fn_factors)
            #collect all tp_factor lists into one master list
            all_tp_factors.append(tp_factors) # 3D: Timepoints X Replicates X Factors
        
        #empty ICs from factors into convenient IC-only list
        for tp in all_tp_factors:
            tp_ic_buffer = [] #flattens fn lists into single tp list
            for fn in tp:
                for factor in fn:
                    for ic in factor.isotope_clusters:
                        tp_ic_buffer.append(ic)
            all_tp_clusters.append(tp_ic_buffer) #2D: timepoints X ICs

        #Save Factors and ICs as instance attributes
        self.all_tp_factors, self.all_tp_clusters = all_tp_factors, all_tp_clusters

###
### Class - PathOptimizer:
### Optimizes HDX timeseries of Isotopic Clusters by generating a set of starting timeseries based on possible trajectory though the integrated mz dimension,
### Then iteratively uses a set of timeseries scoring functions to make the best single substitution until each starting series is optimized for score.
### Timeseries with best score at the end of all minimizations is selected as the winning path, which is output along with the alternatives for each timepoint.
###

class PathOptimizer:
    """
    Generates sample 'paths' - trajectories through HDX timeseries - optimizes 'path' through hdx timeseries for all sample paths and selects an overall winning path. 
    
    all_tp_clusters = <list> of <lists> of <TA.isotope_cluster>s for each HDX timepoint,
    prefilter = 0 =< <float> < 1, proportion of total isotope_clusters to keep after filtering, if prefilter=0 prefiltering will be skipped
    n_samples = <int> number of semi-randomly initialized sample paths to generate
    
    """

    #TODO: add info_tuple-like struct, eventually change IC, PO, and bokeh related scoring systems to use dicts? Dicts would make changing column names simpler.

    def __init__(
        self, 
        name, 
        all_tp_clusters, 
        library_info, 
        timepoints,
        n_undeut_runs,
        prefilter = 0,
        old_data_dir = None,
        **kwargs
        ):
 
        #Set score weights
        self.int_mz_std_rmse_weight = 1
        self.baseline_peak_error_weight = 100
        self.delta_mz_rate_weight = 1.65 #formerly 2
        self.int_mz_rot_fit_weight = 5
        self.dt_ground_fit_weight = 25
        self.rt_ground_fit_weight = 5
        self.rt_ground_rmse_weight = 10
        self.dt_ground_rmse_weight = 10
        self.auc_ground_rmse_weight = 5

        #Set internal variables
        self.name = name
        self.all_tp_clusters = all_tp_clusters
        self.library_info = library_info
        self.prefilter = prefilter
        self.timepoints = timepoints
        self.n_undeut_runs = n_undeut_runs
        self.max_peak_center = len(self.library_info.loc[self.library_info['name'] == self.name]['sequence'].values[0]) #TODO: Fix. Bad solution, keep an eye on this, could break if only one charge in RT-group? Maybe always list, double check
        
        self.old_data_dir = old_data_dir
        self.old_files = None

        self.rt_com_cv = None
        self.dt_com_cv = None
        self.rt_error_rmse = None
        self.dt_error_rmse = None

        self.gather_old_data()
        self.select_undeuterated()
        self.precalculate_fit_to_ground()

        #TODO: Refine prefilter after exploring ic attribute correlations to downstream outcome metrics, find best predictors for outcomes. 
        self.prefilter_ics()
        self.generate_sample_paths()

    def gather_old_data(self):
        if self.old_data_dir is not None:
            self.old_files = sorted([fn for fn in glob.iglob(self.old_data_dir+"*.pickle") if "_".join(self.name.split("/")[-1].split("_")[:3]) in fn])
            self.old_data = []
            for fn in self.old_files:
                ts = pickle.load(open(fn, 'rb'))
                ts['charge'] = int(fn.split('.')[-3][-1])
                ts['delta_mz_rate'] = self.gabe_delta_mz_rate(ts['centroid'])
                ts['major_species_widths'] = [len(np.nonzero(x)[0]) for x in ts['major_species_integrated_intensities']]
                self.old_data.append(ts)
            
    def select_undeuterated(
        self, 
        all_tp_clusters = None, 
        library_info = None, 
        name = None, 
        n_undeut_runs = None
        ):

        """
        Selects undeuterated isotopic cluster which best matches theoretically calculated isotopic distribution for POI sequence, for each observed charge state of the POI
        all_tp_clusters = TensorG enerator attribute, e.g. T1 = TensorGenerator(...)\n select_undeuterated(T1.all_tp_clusters)
        n_undeut_runs = number of undeuterated HDX runs included in the 'library_info' master csv
        """

        if all_tp_clusters is None:
            all_tp_clusters = self.all_tp_clusters

        if name is None:
            name = self.name

        if library_info is None:
            library_info = self.library_info

        if n_undeut_runs is None:
            n_undeut_runs = self.n_undeut_runs

        HX1 = hxtools.hxprot(seq=library_info.loc[library_info['name']==name]['sequence'].values[0])

        if self.old_data_dir is not None: #if comparing to old data, save old-data's fits in-place CONSIDER OUTPUTTING TO SNAKEMAKE DIR
            #open first three (undeut) dicts in list, store fit to theoretical dist
            for charge_dict in self.old_data:
                undeut_amds = [{'major_species_integrated_intensities': charge_dict['major_species_integrated_intensities'][i]} for i in range(3)]#hardcode for gabe's undeut idxs in list   
                charge_dict['fit_to_theo_dist'] = max([HX1.calculate_isotope_dist(d) for d in undeut_amds])

        undeuts = []
        for undeut in all_tp_clusters[:n_undeut_runs]:
            for ic in undeut:
                undeuts.append(ic)
        dot_products = []

        #make list of all normed dot products between an undeuterated IC and the theoretical distribution, calculate_isotope_dist requires Pandas input
        for ic in undeuts:
            df = pd.DataFrame(ic.baseline_integrated_mz, columns = ['major_species_integrated_intensities'])
            fit = HX1.calculate_isotope_dist(df)
            ic.undeut_ground_dot_product = fit
            dot_products.append((fit, ic.charge_states))

        #Append final (0, 0) to be called by charge_idxs which are not in the charge group for a single loop iteration
        dot_products.append((0,0))
        charges = list(set(np.concatenate([ic.charge_states for ic in undeuts])))
        out = dict.fromkeys(charges)
        charge_fits = dict.fromkeys(charges)
        for charge in charges:
            #print(charge)
            #Create sublist of undeuts with single charge state and same shape as undeuts, use -1 for non-matches to retain shape of list
            #[charge] == dot_products[i][1] ensures we only pick undeut_grounds from unconcatenated DataTensors, this saves trouble in undeut comparisons
            charge_idxs = []
            for i in range(len(dot_products)):
                if [charge] == dot_products[i][1]:
                    charge_idxs.append(i)
                else:
                    charge_idxs.append(-1)
            #print(charge_idxs)
            #print(np.asarray(dot_products)[charge_idxs])
            
            #Select best fit of charge state, append to output
            best = undeuts[charge_idxs[np.argmax(np.asarray(dot_products)[charge_idxs][:, 0])]]
            
            out[charge] = best
            charge_fits[charge] = max(np.asarray(dot_products)[charge_idxs][:, 0])
        
        self.undeut_grounds = out
        self.undeut_ground_dot_products = charge_fits
    
    #THIS CAN PROBABLY BE REMOVED? MAYBE PREFILTER BY BPE? TODO
    def prefilter_ics(
        self, 
        all_tp_clusters = None, 
        prefilter = None
        ):

        if all_tp_clusters is None:
            all_tp_clusters = self.all_tp_clusters
        if prefilter is None:
            prefilter = self.prefilter

        if self.prefilter == 0:
            outer = []
            uns = []
            for tp in self.all_tp_clusters[:self.n_undeut_runs]:
                for ic in tp:
                    uns.append(ic)
            outer.append(uns)
            for tp in self.all_tp_clusters[self.n_undeut_runs:]:
                outer.append(tp)
            self.prefiltered_ics = outer
            
        else:
            #prefiltered_ics has same shape as all_tp_clusters, but collapses 3 UN timepoints into 1
            prefiltered_ics = []

            #pre_scoring has same order as a single tp of all_tp_clusters
            pre_scoring = []

            #idx_scores contains the sum of the weighted BPE/BGAR sorted indices of ics in prefiltered_ics, with pre_scoring order
            idx_scores = []

            #Prefilitering UN clusters by integration box alignment to low mz
            prefiltered_ics.append([])
            for tp in all_tp_clusters[:3]:
                for ic in tp:
                    if np.argmax(ic.baseline_integrated_mz) <= 5:
                        prefiltered_ics[0].append(ic)

            #Further filter each tp by BGAR and PeakError

            #Decorate
            for ic in prefiltered_ics[0]:
                pre_scoring.append((ic.baseline_auc, ic.baseline_peak_error))

            #Sort, add weights here for tuning
            AUC_sorted = sorted(pre_scoring, key = lambda x: x[0], reverse = True)
            BPE_sorted = sorted(pre_scoring, key = lambda x: x[1])

            #Undecorate
            for tup in pre_scoring:
                idx_scores.append(sum([AUC_sorted.index(tup), BPE_sorted.index(tup)]))
            idx_scores_sorted = sorted(idx_scores, key = lambda x: x)
            score_filter = idx_scores_sorted[math.ceil(len(idx_scores_sorted)*prefilter)]

            undec_buffer = []
            for i in range(len(idx_scores)):
                if idx_scores[i] <= score_filter:
                    undec_buffer.append(prefiltered_ics[0][i])

            #Reset prefiltered UN tp to BPE and BGAR filtered ics
            prefiltered_ics[0] = undec_buffer

            #Repeat process for remaining timepoints
            for tp in range(self.n_undeut_runs, len(all_tp_clusters)):
                #Reset tp vars
                pre_scoring = []
                idx_scores = []

                #Decorate
                for ic in all_tp_clusters[tp]:
                    pre_scoring.append((ic.baseline_auc, ic.baseline_peak_error))

                #Sort, add weights for tuning
                AUC_sorted = sorted(pre_scoring, key = lambda x: x[0], reverse = True)
                BPE_sorted = sorted(pre_scoring, key = lambda x: x[1])

                #Undecorate
                for tup in pre_scoring:
                    idx_scores.append(sum([AUC_sorted.index(tup), BPE_sorted.index(tup)]))
                idx_scores_sorted = sorted(idx_scores, key = lambda x: x)
                score_filter = idx_scores_sorted[math.ceil(len(idx_scores_sorted)*prefilter)]

                tp_buffer = []
                for i in range(len(idx_scores)):
                    if idx_scores[i] <= score_filter:
                        tp_buffer.append(all_tp_clusters[tp][i])

                prefiltered_ics.append(tp_buffer)

            self.prefiltered_ics = prefiltered_ics   


    def precalculate_fit_to_ground(
        self, 
        all_tp_clusters = None, 
        undeut_grounds = None
        ):

        if all_tp_clusters is None:
            all_tp_clusters = self.all_tp_clusters
        if undeut_grounds is None:
            undeut_grounds = self.undeut_grounds

        for timepoint in all_tp_clusters:
            for ic in timepoint:
                dt_ground_errs = []
                rt_ground_errs = []
                dt_ground_fits = []
                rt_ground_fits = []

                #TODO: Rework this to make better sense

                for i in range(ic.n_concatenated):
                    
                    undeut = undeut_grounds[ic.charge_states[i]]
                    
                    rt_ground_errs.append(abs(ic.rt_com - undeut.rt_com))
                    
                    #Trys all dt_coms (for concatenated tensors) keeps best error
                    dt_ground_errs.append(min([abs(ic.dt_coms[i] - undeut.dt_coms[j]) for j in range(len(undeut.dt_coms))])) 


                    diff = len(ic.dt_norms[i]) - len(undeut_grounds[ic.charge_states[i]].dt_norms[0])
                    if diff == 0:
                        dt_ground_fits.append(np.dot(ic.dt_norms[i], undeut_grounds[ic.charge_states[i]].dt_norms[0]))
                    else:
                        if diff > 0:
                            u_buff = copy.copy(undeut_grounds[ic.charge_states[i]].dt_norms[0])
                            u_buff = np.append(u_buff, [0]*diff)
                            dt_ground_fits.append(np.dot(ic.dt_norms[i], u_buff))
                        else:
                            dt_ground_fits.append(np.dot(ic.dt_norms[i], undeut_grounds[ic.charge_states[i]].dt_norms[0][:diff]))

                    diff = len(ic.rt_norm) - len(undeut_grounds[ic.charge_states[0]].rt_norm)
                    if diff == 0:
                        rt_ground_fits.append(np.dot(ic.rt_norm, undeut_grounds[ic.charge_states[0]].rt_norm))
                    else:
                        if diff > 0:
                            u_buff = copy.copy(undeut_grounds[ic.charge_states[0]].rt_norm)
                            u_buff = np.append(u_buff, [0]*diff)
                            rt_ground_fits.append(np.dot(ic.rt_norm, u_buff))
                        else:
                            rt_ground_fits.append(np.dot(ic.rt_norm, undeut_grounds[ic.charge_states[0]].rt_norm[:diff]) )

                ic.dt_ground_err = np.average(dt_ground_errs)
                ic.rt_ground_err = np.average(rt_ground_errs)
                ic.dt_ground_fit = gmean(dt_ground_fits)
                ic.rt_ground_fit = gmean(rt_ground_fits)

    def generate_sample_paths(self):
        starts = np.linspace(0,0.7,8)
        sample_paths = []
        for start in starts:
            slopes = np.logspace(-0.5,1.8,int(16*(1.0 - start)))
            for slope in slopes:
                sample_paths.append(self.clusters_close_to_line(start, slope))
  
        self.sample_paths = [list(path) for path in set(sample_paths)]

    
    def clusters_close_to_line(
        self, 
        start, 
        slope, 
        undeut_grounds = None, 
        prefiltered_ics = None, 
        max_peak_center = None
        ):

        if undeut_grounds is None:
            undeut_grounds = self.undeut_grounds
        if prefiltered_ics is None:
            prefiltered_ics = self.prefiltered_ics
        if max_peak_center is None:
            max_peak_center = self.max_peak_center

        path = []
        #randomly select UN from PO.undeut_grounds by taking randint within range(len(PO.undeut_grounds.keys()))
        path.append(undeut_grounds[list(undeut_grounds.keys())[np.random.randint(0, high = len(undeut_grounds.keys()))]])
        
        #Relies on use of prefiltered_ics naming convention, when prefilter = 0, all_tp_clusters should have the first n_undeut_runs collapsed into a single list and be named and passed as prefiltered_ics
        xs = np.arange(len(prefiltered_ics))
        expected_centers = 3 + ((start + (1.0-start) * (1-np.exp(-xs * slope / len(xs)))) * (max_peak_center-3) * 0.85)
        
        #prefiltered_ics always has length n_deut_runs+1, all undeuterated are collapsed into PO.prefiltered_ics[0]
        for tp in range(1, len(xs)):
            try:
                peak_dist_rankings = []
                for ic in prefiltered_ics[tp]:
                    if ic.baseline_integrated_mz_com > max_peak_center: continue
                    peak_dist_rankings.append(((abs(ic.baseline_integrated_mz_com-expected_centers[tp])), ic))
                peak_dist_rankings = sorted(peak_dist_rankings, key = lambda x: x[0])
                path.append(peak_dist_rankings[0][1])
            except:
                if len(peak_dist_rankings) == 0:
                    l = len(prefiltered_ics[tp])
                    if l > 0:
                        if l == 1:
                            #pick 0
                            path.append(prefiltered_ics[tp][0])
                        else:
                            #pick lowest mz
                            path.append(sorted(prefiltered_ics[tp], key = lambda ic: ic.baseline_integrated_mz_com)[0])
                    else:
                        #oh shid wutz goin on
                        import os
                        import sys
                        sys.stdout.write("len(PO.prefiltered_ics["+str(tp)+"]) == 0")
                        os.kill(os.getpid())

        
        return tuple(path)
        
    def optimize_paths(self, sample_paths = None, prefiltered_ics = None):
        #Main function of PO, returns the best-scoring HDX IC time-series 'path' of a set of bootstrapped paths. 

        if sample_paths is None:
            sample_paths = self.sample_paths
        if prefiltered_ics is None:
            prefiltered_ics = self.prefiltered_ics

                        
        final_paths = []
        for sample in sample_paths:
            current = copy.deepcopy(sample)
            edited = True
            while edited:

                edited = False
                n_changes = 0
                ic_indices = []
                alt_paths = []

                for tp in range(1,len(current)):
                    for ic in prefiltered_ics[tp]:
                        buffr = copy.copy(current)
                        buffr[tp] = ic
                        alt_paths.append(buffr)

                #Decorate alt_paths
                combo_scoring = []
                for pth in alt_paths:
                    combo_scoring.append(self.combo_score(pth))

                if min(combo_scoring) < self.combo_score(current):
                    current = alt_paths[combo_scoring.index(min(combo_scoring))]
                    n_changes += 1
                    edited = True

                current_score = self.combo_score(current)

                if edited == False:
                    final_paths.append(current)
        final_scores = []
        for pth in final_paths:
            final_scores.append(self.combo_score(pth))

        #This order must be maintained, self.winner must exist before calling find_runners; winner and runners are both needed for set_bokeh tuple
        self.winner = final_paths[final_scores.index(min(final_scores))]
        self.winner_scores = self.report_score(self.winner)
        self.find_runners()
        self.set_bokeh_tuples()
        self.filter_runners()
        self.rt_com_cv = (np.var([ic.rt_com for ic in self.winner])**0.5)/np.mean([ic.rt_com for ic in self.winner])
        self.dt_com_cv = (np.var([np.mean(ic.dt_coms) for ic in self.winner])**0.5)/np.mean([np.mean(ic.dt_coms) for ic in self.winner])
        #Doesn't return, only sets PO attributes

    def find_runners(self):
        # sets self.runners atr. sorts 'runner-up' single substitutions for each tp by score, lower is better.
        winner = self.winner
        prefiltered_ics = self.prefiltered_ics

        runners = []
        for tp in range(len(winner)):
            
            alt_paths = []
            for ic in prefiltered_ics[tp]:
                if ic is not winner[tp]:
                    buffr = copy.copy(winner)
                    buffr[tp] = ic
                    alt_paths.append(buffr)

            combo_scoring = []
            for pth in alt_paths:
                combo_scoring.append(self.combo_score(pth))

            out_buffer = []
            for i in range(len(combo_scoring)):
                min_idx = combo_scoring.index(min(combo_scoring))
                out_buffer.append(alt_paths[min_idx][tp])
                alt_paths.pop(min_idx)
                combo_scoring.pop(min_idx)
            runners.append(out_buffer)

        self.runners = runners

    def set_bokeh_tuples(self):
        #Sets IC.bokeh_tuple to be passed to bokeh for display through the HoverTool
        #Winners store the full values of the winning series scores
        #Runners store the differences between the winning scores and the score if they were to be substituted

        def score_dict(series):
            return {
            "int_mz_std_rmse": self.int_mz_std_rmse(series) * self.int_mz_std_rmse_weight,
            "delta_mz_rate": self.delta_mz_rate(series) * self.delta_mz_rate_weight,
            "dt_ground_rmse": self.dt_ground_rmse(series) * self.dt_ground_rmse_weight, 
            "rt_ground_rmse": self.rt_ground_rmse(series) * self.rt_ground_rmse_weight,
            "dt_ground_fit": self.dt_ground_fit(series) * self.dt_ground_fit_weight,
            "rt_ground_fit": self.rt_ground_fit(series) * self.rt_ground_fit_weight,
            "baseline_peak_error": self.baseline_peak_error(series) * self.baseline_peak_error_weight,
            "auc_ground_rmse": self.auc_ground_rmse(series) * self.auc_ground_rmse_weight
            }

        def score_diff(winner_scores, substituted_scores):
            return (
            winner_scores["int_mz_std_rmse"] - substituted_scores["int_mz_std_rmse"],
            winner_scores["delta_mz_rate"] - substituted_scores["delta_mz_rate"],
            winner_scores["dt_ground_rmse"] - substituted_scores["dt_ground_rmse"],
            winner_scores["rt_ground_rmse"] - substituted_scores["rt_ground_rmse"],
            winner_scores["dt_ground_fit"] - substituted_scores["dt_ground_fit"],
            winner_scores["rt_ground_fit"] - substituted_scores["rt_ground_fit"],
            winner_scores["baseline_peak_error"] - substituted_scores["baseline_peak_error"],
            winner_scores["auc_ground_rmse"] - substituted_scores["auc_ground_rmse"],
            sum([winner_scores[key] for key in winner_scores.keys()]) - sum([substituted_scores[key] for key in substituted_scores.keys()])
            )


        winner = self.winner
        runners = self.runners
        winner_scores = score_dict(winner)

        #Winners store absolute values of scores
        for ic in winner:
            ic.bokeh_tuple = ic.info_tuple + (
                ic.rt_ground_err,
                ic.dt_ground_err,
                winner_scores['int_mz_std_rmse'],
                winner_scores['delta_mz_rate'],
                winner_scores['dt_ground_rmse'],
                winner_scores['rt_ground_rmse'],
                winner_scores['dt_ground_fit'],
                winner_scores['rt_ground_fit'],
                winner_scores['baseline_peak_error'],
                winner_scores['auc_ground_rmse'],
                0
                )

        #Runners store the differences between the winning scores and their scores when substituted
        #Negative values here mean the Runner causes that score to be worse, as all scores are low-is-better and the diffs are calculated winner-substituted
        for tp in range(len(runners)):
            for ic in runners[tp]:
                substituted_series = copy.copy(winner)
                substituted_series[tp] = ic
                substituted_scores = score_dict(substituted_series)
                ic.bokeh_tuple = ic.info_tuple + (ic.rt_ground_err, ic.dt_ground_err) + score_diff(winner_scores, substituted_scores)

    def filter_runners(self, n_runners = 10):
        filtered_runners = []
        for tp in self.runners:
            if len(tp) > n_runners:
                #Sum score diffs and sort index array in descending order by score value
                #Put first n_runners ics into filtered_runners
                tp_scores = [sum(ic.bokeh_tuple[-8:]) for ic in tp]
                tp_idxs = list(range(len(tp)))
                hi_to_lo = sorted(tp_idxs, key=lambda idx: tp_scores[idx], reverse = True)
                filtered_runners.append([tp[idx] for idx in hi_to_lo if hi_to_lo.index(idx) < n_runners])
            else:
                #if len(tp)<n_runners just append tp to filtered_runners
                filtered_runners.append(tp)
        self.filtered_runners = filtered_runners

        
##########################################################################################################################################################################################################################################
### Scoring Functions for PathOptimizer ##################################################################################################################################################################################################
##########################################################################################################################################################################################################################################



    def int_mz_std_rmse(self, ics):
        #calculates the difference in standard deviation from the mean from timepoint i-1 to i for i in [2, len(ics)]
        sd = 0
        for i in range(2, len(ics)):
            sd += (ics[i].baseline_integrated_mz_std-ics[i-1].baseline_integrated_mz_std)**2.0

        return math.sqrt(sd)

    def gabe_delta_mz_rate(self, major_species_centroids, timepoints=None):
        #reproduce logic of delta_mz_rate for Gabe's old data

        if timepoints is None:
            timepoints = self.timepoints

        #take mean of undeuts, fix for gabe's data
        major_species_centroids[0] = np.mean([major_species_centroids.pop(0) for i in range(3)])
        sd = 0
        previous_rate = (major_species_centroids[1] - major_species_centroids[0]) / (timepoints[1] - timepoints[0])
        for i in range(2, len(major_species_centroids)):
            #if previous_rate == 0: diagnostic for /0 error
            new_com = major_species_centroids[i]
            if new_com < major_species_centroids[i-1]: #if we went backwards
                sd += 100*(new_com - major_species_centroids[i-1]) ** 2.0 #penalize for going backwards
                new_com = major_species_centroids[i-1] + 0.01 #pretend we went forwards for calculating current rate
            current_rate = max([(new_com - major_species_centroids[i-1]), 0.1]) / (timepoints[i] - timepoints[i-1])
            if (current_rate / previous_rate) > 1.2:
                sd += (current_rate / previous_rate) ** 2.0
            previous_rate = current_rate
        return sd/len(major_species_centroids)

    def delta_mz_rate(self, ics, timepoints = None):

        if timepoints is None:
            timepoints = self.timepoints

        sd = 0
        previous_rate = (ics[1].baseline_integrated_mz_com - ics[0].baseline_integrated_mz_com) / (timepoints[1] - timepoints[0])
        for i in range(2, len(ics)):
            #if previous_rate == 0: diagnostic for /0 error
            new_com = ics[i].baseline_integrated_mz_com
            if new_com < ics[i-1].baseline_integrated_mz_com: #if we went backwards
                sd += 100*(new_com - ics[i-1].baseline_integrated_mz_com) ** 2.0 #penalize for going backwards
                new_com = ics[i-1].baseline_integrated_mz_com + 0.01 #pretend we went forwards for calculating current rate
            current_rate = max([(new_com - ics[i-1].baseline_integrated_mz_com), 0.1]) / (timepoints[i] - timepoints[i-1])
            if (current_rate / previous_rate) > 1.2:
                sd += (current_rate / previous_rate) ** 2.0
            previous_rate = current_rate
        return sd/len(ics)
    
    def int_mz_rot_fit(self, ics):
        #Compares i to i-1 from ics[2]
        errors = []
        for i in range(2, len(ics)):
            i_mz, j_mz = ics[i].baseline_integrated_mz_norm, ics[i-1].baseline_integrated_mz_norm

            new_indices_i =  np.nonzero(i_mz)[0] - np.argmax(i_mz)
            new_indices_j =  np.nonzero(j_mz)[0] - np.argmax(j_mz)
            
            concat_indices = np.concatenate([new_indices_i, new_indices_j])
            common_low_index = min(concat_indices)
            common_high_index = max(concat_indices)
            
            new_array_i, new_array_j = np.zeros((common_high_index - common_low_index + 1)), np.zeros((common_high_index - common_low_index + 1))
            
            new_indices_i -= common_low_index
            new_indices_j -= common_low_index
             
            new_array_i[new_indices_i] = i_mz[np.nonzero(i_mz)[0]]
            new_array_j[new_indices_j] = j_mz[np.nonzero(j_mz)[0]]
            
            errors.append(np.dot(new_array_i, new_array_j))
        return -np.average(errors)

    def dt_ground_rmse(self, ics): #rmse penalizes strong single outliers, score is minimized - lower is better
        return math.sqrt(sum([ic.dt_ground_err**2 for ic in ics])/len(ics))

    def rt_ground_rmse(self, ics): 
        return math.sqrt(sum([ic.rt_ground_err**2 for ic in ics])/len(ics))

    def dt_ground_fit(self, ics): #This is incorrect, TODO
         return sum([(1.0 / ic.dt_ground_fit) for ic in ics])

    def rt_ground_fit(self, ics): #Both, TODO
        return sum([(1.0 / ic.rt_ground_fit) for ic in ics])

    def baseline_peak_error(self, ics): # Use RMSE instead TODO
        #returns avg of peak_errors from baseline subtracted int_mz -> minimize score
        return np.average([ic.baseline_peak_error for ic in ics])

    def auc_ground_rmse(self, ics, undeut_grounds = None): #TODO make this better
        #find corresponding charge state to each ic, compute AUC error, return avg err
        if undeut_grounds is None:
            undeut_grounds = self.undeut_grounds

        sd = 0
        for ic in ics:
            for key in undeut_grounds.keys():
                if set(undeut_grounds[key].charge_states).intersection(set(ic.charge_states)):
                    sd += (ic.log_baseline_auc-undeut_grounds[key].log_baseline_auc)**2
        return math.sqrt(np.mean(sd))
    
    #Eventually put defaults here as else statements
    def set_score_weights(
        self, 
        int_mz_std_rmse_weight = None, 
        baseline_peak_error_weight = None, 
        delta_mz_rate_weight = None, 
        int_mz_rot_fit_weight = None, 
        rt_ground_fit_weight = None, 
        dt_ground_fit_weight = None, 
        rt_ground_rmse_weight = None, 
        dt_ground_rmse_weight = None, 
        auc_ground_rmse_weight = None
        ):
        
        if int_mz_std_rmse_weight != None:
            self.int_mz_std_rmse_weight = int_mz_std_rmse_weight
        if baseline_peak_error_weight != None:
            self.baseline_peak_error_weight = baseline_peak_error_weight
        if delta_mz_rate_weight != None:
            self.delta_mz_rate_weight = delta_mz_rate_weight
        if int_mz_rot_fit_weight != None:
            self.int_mz_rot_fit_weight = int_mz_rot_fit_weight
        if rt_ground_fit_weight != None: 
            self.rt_ground_fit_weight = rt_ground_fit_weight
        if dt_ground_fit_weight != None: 
            self.dt_ground_fit_weight = dt_ground_fit_weight
        if rt_ground_rmse_weight != None:
            self.rt_ground_rmse_weight = rt_ground_rmse_weight
        if dt_ground_rmse_weight != None:
            self.dt_ground_rmse_weight = dt_ground_rmse_weight
        if auc_ground_rmse_weight != None:
            self.auc_ground_rmse_weight = auc_ground_rmse_weight
    
    def combo_score(self, ics):
        
        return sum([self.int_mz_std_rmse_weight*self.int_mz_std_rmse(ics),
                self.baseline_peak_error_weight*self.baseline_peak_error(ics),
                self.delta_mz_rate_weight*self.delta_mz_rate(ics), 
                #self.int_mz_rot_fit_weight*self.int_mz_rot_fit(ics), 
                self.dt_ground_rmse_weight*self.dt_ground_rmse(ics), 
                self.dt_ground_fit_weight*self.dt_ground_fit(ics), 
                self.rt_ground_fit_weight*self.rt_ground_fit(ics), 
                self.rt_ground_rmse_weight*self.rt_ground_rmse(ics),  
                self.auc_ground_rmse_weight*self.auc_ground_rmse(ics)])

    def report_score(self, ics):
        
        return {"int_mz_std_rmse": (self.int_mz_std_rmse_weight, self.int_mz_std_rmse(ics)),
                "baseline_peak_error":(self.baseline_peak_error_weight, self.baseline_peak_error(ics)),
                "delta_mz_rate":(self.delta_mz_rate_weight, self.delta_mz_rate(ics)), 
                #self.int_mz_rot_fit_weight*self.int_mz_rot_fit(ics), 
                "dt_ground_rmse":(self.dt_ground_rmse_weight, self.dt_ground_rmse(ics)), 
                "dt_ground_fit":(self.dt_ground_fit_weight, self.dt_ground_fit(ics)), 
                "rt_ground_fit":(self.rt_ground_fit_weight, self.rt_ground_fit(ics)), 
                "rt_ground_rmse":(self.rt_ground_rmse_weight, self.rt_ground_rmse(ics)),  
                "auc_ground_rmse":(self.auc_ground_rmse_weight, self.auc_ground_rmse(ics))}

    def bokeh_plot(self, outpath):


        def manual_cmap(value, low, high, palette):
            interval = (high-low)/len(palette)
            n_colors = len(palette)
            if value <= interval:
                return palette[0]
            else:
                if value > (n_colors-1)*interval:
                    return palette[n_colors-1]
                else:
                    for i in range(1, n_colors-2):
                        if value > interval*i and value <= interval*i+1:
                            return palette[i]
                        
        def winner_added_mass_plotter(source, tooltips, old_source = None):
            p = figure(title = 'Winning Timeseries Mean Added-Mass, Colored by RTxDT Error in ms', plot_height = 400, plot_width = 1275, background_fill_color = 'whitesmoke', x_range = (-1, max([int(tp) for tp in source.data['timepoint']])+1), tools = 'pan,wheel_zoom,reset,help', tooltips = tooltips)
            winner_view = CDSView(source = source, filters = [GroupFilter(column_name = "winner_or_runner", group = str(0))])
            err_mapper = linear_cmap(field_name = 'rtxdt_err', palette = Spectral6, low = 0, high = 1)
            color_bar = ColorBar(color_mapper = err_mapper['transform'], width = 10,  location = (0,0))
            
            #Get mean value from source and map value to Spectral6
            mean_rtxdt_err = source.data['rtxdt_err'][0]
            mean_color = manual_cmap(mean_rtxdt_err, 0, 2, Spectral6)
            p.multi_line(xs = 'whisker_x', ys = 'whisker_y', source = source, view = winner_view, line_color = 'black', line_width = 1.5)
            p.line(x = 'timepoint', y = 'baseline_integrated_mz_com', line_color = mean_color, source = source, view = winner_view,  line_width = 3)
            p.circle(x = 'timepoint', y = 'baseline_integrated_mz_com', source = source, view = winner_view, line_color = err_mapper, color = err_mapper, fill_alpha = 1, size = 12)
            #p.add_layout(Whisker(source = source, base = "timepoint", view = winner_view, upper = "upper_added_mass", lower = "lower_added_mass"))

            if old_source is not None: #plot added-masses of all charges of protein
                old_hover = HoverTool(tooltips = [
                    ("Charge", "@charge"),
                    ("Delta MZ Rate Score", "@delta_mz_rate"),
                    ("Fit of Undeuterated Added-Mass Distribution to Theoretical Distribution", "@fit_to_theo_dist")
                    ],
                    names = ['old'])

                old_ics = MultiLine(xs = 'added_mass_xs', ys = 'centroid', line_color = 'wheat', line_width = 1.5)
                old_tp_view = CDSView(source = old_source, filters = [GroupFilter(column_name = 'type', group = 'ts')])
                old_renderer = p.add_glyph(old_source, old_ics, view = old_tp_view, name = 'old')
                #p.add_layout(Whisker(source = old_source, base = "added_mass_xs", upper = "uppers", lower = "lowers")) #TODO Change alt color
                #p.add_tools(old_hover)

            p.xaxis.axis_label = "Timepoint Index"
            p.yaxis.axis_label = "Mean Added-Mass Units"
            p.min_border_top = 100
            p.min_border_left = 100
            p.min_border_right = 100
            p.add_layout(color_bar, 'right')
            return p
        
        
        def winner_rtdt_plotter(source, tooltips):
            #set top margin
            p = figure(title = "Winning Timeseries RT and DT Center-of-Mass Error to Undeuterated Isotopic Cluster", plot_height = 300, plot_width = 1275, x_range = (-3.5, 70), y_range = (-0.1, 2), background_fill_color = 'whitesmoke', tools = 'pan,wheel_zoom,box_zoom,hover,reset,help', tooltips = tooltips)
            winner_view = CDSView(source = source, filters = [GroupFilter(column_name = "winner_or_runner", group = str(0))])
            p.x(x = 'rt_ground_err', y = 'dt_ground_err', source = source, view = winner_view, fill_alpha = 1, size = 5, color = "black")
            p.xaxis.axis_label = "RT Error (ms)"
            p.yaxis.axis_label = "DT Error (ms)"
            p.min_border_left = 100
            p.min_border_right = 100
            
            glyph = Text(x = "rt_ground_err", y = "dt_ground_err", text = "timepoint")
            p.add_glyph(source, glyph, view = winner_view)
            return p

        
        def winner_plotter(source, i, tooltips, old_source = None):
            if i == max([int(tp) for tp in source.data['timepoint']]):
                p = figure(title = "Timepoint "+str(i)+": Winning Isotopic-Cluster Added-Mass Distribution", plot_height = 400, plot_width = 450, y_range = (0, max_intensity), background_fill_color = 'whitesmoke', tools = 'pan,wheel_zoom,hover,reset,help')
                p.min_border_bottom = 100
            else:
                p = figure(title = "Timepoint "+str(i)+": Winning Isotopic Cluster Added-Mass Distribution", plot_height = 300, plot_width = 450, y_range = (0, max_intensity), background_fill_color = 'whitesmoke', tools = 'pan,wheel_zoom,hover,reset,help')
            p.title.text_font_size = "8pt"

            #Have a figure by here, use glyph plotting from here
            new_hover = HoverTool(tooltips = tooltips, names = ['new'])
            index_view = CDSView(source = source, filters = [IndexFilter(indices = [i])])
            new_ics = MultiLine(xs = 'int_mz_x', ys = 'int_mz_y', line_color = "blue", line_width = 1.5)
            new_ics_hover = MultiLine(xs = 'int_mz_x', ys = 'int_mz_y', line_color = "red", line_width = 1.5)
            new_renderer = p.add_glyph(source, new_ics, view = index_view, name = 'new', hover_glyph=new_ics_hover)
            p.add_tools(new_hover)

            if old_source is not None: #plot ics matching the timepoint from old data
                old_hover = HoverTool(tooltips = [("Charge", '@charge'), ("Added-Mass Distribution Centroid", "@"),("Width", '@width')], names = ['old'])
                old_ics = MultiLine(xs = 'int_mz_xs', ys = 'int_mz_ys', line_color = 'wheat', line_width = 1.5)
                old_ics_hover = MultiLine(xs = 'int_mz_xs', ys = 'int_mz_ys', line_color = 'red', line_width = 1.5)
                old_tp_view = CDSView(source = old_source, filters = [GroupFilter(column_name = 'type', group = 'ic'), GroupFilter(column_name = 'timepoint', group = str(i))])
                old_renderer = p.add_glyph(old_source, old_ics, view = old_tp_view, hover_glyph = old_ics_hover, name = 'old')
                p.add_tools(old_hover)

            p.xaxis.axis_label = "Added-Mass Units"
            p.yaxis.axis_label = "Relative Intensity"
            p.min_border_left = 100
            return p


        def runner_plotter(source, i, tooltips):
            if i == max([int(tp) for tp in source.data['timepoint']]):
                p = figure(title = "Runner-Up Isotopic Cluster Added-Mass Distributions", plot_height = 400, plot_width = 375, y_range = (0,max_intensity), background_fill_color = 'whitesmoke', tools = 'pan,wheel_zoom,hover,reset,help', tooltips = tooltips)
                p.min_border_bottom = 100
            else:
                p = figure(title = "Runner-Up Isotopic Cluster Added-Mass Distributions", plot_height = 300, plot_width = 375, y_range = (0,max_intensity), background_fill_color = 'whitesmoke', tools = 'pan,wheel_zoom,hover,reset,help', tooltips = tooltips)
            p.title.text_font_size = "8pt"
            runner_timepoint_view = CDSView(source = source, filters = [GroupFilter(column_name = "timepoint", group = str(i)), GroupFilter(column_name = "winner_or_runner", group = str(1))])
            p.multi_line(xs = 'int_mz_x', ys = 'int_mz_y', source = source, view = runner_timepoint_view, line_color = "blue", alpha = 0.5, hover_color = 'red', hover_alpha = 1, line_width = 1.5)
            p.xaxis.axis_label = "Added-Mass Units"
            p.yaxis.axis_label = "Relative Intensity"
            return p


        def rtdt_plotter(source, i, tooltips):
            if i == max([int(tp) for tp in source.data['timepoint']]):
                p = figure(title = "RT and DT Error from Undeuterated", plot_height = 400, plot_width = 450, background_fill_color = 'whitesmoke', x_range = (-3.5, 70), y_range = (-0.1,2), tools = 'pan,wheel_zoom,hover,reset,help', tooltips = tooltips)
                p.min_border_bottom = 100
            else:
                p = figure(title = "Retention and Drift Center-of-Mass Error to Undeuterated", plot_height = 300, plot_width = 450, background_fill_color = 'whitesmoke', x_range = (-3.5, 70), y_range = (-0.1,2), tools = 'pan,wheel_zoom,hover,reset,help', tooltips = tooltips)
            p.title.text_font_size = "8pt"
            timepoint_runner_view = CDSView(source=source, filters = [GroupFilter(column_name = "timepoint", group = str(i)), GroupFilter(column_name = 'winner_or_runner', group = str(1))])
            p.circle(x = "rt_ground_err", y = 'dt_ground_err', source = source, view = timepoint_runner_view, line_color = 'blue', hover_color = 'red', alpha = 0.25, hover_alpha = 1, size = 5)
            timepoint_winner_view = CDSView(source = source, filters = [GroupFilter(column_name = "timepoint", group = str(i)), GroupFilter(column_name = 'winner_or_runner', group = str(0))])
            p.circle(x = "rt_ground_err", y = 'dt_ground_err', source = source, view = timepoint_winner_view, line_color = 'black', fill_color = "black", hover_color = 'red', size = 5)
            p.xaxis.axis_label = "RT COM Error from Ground (ms)"
            p.yaxis.axis_label = "DT COM Error from Ground (ms)"
            p.min_border_right = 100 
            return p


        output_file(outpath, mode = 'inline')

        #Start old_data source creation
        if self.old_data_dir is not None:

            #create bokeh datasource from gabe's hx_fits

            #divide source by data-types: single-tp ic-level data and charge-state time-series level data
            #old_charges will be used for plotting added-mass and time-series stats

            #init dicts with columns for plotting
            old_ics = dict.fromkeys(['timepoint', 'added_mass_centroid', 'added_mass_width', 'int_mz_ys', 'int_mz_xs', 'type'])
            for key in old_ics.keys():
                old_ics[key] = []

            old_charges = dict.fromkeys(['major_species_integrated_intensities', 'centroid', 'major_species_widths', 'fit_to_theo_dist', 'delta_mz_rate', 'added_mass_xs', 'lowers', 'uppers', 'type', 'charge'])
            for key in old_charges.keys():
                old_charges[key] = []

            #set switch to pull values that only need to be computed once, append each old_file's values to old_data{}
            old_switch = None
            for ts in self.old_data:

                
                int_mz_xs = list(range(len(ts['major_species_integrated_intensities'][0])))
                timepoints = list(range(len(ts['centroid'])))
               

                #Add line to old_charges for each charge file
                for key in ['major_species_integrated_intensities', 'centroid', 'fit_to_theo_dist', 'delta_mz_rate', 'charge']:
                    old_charges[key].append(ts[key])

                old_charges['added_mass_xs'].append(timepoints)
                old_charges['major_species_widths'].append([len(np.nonzero(ic)[0]) for ic in ts['major_species_integrated_intensities']])
                old_charges['lowers'].append([ts['centroid'][tp]-(ts['major_species_widths'][tp]/2) for tp in timepoints])
                old_charges['uppers'].append([ts['centroid'][tp]+(ts['major_species_widths'][tp]/2) for tp in timepoints])
                old_charges['type'].append('ts')

                #Add line to old_ics for each hdx timepoint in each charge file
                for tp in timepoints:
                    old_ics['timepoint'].append(str(tp))
                    old_ics['added_mass_centroid'].append(ts['centroid'][tp])
                    old_ics['added_mass_width'].append(len(np.nonzero(ts['major_species_integrated_intensities'][tp])[0]))
                    old_ics['int_mz_ys'].append(ts['major_species_integrated_intensities'][tp])
                    old_ics['int_mz_xs'].append(int_mz_xs)
                    old_ics['type'].append("ic")

            ts_df = pd.DataFrame.from_dict(old_charges) #len = number of identified charge states for given protein name
            ic_df = pd.DataFrame.from_dict(old_ics) #len = n_charges * n_hdx_timepoints

            self.old_undeut_ground_dot_products = dict.fromkeys(ts_df['charge'].values)
            for charge in self.old_undeut_ground_dot_products.keys():
                self.old_undeut_ground_dot_products[charge] = ts_df.loc[ts_df['charge'] == charge]['fit_to_theo_dist'].values

            old_df = pd.concat([ts_df, ic_df])

            #make cds from df
            gabe_cds = ColumnDataSource(old_df)
        #End old_data source creation

        #TODO: Eventually all the info tuples should be dicts for ease of use in pandas and bokeh, change this style of source construction to dicts of lists
        #This has also just become a clusterfuck and needs to be cleaned up
        data = []
        winner_rtxdt_rmse = np.sqrt(np.mean([((ic.bokeh_tuple[18]*0.07)*ic.bokeh_tuple[19])**2 for ic in self.winner]))
        for tp in range(len(self.winner)):
            edit_buffer = copy.copy(self.winner[tp].bokeh_tuple)
            edit_buffer = edit_buffer[:18]+(edit_buffer[18]*0.07,)+edit_buffer[19:]+(str(tp), np.nonzero(edit_buffer[17])[0][-1], np.nonzero(edit_buffer[17])[0][0], "0", ((edit_buffer[18]*0.07)*edit_buffer[19]), winner_rtxdt_rmse, np.asarray([tp, tp]), np.asarray([np.nonzero(edit_buffer[17])[0][0], np.nonzero(edit_buffer[17])[0][-1]])) #0.07 is adjustment from bins to ms
            data.append(edit_buffer)

        for tp in range(len(self.filtered_runners)):
            for ic in self.filtered_runners[tp]:
                edit_buffer = copy.copy(ic.bokeh_tuple)
                edit_buffer = edit_buffer[:18]+(edit_buffer[18]*0.07,)+edit_buffer[19:]+(str(tp), np.nonzero(edit_buffer[17])[0][-1], np.nonzero(edit_buffer[17])[0][0], "1", ((edit_buffer[18]*0.07)*edit_buffer[19]), "NA", np.asarray([tp, tp]), np.asarray([np.nonzero(edit_buffer[17])[0][0], np.nonzero(edit_buffer[17])[0][-1]])) #0.07 is adjustment from bins to ms
                data.append(edit_buffer)


        source_frame = pd.DataFrame(data, columns = [
                                        "source_file",
                                        "tensor_idx", 
                                        "n_factors", 
                                        "factor_idx", 
                                        "cluster_idx", 
                                        "charge_states",
                                        "n_concatenated",
                                        "mz_bin_low", #from factor bins
                                        "mz_bin_high", #from factor bins
                                        "baseline_subtracted_area_under_curve", 
                                        "baseline_subtracted_grate_sum",  
                                        "baseline_subtracted_peak_error",
                                        "baseline_integrated_mz_com", 
                                        "abs_mz_com",
                                        "rt", 
                                        "dt", 
                                        "int_mz_x",
                                        "int_mz_y",
                                        "rt_ground_err",
                                        "dt_ground_err",
                                        "int_mz_std_rmse",
                                        "delta_mz_rate",
                                        "dt_ground_rmse_score",
                                        "rt_ground_rmse_score",
                                        "dt_ground_fit",
                                        "rt_ground_fit",
                                        "baseline_peak_error",
                                        "auc_ground_rmse",
                                        "net_score_difference",
                                        "timepoint",
                                        "upper_added_mass",
                                        "lower_added_mass",
                                        "winner_or_runner",
                                        "rtxdt_err",
                                        "rtxdt_rmse",
                                        "whisker_x",
                                        "whisker_y"
                                        ]
                                    )
        
        cds = ColumnDataSource(source_frame)

        max_intensity = max([max(int_mz) for int_mz in source_frame['int_mz_y'].values])

        #HoverToolTips, determines information to be displayed when hovering over a glyph
        winner_tts = [
            ("Tensor Index", "@tensor_idx"),
            ("Charge State(s)", "@charge_states"),
            ("Timepoint", "@timepoint"),
            ("Peak Error", "@baseline_subtracted_peak_error"),
            ("Center of Mass in Added-Mass_Units", "@baseline_integrated_mz_com"),
            ("Center of Mass in M/Z", "@abs_mz_com"),
            ("Retention Time COM Error to Ground", "@rt_ground_err"),
            ("Drift Time COM Error to Ground", "@dt_ground_err"),
            ("int_mz_std_rmse", "@int_mz_std_rmse"),
            ("delta_mz_rate", "@delta_mz_rate"),
            ("dt_ground_rmse", "@dt_ground_rmse_score"),
            ("rt_ground_rmse", "@rt_ground_rmse_score"),
            ("dt_ground_fit", "@dt_ground_fit"),
            ("rt_ground_fit", "@rt_ground_fit"),
            ("baseline_peak_error", "@baseline_peak_error"),
            ("auc_ground_rmse", "@auc_ground_rmse"),
            ("Δ scores", "positive is better"),
            ("Δ_net_score", "@net_score_difference")
            ]

        runner_tts = [
            ("Tensor Index", "@tensor_idx"),
            ("Charge State(s)", "@charge_states"),
            ("Timepoint", "@timepoint"),
            ("Peak Error", "@baseline_subtracted_peak_error"),
            ("Center of Mass in Added-Mass_Units", "@baseline_integrated_mz_com"),
            ("Center of Mass in M/Z", "@abs_mz_com"),
            ("Retention Time COM Error to Ground", "@rt_ground_err"),
            ("Drift Time COM Error to Ground", "@dt_ground_err"),
            ("Δ scores", "positive is better"),
            ("Δ int_mz_std_err", "@int_mz_std_err"),
            ("Δ delta_mz_rate", "@delta_mz_rate"),
            ("Δ dt_ground_rmse", "@dt_ground_rmse_score"),
            ("Δ rt_ground_rmse", "@rt_ground_rmse_score"),
            ("Δ dt_ground_fit", "@dt_ground_fit"),
            ("Δ rt_ground_fit", "@rt_ground_fit"),
            ("Δ baseline_peak_error", "@baseline_peak_error"),
            ("Δ auc_ground_rmse", "@auc_ground_rmse"),
            ("Δ_net_score", "@net_score_difference")
            ]
        
        mass_added_tts = [
            ("Timepoint", "@timepoint"),
            ("Charge State(s)", "@charge_states"),
            ("RT COM Error (ms)", "@rt_ground_err"),
            ("DT COM Error (ms)","@dt_ground_err"),
            ("DTxRT Error", "@rtxdt_err")
        ]
        
        winner_rtdt_tts = [
            ("Timepoint","@timepoint"),
            ("Charge State(s)", "@charge_states"),
            ("RT COM Error (ms)", "@rt_ground_err"),
            ("DT COM Error (ms)","@dt_ground_err")
        ]
        
        n_timepoints = len(self.winner)
        if self.old_data_dir is not None:
            winner_plots = [winner_plotter(cds, i, winner_tts, old_source=gabe_cds) for i in range(n_timepoints)]    
        
        else:
            winner_plots = [winner_plotter(cds, i, winner_tts) for i in range(n_timepoints)]  
            
        runner_plots = [runner_plotter(cds, i, runner_tts) for i in range(n_timepoints)]    
        rtdt_plots = [rtdt_plotter(cds, i, runner_tts) for i in range(n_timepoints)]

        rows = []
        if self.old_data_dir is not None:
            rows.append(gridplot([winner_added_mass_plotter(cds, mass_added_tts, old_source=gabe_cds)], sizing_mode='fixed', toolbar_location="left", ncols=1))
            
        else: 
            rows.append(gridplot([winner_added_mass_plotter(cds, mass_added_tts)], sizing_mode='fixed', toolbar_location="left", ncols=1))
        rows.append(gridplot([winner_rtdt_plotter(cds, winner_rtdt_tts)], sizing_mode='fixed', toolbar_location="left", ncols=1))
        [rows.append(gridplot([winner_plots[i], runner_plots[i], rtdt_plots[i]], sizing_mode='fixed', toolbar_location="left", ncols=3)) for i in range(n_timepoints)]

        if self.old_files is not None:
            final = column(Div(text='''<h1 style='margin-left: 300px'>HDX Timeseries Plot for '''+self.name+'''</h1>'''), Div(text="<h3 style='margin-left: 300px'>New Undeuterated-Ground Fits to Theoretical MZ Distribution: "+str(self.undeut_ground_dot_products)+"</h3>"), Div(text="<h3 style='margin-left: 300px'>Old Undeuterated-Ground Fits to Theoretical MZ Distribution: "+str(self.old_undeut_ground_dot_products)+"</h3>"), gridplot(rows, sizing_mode='fixed', toolbar_location=None, ncols=1))
        
        else:
            final = column(Div(text='''<h1 style='margin-left: 300px'>HDX Timeseries Plot for '''+self.name+'''</h1>'''), Div(text="<h3 style='margin-left: 300px'>Undeuterated-Ground Fits to Theoretical MZ Distribution: "+str(self.undeut_ground_dot_products)+"</h3>"), gridplot(rows, sizing_mode='fixed', toolbar_location=None, ncols=1))

        save(final)

        return source_frame, old_df



##########################################################################################################################################################################################################################################
### Plotting Functions ###################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################

#CAN PROBABLY REMOVE THESE? OR KEEP FOR ALTERNATIVE TO BOKEH? TODO

def example_plot(ax, fontsize = 12):
    ax.plot([1, 2])
    ax.locator_params(nbins = 3)
    ax.set_xlabel('x-label', fontsize = fontsize)
    ax.set_ylabel('y-label', fontsize = fontsize)
    ax.set_title('Title', fontsize = fontsize)
    
def winner_plot(ax, ic, olds = None):
    if olds is not None:
        for old in olds:
            ax.plot(old / np.max(old) * np.max(ic.baseline_integrated_mz), color = 'grey')
    ax.plot(ic.baseline_integrated_mz)
    ax.set_xlabel('Added-Mass Units', fontsize = 12)
    ax.set_ylabel('Intensity', fontsize = 12)
    ax.set_title('Winning Timepoint Isotope Cluster', fontsize = 12)

#save runners up within the path_opt fxn, pass tp-wise in the plotting loop
def runner_up_plot(ax, ics):
    for ic in ics:
        ax.plot(ic.baseline_integrated_mz)
    ax.set_xlabel('Added-Mass Units', fontsize = 12)
    ax.set_ylabel('Intensity', fontsize = 12)
    ax.set_title('5 Highest Area Timepoint Isotope Clusters', fontsize = 12)

def winner_rtdt_plot(ax, undeuts, winners, runners):

    """
    Plots rt/dt differences between winning path ICs of PathOptimizer and the undeuterated IC with matching charge
    undeut_grounds = undeut_grounds (self = instance of PathOptimizer)
    winners = PO.winner (PO = instance of PathOptimizer)
    """
    #DT is now X axis, RT is Y
    runner_tp_rts, runner_tp_dts = [[ic.rt_com for ic in tp] for tp in runners], [[ic.dt_coms for ic in tp] for tp in runners]
    win_rts, win_dts = [(ic.rt_com, ic.charge_states) for ic in winners], [(ic.dt_coms, ic.charge_states)  for ic in winners]
    un_rts, un_dts = [undeuts[key].rt_com for key in undeuts.keys()], [undeuts[key].dt_coms for key in undeuts.keys()]
    diffs = []
    for rtup, dtup in zip(win_rts, win_dts):
        if rtup[1] == dtup[1]:
            diffs.append((rtup[0]-undeuts[rtup[1][0]].rt_com, dtup[0]-undeuts[rtup[1][0]].dt_coms[0]))
    
    #for rt_tp, dt_tp in zip(runner_tp_rts, runner_tp_dts):
        #ax.scatter(rt_tp, dt_tp, color='orange', alpha=0.1)
    for ic in diffs:
        #Swapped axes
        ax.scatter(ic[1], ic[0])
    for i in range(len(diffs)):
        #Swapped axes
        ax.annotate(str(i), (diffs[i][1], diffs[i][0]))
    ax.set_ylabel('Retention Time Error Bins', fontsize = 12)
    ax.set_xlabel('Drift Time Error Bins', fontsize = 12)
    ax.set_title('DT_Undeut_Error/RT_Undeut_Error of Winning Isotope Clusters by Timepoint', fontsize = 12)
    ax.set_xlim(-10,10)
    ax.set_ylim(-25,25)
    ax.grid()
    
def runner_up_rtdt_plot(ax, undeuts, winners, runners):
    #DT is now X axis, RT is Y
    win_rts, win_dts = [ic.rt_com for ic in winners], [ic.dt_coms for ic in winners]
    runner_tp_rts, runner_tp_dts = [[(ic.rt_com, ic.charge_states) for ic in tp] for tp in runners], [[(ic.dt_coms, ic.charge_states) for ic in tp] for tp in runners]
    un_rts, un_dts = [undeuts[key].rt_com for key in undeuts.keys()], [undeuts[key].dt_coms for key in undeuts.keys()]
    diffs = []
    for rtp, dtp in zip(runner_tp_rts, runner_tp_dts):
        rt_tp_buffer = []
        dt_tp_buffer = []
        for tup in rtp:
            rt_tp_buffer.append(tup[0]-undeuts[tup[1][0]].rt_com)
        for tup in dtp:
            dt_tp_buffer.append(tup[0]-undeuts[tup[1][0]].dt_coms[0])
        diffs.append((rt_tp_buffer, dt_tp_buffer))
        
    #for rt, dt in zip(win_rts, win_dts):
        #ax.scatter(rt,dt, color='white')
    for rt_tp, dt_tp in diffs:
        #Swapped axes
        #dt here is always a list and may have len>1, rt is always a float
        xs = []
        ys = []
        for ddts, drt in zip(np.asarray(dt_tp).ravel(), rt_tp):
            if type(ddts) == list:
                for ddt in ddts:
                    xs.append(ddt)
                    ys.append(drt)
            else:
                xs.append(ddts)
                ys.append(drt)
        ax.scatter(xs,ys)
    for tp in range(len(diffs)):
        #Swapped axes
        rt_tp, dt_tp = diffs[tp]
        for i in range(len(rt_tp)):
            ax.annotate("tp: "+str(tp)+", #: "+str(i), (dt_tp[i], rt_tp[i]))
    ax.set_ylabel('Retention Time Bins', fontsize = 12)
    ax.set_xlabel('Drift Time Bins', fontsize = 12)
    ax.set_title('DT/RT Error of Runner-Up Isotope Clusters by Timepoint', fontsize = 12)
    ax.set_xlim(-10,10)
    ax.set_ylim(-25,25)
    ax.grid()

def stg2_plot(winners, runners, undeuts, name, path = None, compare_old_path = None):
    
    strip_name = "_".join(name.split('_')[:3])

    if compare_old_path is not None:
        old_files = glob.glob(compare_old_path+"*"+strip_name+"*.pickle")
        gabe_charge_states = []
        for file in old_files:
            with open(file, 'rb') as charge_series:
                gabe_charge_states.append(pickle.load(charge_series))
            
    fig1 = plt.figure(constrained_layout = True, figsize = (30, 40))
    fig1.suptitle(name)
    spec1 = fig1.add_gridspec(ncols = 3, nrows = len(winners))
    for col in range(3):
        for row in range(len(winners)):
            ax = fig1.add_subplot(spec1[row, col])
            if compare_old_path is not None:
                gabe_tp_int_mzs = [cs[1][row+2]['comp_filtered'] for cs in gabe_charge_states]
            if col == 0:
                if compare_old_path is not None:
                    winner_plot(ax, winners[row], gabe_tp_int_mzs)
                else:
                    winner_plot(ax, winners[row])
            else:
                runner_up_plot(ax, runners[row])

    ax1 = fig1.add_subplot(spec1[:len(winners)//2,  2])
    ax2 = fig1.add_subplot(spec1[len(winners)//2:, 2])
    
    winner_rtdt_plot(ax1, undeuts, winners, runners)
    runner_up_rtdt_plot(ax2, undeuts, winners, runners)
    spec1.tight_layout(fig1)
    if path != None:
        fig1.savefig(path, format = "pdf")
    plt.show()
    plt.close()



def limit_read(path):
   return cpickle.loads(zlib.decompress(open(path, 'rb').read()))

def limit_write(obj, outpath):
    with open(outpath, 'wb') as file:
        file.write(zlib.compress(cpickle.dumps(obj)))