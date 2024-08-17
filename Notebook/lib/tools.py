import numpy as np
import scipy as scipy
from scipy.stats import gaussian_kde
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import matplotlib.pyplot as plt
from lib.cytoreactors.analysis.analysis import ReactorData
from IPython.display import HTML
import json
import analysis
from matplotlib.pyplot import figure
import matplotlib
from datetime import datetime, timedelta
from tkinter import Tk
from tkinter import filedialog
from collections import defaultdict
import pickle
import fcsparsing
from sklearn.neighbors import KernelDensity
import seaborn as sns
import os, sys, stat
from pathlib import Path
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import least_squares as lsq
from scipy.integrate import solve_ivp

#This function takes in a single well data frame from a time point data frame dict to apply a gating metric to each cell and outputs the updated
#data frame

def compute_ks_FSC_SSC_gating_metric(orig_data, grid_size=100):
    data = orig_data.copy()
    FSC_lin = data['FSC-HLin']
    SSC_lin = data['SSC-HLin']
    I_pos_fsc_ssc = np.logical_and(FSC_lin > 0, SSC_lin > 0)
    FSC_log = np.log10(FSC_lin[I_pos_fsc_ssc])
    SSC_log = np.log10(SSC_lin[I_pos_fsc_ssc])
    kernel = gaussian_kde([FSC_log, SSC_log])
    extent= [np.min(FSC_log),
             np.max(FSC_log),
             np.min(SSC_log),
             np.max(SSC_log)]
    FSC_grid, SSC_grid = np.mgrid[extent[0]:extent[1]:grid_size*1j,
                                  extent[2]:extent[3]:grid_size*1j]
    positions = np.vstack([FSC_grid.ravel(), SSC_grid.ravel()])
    P = np.reshape(kernel(positions).T, FSC_grid.shape)
    FSC_bin = np.floor((FSC_log-extent[0])*grid_size/(extent[1]-extent[0]))
    SSC_bin = np.floor((SSC_log-extent[2])*grid_size/(extent[3]-extent[2]))
    FSC_bin[FSC_bin == grid_size] = grid_size - 1
    SSC_bin[SSC_bin == grid_size] = grid_size - 1
    P_max = np.max(P)
    gating_metric = np.array([P[fsc_bin, ssc_bin] for fsc_bin,ssc_bin in zip(FSC_bin.astype(int),SSC_bin.astype(int))]) / P_max
    data['gating-metric'] = 0.0
    data['gating-metric'][I_pos_fsc_ssc] = gating_metric
    return data
#This function takes in a single well data frame from a time point data frame dict to apply a doublet metric to each cell and outputs the updated
    #data frame
def compute_doublet_metric(orig_data):
    data = orig_data.copy()
    fsch_to_fsca_fit = np.polyfit(data['FSC-HLin'], data['FSC-ALin'], 1)
    pred_fsca = np.polyval(p=fsch_to_fsca_fit, x=data['FSC-HLin'])
    deviation = pred_fsca - data['FSC-ALin']
    std_deviation = np.std(deviation)
    data['doublet-metric'] = np.abs(deviation) / std_deviation
    return data

def infer_FP_amounts(data, AF, fps):
    # build fluo channels names
    ch_colors = list(fps['mCerulean'].index)
    ch_color_names = [ch for ch in ch_colors]
    # build the matrix with with fluorescence vector as lines and cells as columns
    signal = (data[ch_color_names] - AF[ch_color_names]).transpose()
    # build the FP signature matrix
    A = pd.concat([fps[fp] for fp in fps], axis=1)
    # solve the linear system
    result = np.linalg.lstsq(A,signal,rcond=None)
    # store result in the data frame and return
    N = result[0]
    R = result[1]
    for i,fp in enumerate(fps):
        data[fp] = N[i]
    data['deconv_residual'] = R
    return data
def deconvolve_turbi(data,AF,FP_signatures):
    Y2 = infer_FP_amounts(data, AF, FP_signatures)

    return(Y2)

def do_gating(data):

    cols = list()
    cols = list(data.cells.columns)
    cols.append('gating-metric')
    cols.append('doublet-metric')
    gated_cells = pd.DataFrame(columns=cols)
    gk = data.cells.groupby('rel_time_hrs')
    t= np.unique(np.array(data.cells['rel_time_hrs']))
    for jj in range(np.size(t)):
        lk = gk.get_group(t[jj])
        tlk = get_singlets(lk)
        gated_cells=gated_cells.append(tlk)

    data.cells =  pd.merge(data.cells, gated_cells)
    print('done')
    return(data.cells)

def get_singlets(data):
    data = compute_ks_FSC_SSC_gating_metric(data, grid_size=100)
    data = compute_doublet_metric(data)
    return(data)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#Load fluorescence signatures and AF values for deconvolution

with open('FP_signatures_turbi_chosen.pickle', 'rb') as f:   #Load signatures
    FP_signatures = pickle.load(f)

with open('AF_chosen.pickle', 'rb') as f:    #Load Autofluorescence
    AF = pickle.load(f)

#Plotting Parameters

params = {'legend.fontsize': 10,
          'legend.loc' :'lower right',
          'legend.markerscale'   : 1,
          'axes.labelsize' : 10,
          'xtick.labelsize' : 10,
          'ytick.labelsize' : 10,
          'axes.titlesize' : 12,

          'font.family':  'sans-serif',
          'font.style':   'normal',
          'font.variant': 'normal',
          'font.weight':  600,
          'font.stretch': 'normal',
          'font.size':    10.0,
          }

plt.rcParams.update(params)

cgfont = {'fontname':'Century gothic'}
