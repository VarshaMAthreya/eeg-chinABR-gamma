# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 04:27:58 2023

@author: vmysorea
"""
import sys
sys.path.append('C:/Users/vmysorea/Documents/mne-python/')
sys.path.append('C:/Users/vmysorea/Documents/ANLffr/')
import warnings
import mne
import numpy as np
from anlffr.helper import biosemi2mne as bs
from matplotlib import pyplot as plt
import os
import fnmatch
from scipy.io import savemat
from mne.time_frequency import tfr_multitaper
# from itertools import zip_longest
# from scipy.stats import sem 

plt.switch_backend('QT5Agg')  # Making the plots interactive (Scrollable)
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
# Defining the dimensions and quality of figures
plt.rcParams['figure.dpi'] = 120
plt.rcParams["figure.figsize"] = (5.5, 5)
plt.tight_layout

# %% Loading subjects, reading data, mark bad channels
froot = 'D:/PhD/Data/Chin_Data/Awake/'  # file location
save_loc = 'D:/PhD/Data/Chin_Data/AnalyzedABR_Gamma/'

#Q410 has only 3 RR
#Q419 and Q427 doesn't have ABR saved

subjlist = ['Q417']

for subj in subjlist:
    # Load data and read event channel
    fpath = froot + subj + '/'
    bdfs = fnmatch.filter(os.listdir(fpath), subj +'_Awake_ABR.bdf')
    
    print('LOADING! ' + subj +' raw data')

    # Load data and read event channel
    rawlist = []
    evelist = []

    for k, rawname in enumerate(bdfs):
        rawtemp, evestemp = bs.importbdf(fpath + rawname, verbose='DEBUG',refchans=['EXG1', 'EXG2'])
        rawlist += [rawtemp, ]
        evelist += [evestemp, ]
    raw, eves = mne.concatenate_raws(rawlist, events_list=evelist)
    raw.set_channel_types({'EXG3':'eeg'})           #Mastoid -34
    raw.set_channel_types({'EXG4':'eeg'})           #Vertex -35
    raw.set_channel_types({'EXG5':'eeg'})           #Ground -36
    raw.info['bads'].append('EXG6') 
    raw.info['bads'].append('EXG7')
    raw.info['bads'].append('EXG8')
    raw.info['bads'].append('A12') 
    raw.info['bads'].append('A26') 
    raw.info['bads'].append('A27') 
    raw.info['bads'].append('A20') 
    raw.info['bads'].append('A17')
    raw.info['bads'].append('A1')
    raw.info['bads'].append('EXG5')
    raw.info['bads'].append('EXG3')
    raw.info['bads'].append('EXG4')
    # raw.info['bads'].append('A32')
    
    # raw, eves = raw.resample(4096, events=eves)

    # To check and mark bad channels
    # raw.plot(duration=25.0, n_channels=41, scalings=dict(eeg=100e-6))
    
#%% Bad channels for each subject
  
    if subj == ['Q364', 'Q404']:
       raw.info['bads'].append('A13')
       raw.info['bads'].append('A19')
              
    if subj == ['Q412', 'Q424']:
       raw.info['bads'].append('A1')
       
    if subj == ['Q422']:
        raw.info['bads'].append('EXG5')
        raw.info['bads'].append('A21')
        
    if subj == ['Q426']:
        raw.info['bads'].append('A5')
        raw.info['bads'].append('A31')
               
# %% Filtering
    raw.filter(0.1, 90.)
    raw.info
    # raw.plot(duration=25.0, n_channels=41, scalings=dict(eeg=100e-6))

    raw.notch_filter(np.arange(60, 241, 60), filter_length='auto', phase='zero')
    
#%% Creating manual events for frequency analysis -- 5 second duration
    fs = raw.info['sfreq']
    a = np.int64((raw._data.shape[1]/fs) * 5)      #To identify the number of 5 second events required 
    eves_manual = np.zeros((a, 3))                              #Creating manual events here 
    
    eves_manual [:, 2] = 2      #Setting the third column with trigger value 2 (random, just want it to be not 1)
    
    for k in range(a):
        eves_manual[k, 0] = (k+1)*5*fs     #Creating 5 second rows 
    
    eves_manual = np.int64(eves_manual)

#%% Epoching 

    epochs = mne.Epochs(raw, eves_manual, event_id=[2], baseline=None, proj=True,tmin=0, tmax=5, reject=dict(eeg=200e-6), preload=False)
    t_full = epochs.times
    evoked = epochs.average()
    
    # evoked.plot()
    
#%%% Power and ITC Analysis 
    freqs = np.arange(0.1, 90., 2.)
    n_cycles = freqs * 0.2
    
    # epochs_i = epochs.copy().subtract_evoked()
    
    picks = (6, 7, 8,21, 22, 23, 28, 29, 13)
    power = tfr_multitaper(epochs, freqs, n_cycles, 
                           time_bandwidth=4.0, n_jobs=-1, return_itc=False)
    
    # power.plot([2],mode='mean',title='Power', vmin=-2*1e-6, vmax=2*1e-6)
    
    # itc.plot(title='Gap duration of 16 ms - Intertrial Coherence (' + subj + ')',  baseline=(-0.1,0), combine='mean')
   
    # c = (power.data.mean(axis=0))*1e8
    # plt.plot(freqs, c.mean(axis=1))
    # plt.xlim([0,80])
    # plt.show()
    
    psd = epochs.compute_psd(fmin=0.1, fmax=80.0)
    
    # psd = epochs.compute_psd(fmin=0.1, fmax=80.0).plot(average=True, picks=picks)
    
    mat_ids = dict (power = power.data, picks=picks, freqs=freqs, n_cycles=n_cycles, t_full=t_full, psd = psd._data)
    savemat(save_loc + subj + '_Awake_ABRpower.mat', mat_ids)
    
    print('WOOOHOOOO! Saved ' + subj)
    
    del (epochs, evoked, power)
