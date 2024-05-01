# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 06:53:51 2023

@author: vmysorea
"""
import sys
sys.path.append('C:/Users/vmysorea/Documents/mne-python/')
import warnings
from matplotlib import pyplot as plt
from scipy import io
import numpy as np
from scipy.stats import sem
from mne.viz import centers_to_edges

plt.switch_backend('QT5Agg')  # Making the plots interactive (Scrollable)
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
# Defining the dimensions and quality of figures
plt.rcParams["figure.figsize"] = (5.5,5)
plt.rcParams['figure.dpi'] = 120
#%%Setting up stuff
save_loc='C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/FinalThesis/ForPPT/'
data_loc = 'D:/PhD/Data/Chin_Data/AnalyzedABR_Gamma/'

# subjlist_y = ['Q412', 'Q422', 'Q424', 'Q426', 'Q428']

# subjlist_m = ['Q351', 'Q363', 'Q364', 'Q365', 'Q368']

# subjlist_tts = ['Q402', 'Q404', 'Q406', 'Q407']

#%% Loading files for ABR frequency analysis across the groups  

#Young Chins
for subj in range(len(subjlist_y)):
    sub = subjlist_y [subj]
    dat = io.loadmat(data_loc + sub + '_ABRpower.mat', squeeze_me=True)
    dat.keys()
    freqs = dat['freqs']
    n_cycles = dat['n_cycles']
    picks=dat['picks']
    a = (dat['power']).mean(axis=0)
    b = (dat['psd']).mean(axis=0)
    t=dat['t_full']

c = freqs.size
d = t.size
power_y = np.zeros((len(subjlist_y),c, d))

for subj in range(len(subjlist_y)):
    sub = subjlist_y [subj]
    power_y[subj, :] = a

power_y_all=(power_y.mean(axis=0))*1e6
power_y_sem = sem(power_y)

#MNH Chins

for subj in range(len(subjlist_m)):
    sub = subjlist_m [subj]
    dat = io.loadmat(data_loc + sub + '_ABRpower.mat', squeeze_me=True)
    dat.keys()
    freqs = dat['freqs']
    n_cycles = dat['n_cycles']
    picks=dat['picks']
    a = (dat['power']).mean(axis=0)
    b = (dat['psd']).mean(axis=0)
    t=dat['t_full']

c = freqs.size
d = t.size
power_m = np.zeros((len(subjlist_m),c, d))

for subj in range(len(subjlist_m)):
    sub = subjlist_m [subj]
    power_m[subj, :] = a

power_m_all=(power_m.mean(axis=0))
power_m_sem = sem(power_m)

#TTS Chins
for subj in range(len(subjlist_tts)):
    sub = subjlist_tts [subj]
    dat = io.loadmat(data_loc + sub + '_ABRpower.mat', squeeze_me=True)
    dat.keys()
    freqs = dat['freqs']
    n_cycles = dat['n_cycles']
    picks=dat['picks']
    a = (dat['power']).mean(axis=0)
    b = (dat['psd']).mean(axis=0)
    t=dat['t_full']

c = freqs.size
d = t.size
power_tts = np.zeros((len(subjlist_tts),c, d))

for subj in range(len(subjlist_tts)):
    sub = subjlist_tts [subj]
    power_tts[subj, :] = a

power_tts_all=(power_tts.mean(axis=0))
power_tts_sem = sem(power_tts)

#%%Plotting power across the groups 
vmin=-0.5
vmax=0.5

fig, ax = plt.subplots()
x, y = centers_to_edges(t, freqs)
mesh = ax.pcolormesh(x, y, power_tts_all, cmap='RdBu_r', vmin=vmin, vmax=vmax)
ax.set_title("Gamma NB (N=" + str(len(subjlist_y)) + ")", y=1.03)
ax.set(ylim=freqs[[0, -1]], xlabel="Time (s)")
fig.colorbar(mesh)
plt.tight_layout()
plt.show()

# plt.savefig(save_loc + 'ONH_highGamma_NB.png', dpi=500)

#Trial whatever -- Plot in spectrum 

plt.plot(freqs, power_y_all.mean(axis=1), label = 'YNH')
plt.plot(freqs, power_m_all.mean(axis=1), label = 'MNH')
plt.plot(freqs, power_tts_all.mean(axis=1), label = 'TTS')
plt.xlim(0,40)
plt.legend()
plt.show()

#%% Looking across the sedation conditions 

save_loc='C:/Users/vmysorea/Desktop/PhD/Stim_Analysis/MTB_Analysis/FinalThesis/'
data_loc = 'D:/PhD/Data/Chin_Data/AnalyzedABR_Gamma/'

subjlist_awake = ['Q419', 'Q415', 'Q414']

subjlist_anes = ['Q419', 'Q415', 'Q414']

subjlist_ls = ['Q351', 'Q363', 'Q364', 'Q365', 'Q368', 'Q402', 'Q404', 'Q406', 'Q407', 'Q412', 'Q422', 'Q424', 'Q426', 'Q428']

#%% Loading files for ABR frequency analysis across the age groups  

#Young Chins
for subj in range(len(subjlist_awake)):
    sub = subjlist_awake [subj]
    dat = io.loadmat(data_loc + sub + '_Awake_ABRpower.mat', squeeze_me=True)
    dat.keys()
    freqs = dat['freqs']
    n_cycles = dat['n_cycles']
    picks=dat['picks']
    a = (dat['power']).mean(axis=0)
    b = (dat['psd']).mean(axis=0)
    t=dat['t_full']

c = freqs.size
d = t.size
power_awake = np.zeros((len(subjlist_awake),c, d))

for subj in range(len(subjlist_awake)):
    sub = subjlist_awake [subj]
    power_awake[subj, :] = a

power_awake_all=(power_awake.mean(axis=0))*1e6
power_awake_sem = sem(power_awake)

#MNH Chins

for subj in range(len(subjlist_anes)):
    sub = subjlist_anes [subj]
    dat = io.loadmat(data_loc + sub + '_Anesthetized_ABRpower.mat', squeeze_me=True)
    dat.keys()
    freqs = dat['freqs']
    n_cycles = dat['n_cycles']
    picks=dat['picks']
    a = (dat['power']).mean(axis=0)
    b = (dat['psd']).mean(axis=0)
    t=dat['t_full']

c = freqs.size
d = t.size
power_anes = np.zeros((len(subjlist_anes),c, d))

for subj in range(len(subjlist_anes)):
    sub = subjlist_anes [subj]
    power_anes[subj, :] = a

power_anes_all=(power_anes.mean(axis=0))
power_anes_sem = sem(power_anes)

# TTS Chins

for subj in range(len(subjlist_ls)):
    sub = subjlist_ls [subj]
    dat = io.loadmat(data_loc + sub + '_ABRpower.mat', squeeze_me=True)
    dat.keys()
    freqs = dat['freqs']
    n_cycles = dat['n_cycles']
    picks=dat['picks']
    a = (dat['power']).mean(axis=0)
    b = (dat['psd']).mean(axis=0)
    t=dat['t_full']

c = freqs.size
d = t.size
power_ls = np.zeros((len(subjlist_ls),c, d))

for subj in range(len(subjlist_ls)):
    sub = subjlist_ls [subj]
    power_ls[subj, :] = a

power_ls_all=(power_ls.mean(axis=0))
power_ls_sem = sem(power_ls)

###Plotting power across the groups 
vmin=-0.5
vmax=0.5

# fig, ax = plt.subplots()
# x, y = centers_to_edges(t, freqs)
# mesh = ax.pcolormesh(x, y, power_ls_all, cmap='RdBu_r', vmin=vmin, vmax=vmax)
# ax.set_title("Frequency (N=" + str(len(subjlist_ls)) + ")", y=1.03)
# ax.set(ylim=freqs[[0, -1]], xlabel="Time (s)")
# fig.colorbar(mesh)
# plt.tight_layout()
# plt.show()

# plt.savefig(save_loc + 'ONH_highGamma_NB.png', dpi=500)

###Trial whatever -- Plot in spectrum 
mean_awake = power_awake_all.mean(axis=1)
sem_awake = power_awake_sem.mean(axis=1)

mean_anes = power_anes_all.mean(axis=1)*1e5
sem_anes = power_anes_sem.mean(axis=1)*1e5

mean_ls = power_ls_all.mean(axis=1)*1e5
sem_ls = power_ls_sem.mean(axis=1)*1e5

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
# ax.plot(freqs, mean_awake, label = 'Awake (N=3)', color= '#1f77b4')
# ax.fill_between(freqs, (mean_awake - sem_awake), (mean_awake + sem_awake), color='#1f77b4', alpha=0.3)
ax.plot(freqs, mean_ls, label = 'Light Sedation (N=15)', color= '#ff7f0e')
ax.fill_between(freqs, (mean_ls - sem_ls), (mean_ls + sem_ls), color= '#ff7f0e', alpha=0.3)
ax.plot(freqs, mean_anes, label = 'Anesthetized (N=3)', color= '#2ca02c')
ax.fill_between(freqs, (mean_anes - sem_anes), (mean_anes + sem_anes), color= '#2ca02c', alpha=0.3)

ax.fill_between(x=(0,4), y1=-0, y2=0.02,color= 'limegreen',alpha=0.3)
ax.fill_between(x=(8,12), y1=-0, y2=0.02,color= 'slategray',alpha=0.3)
ax.fill_between(x=(12,20), y1=-0, y2=0.02,color= 'indianred',alpha=0.3)
ax.fill_between(x=(30,90), y1=-0, y2=0.02,color= 'lightsteelblue',alpha=0.6)
ax.text(1.9,0.0175, r'$\delta$', fontsize=18, weight='bold')
ax.text(9,0.0175, r'$\alpha$', fontsize=18, weight='bold')
ax.text(15,0.0175, r'$\beta$', fontsize=18, weight='bold')
ax.text(50,0.0175, r'$\gamma$', fontsize=18, weight='bold')
# plt.vlines(x=(4, 8, 12, 20, 30, 90),ymin=-3.5, ymax=3, color='black', linestyle='--', alpha=1)
# plt.vlines(x=1, ymin=-3.5, ymax=3,color='blue', linestyle='--', alpha=1)
# ax.text(0, 3.1, 'Stim On', va='center', ha='center', fontsize = 12, weight='bold')
# ax.text(1, 3.1, 'Gap', va='center', ha='center', fontsize = 11, weight='bold')
# ax.text(2, 3.1, 'Stim End', va='center', ha='center', fontsize = 12, weight='bold')
plt.xlabel ("Frequency (Hz)", fontsize =16)
plt.ylabel ("Power (\u03bcV$^2$ / Hz)", fontsize=16)
plt.ylim(0, 0.02)
plt.xlim (0,90)
plt.legend(loc='upper right', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

plt.savefig(save_loc + "FreqSpectrum_LSAnes.png", dpi=500, bbox_inches="tight", transparent=True)

