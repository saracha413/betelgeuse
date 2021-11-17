
import numpy as np
import scipy as sp
import math
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from djs_spectra import *
from collections import deque

import os
from os import listdir
from os.path import isfile, join
#find the files in directory
#from https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
path = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/hot_jupiters/'

files = [f for f in listdir(path) if isfile(join(path, f))]

#where to save Djs values to 
newFileName = '/Users/saravannah/AstroSpec/' + 'Djs-density-earth-vs-Jupiter.txt'


base_infile = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/earth.txt'


#add skiprows=2 argument for raw Exo_Transmit data
base_wav, base_flux = np.loadtxt(base_infile, unpack=True, skiprows=2)

#convert to um
base_wav = base_wav * 1e6

#wavelength bins for biomarkers from HITRAN
bins = {'H2O':[0.5, 0.7], 'CO2': [0.42, 0.44], 'CO2_high': [1.4,1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70], 'O2_high':[4.38, 4.85]}


#if file already exists, remove it 
if os.path.exists(newFileName):
	os.remove(newFileName)




djs_dens_bins = {}
for fileName in files:
	if fileName != '.DS_Store':
		exo_wav, exo_flux = np.loadtxt(path+fileName, unpack=True)
		exo_wav = exo_wav * 1e6
		exo_name = fileName.split('_')[0]


		binned_Djs_density = get_binned_Djs_density(base_wav, base_flux, exo_wav, exo_flux, bins)

		binned_wave = get_binned_wave(base_wav, base_flux, exo_wav, exo_flux, bins)

		wavelengths, Djs_dens = np.concatenate(list(binned_wave.values())),  np.concatenate(list(binned_Djs_density.values()))#np.array(binned_wave.values()).flatten(), list(binned_Djs_density.values())

		djs_dens_bins[fileName] = [wavelengths, Djs_dens]

		#for key in binned_wave:


		#	plt.scatter(binned_wave[key], binned_Djs_density[key], label=fileName)#, color=colors[X[key]])

#for key in djs_dens_bins:
#	plt.scatter(djs_dens_bins[key][0], djs_dens_bins[key][1], label=key, s=0.7)
plt.scatter(djs_dens_bins['HAT-P-1b_spectrum2.txt'][0],djs_dens_bins['HAT-P-1b_spectrum2.txt'][1], s=0.9)

wavs = {'H2O': 0.596, 'CO2':0.43, 'O3': 0.955, 'CH4': 0.752, 'O2': 0.643, 'CO2high': 1.471, 'O2high':4.5}

wavs_list = list(wavs.values())

for key in wavs:
	plt.annotate(key, (wavs[key],1e-9))


plt.xlabel('wavelength (um)')
plt.ylabel('Djs density')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

#np.savetext()


