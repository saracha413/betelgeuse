
import numpy as np
import scipy as sp
import math
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from djs_spectra import *
from collections import deque


from os import listdir
from os.path import isfile, join
#find the files in directory
#from https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
path = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/hot_jupiters/'

files = [f for f in listdir(path) if isfile(join(path, f))]

#where to save Djs values to 
newFileName = path + 'Djs-hotJup-vs-Jupiter_v2.txt'


base_infile = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/Jupiter_spectrum_v2.txt'


#add skiprows=2 argument for raw Exo_Transmit data
base_wav, base_flux = np.loadtxt(base_infile, unpack=True)

#convert to um
base_wav = base_wav * 1e6

#wavelength bins for biomarkers from HITRAN
bins = {'H2O':[0.5, 0.7], 'CO2': [0.42, 0.44], 'CO2_high': [1.4,1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70], 'O2_high':[4.38, 4.85]}

data = []

for fileName in files:
	if fileName != '.DS_Store':
		with open(newFileName, 'a') as file:
			exo_wav, exo_flux = np.loadtxt(path+fileName, unpack=True)
			exo_wav = exo_wav * 1e6
			exo_name = fileName.split('_')[0]

			binned_Djs = map(str, get_binned_Djs(base_wav, base_flux, exo_wav, exo_flux, bins).values())

			line = ','.join(binned_Djs)
			line = exo_name+','+line

			#binned_Djs = str(binned_Djs)[1:-1]

			#line = [exo_name, binned_Djs]
			#line = str(line)[1:-1]



			file.write(str(line))
			#np.savetxt()
			file.write('\n')

		#data.append(line)
file.close()

#print('data is ', data)
#np.savetxt(newFileName, data)
