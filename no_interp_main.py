#main_simple.py but without interpolation


import numpy as np
import scipy as sp
import math
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from djs_spectra import *


def modified_binned_djs(earth_wave, earth_flux, exo_wave, exo_flux):

	#bins = {'H2O':[0.52, 0.72], 'CO2': [1.42, 1.58], 'O3': [0.94, 0.99], 'CH4':[0.73, 0.81], 'O2':[0.59, 0.70]}
	bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}

	binned_Djs = {'H2O': np.nan, 'CO2': np.nan, 'O3': np.nan, 'CH4':np.nan, 'O2':np.nan}


	for key in bins:
		##truncate spectra
		beg, end = bins[key]

		beg_idx, end_idx =  (np.abs(earth_wave-beg)).argmin(), (np.abs(earth_wave-end)).argmin()

		bin_earth_flux, bin_exo_flux = earth_flux[beg_idx:end_idx], exo_flux[beg_idx:end_idx]


		#int_earth = interp_bins_earth[key]
		#int_exo = interp_bins_exo[key]

		##new, uniform array of wavelengths to use for all spectra
		#wave_new = np.linspace(beg, end, num=num_steps)


		binned_Djs[key] = djs(bin_earth_flux, bin_exo_flux)
	

	return binned_Djs


if __name__ == '__main__':

	#-------------------------------------------------------------
	#
	#               USER PARAMETERS
	#
	#-------------------------------------------------------------

	#what plot(s) do you want to make?

	plt_plain_spectra = False #just plot plain spectra
	plt_binned_interp = True #plot both spectra with interpolation in bins near biomarkers

	#do you want to normalize the flux of all the spectra to the scale of the flux of Earth? (good for Djs calculations)
	normalize = True

	###change these to the location of your Earth and exoplanet spectra files
	earth_infile = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/earth.txt'
	exo_infile = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/HAT-P-1b_spectrum.txt'


	##########You should not need to edit anything below this line
	#-------------------------------------------------------------
	#
	#                IMPORT SPECTRAL DATA
	#
	#-------------------------------------------------------------
	#W/m^2*Ansgstrom, microns
	earth_wave, earth_flux = np.loadtxt(earth_infile, unpack=True, skiprows=2)
	exo_wave, exo_flux = np.loadtxt(exo_infile, unpack=True, skiprows=2)

	###for exo-transmit only!! convert wavelength from m to um 
	earth_wave = earth_wave * 1e6
	exo_wave = exo_wave * 1e6

	if normalize:
		exo_flux = exo_flux*np.mean(earth_flux/exo_flux)

	#wavelength bins for biomarkers from HITRAN
	bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}

	binned_Djs = modified_binned_djs(earth_wave, earth_flux, exo_wave, exo_flux)
	print('Djs near H20 is ', binned_Djs['H2O'], ', near CO2 is ', binned_Djs['CO2'], ', near O3 is ', binned_Djs['O3'], ', near CH4 is ', binned_Djs['CH4'], ', near O2 is ', binned_Djs['O2'])
