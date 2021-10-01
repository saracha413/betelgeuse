import numpy as np
import scipy as sp
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt

#calculate JSD for two spectra
def djs(spec1, spec2):

	#take Fourier transform of spectra
	pow1 = fft(spec1)
	pow2 = fft(spec2)
 
	#get modal fractions
	p, q = pow1*np.conj(pow1)/sum(pow1*np.conj(pow1)), pow2*np.conj(pow2)/sum(pow2*np.conj(pow2))
	p,q = np.real(p), np.real(q) #force Python to drop the "+0j" part of the "complex" number
	
	#get D_JS of spectra
	r = 1/2 * (p+q)
	p,r,q = np.ma.array(p), np.ma.array(r), np.ma.array(q)
	Djs = 1/2 * np.sum(p*np.log(p/r)) + 1/2 * np.sum(q*np.log(q/r))
 
	return Djs
	
 
#calculate JSD density for two spectra
def djs_density(spec1, spec2):
   pow1 = fft(spec1)
   pow2 = fft(spec2)

   #get modal fractions, |power|^2/sum(|power|^2)
   p, q = pow1*np.conj(pow1)/sum(pow1*np.conj(pow1)), pow2*np.conj(pow2)/sum(pow2*np.conj(pow2))
   p,q = np.real(p), np.real(q) #force Python to drop the "+0j" part
   
   #get D_JS density of spectra
   r = 1/2 * (p+q)
   djs_density = 1/2 * p*np.log(p/r) + 1/2 * q*np.log(q/r)
   
   return djs_density
 
 
	
	
#remove unphysical absorption and emission lines
def clean(wave, flux):
	new_wave, new_flux = [], []
	for i in range(len(wave)-1):
		if abs(flux[i] - flux[i+1]) < 0.5e-10:
			new_flux.append(flux[i])
			new_wave.append(wave[i])
		
	new_wave, new_flux = np.array(new_wave), np.array(new_flux)

	return new_wave, new_flux


#get a
def find_nearest(arr, value):

	idx =  (np.abs(arr-value)).argmin()

	return arr[idx]


#adapted from https://www.codegrepper.com/code-examples/python/python+find+nearest+point+in+array
#get Djs density to match a wavelength value
def find_nearest_Djs(wave_arr, djs_arr, wave_value):

	idx = (np.abs(wave_arr-wave_value)).argmin()

	return wave_arr[idx], djs_arr[idx]


#truncates wavelength and flux arrays from the beginning and ending wavelengths in a wavelength packet
def truncate(beg_wave, end_wave, wave_arr, flux_arr):

	beg_idx, end_idx = (np.abs(wave_arr-beg_wave)).argmin(), (np.abs(wave_arr-end_wave)).argmin()
	wave_arr, flux_arr = wave_arr[beg_idx:end_idx], flux_arr[beg_idx:end_idx]

	return wave_arr, flux_arr


	
#-------------------------------------------------------------
#
#                IMPORT SPECTRAL DATA
#
#-------------------------------------------------------------
#W/m^2*Ansgstrom, microns
#earth_wave, earth_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Earth_Lundock081121_Spec_Sun_HiRes.txt', unpack=True)
#mars_wave, mars_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Mars_McCord1971_Spec_Sun_HiRes.txt', unpack=True)
#jupiter_wave, jupiter_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Jupiter_Lundock080507_Spec_Sun_HiRes.txt', unpack=True)


#earth_wave, earth_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/R8resolution/Sun/Earth_Lundock081121_Spec_Sun_LoRes.txt', unpack=True)
#mars_wave, mars_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/R8resolution/Sun/Mars_McCord1971_Spec_Sun_LoRes.txt', unpack=True)
#jupiter_wave, jupiter_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/R8resolution/Sun/Jupiter_Lundock080507_Spec_Sun_LoRes.txt', unpack=True)


def get_binned_interp(earth_wave, earth_flux, mars_wave, mars_flux, jupiter_wave, jupiter_flux):

	#clean spectra
	earth_wave, earth_flux = clean(earth_wave, earth_flux)
	mars_wave, mars_flux = clean(mars_wave, mars_flux)
	jupiter_wave, jupiter_flux = clean(jupiter_wave, jupiter_flux)


	#bins = {'H2O':[0.52, 0.72], 'CO2': [1.42, 1.58], 'O3': [0.94, 0.99], 'CH4':[0.73, 0.81], 'O2':[0.59, 0.70]}
	bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}

	interp_bins_earth =  {'H2O':[], 'CO2': [], 'O3': [], 'CH4':[], 'O2':[]}
	interp_bins_mars =  {'H2O':[], 'CO2': [], 'O3': [], 'CH4':[], 'O2':[]}
	interp_bins_jupiter =  {'H2O':[], 'CO2': [], 'O3': [], 'CH4':[], 'O2':[]}

	for key in bins:
		#truncate spectra
		beg, end = bins[key]
		earth_wave_trunc, earth_flux_trunc = truncate(beg, end, earth_wave, earth_flux)
		mars_wave_trunc, mars_flux_trunc = truncate(beg, end, mars_wave, mars_flux)
		jupiter_wave_trunc, jupiter_flux_trunc = truncate(beg, end, jupiter_wave, jupiter_flux)


		#interpolate spectra (this returns an interpolation function that's a function of wavelength)
		interp_bins_earth[key] = interp1d(earth_wave_trunc, earth_flux_trunc, kind='nearest', fill_value="extrapolate")
		interp_bins_mars[key] = interp1d(mars_wave_trunc, mars_flux_trunc, kind='nearest', fill_value="extrapolate")
		interp_bins_jupiter[key] = interp1d(jupiter_wave_trunc, jupiter_flux_trunc, kind='nearest', fill_value="extrapolate")



	return interp_bins_earth, interp_bins_mars, interp_bins_jupiter


def get_binned_interp(wave, flux, sagan=False):

	#clean spectra ---- TURN ON FOR SAGAN INSTITUTE DATA or data with unphysical emission/absorption lines
	if sagan:
		wave, flux = clean(earth_wave, earth_flux)


	#bins = {'H2O':[0.52, 0.72], 'CO2': [1.42, 1.58], 'O3': [0.94, 0.99], 'CH4':[0.73, 0.81], 'O2':[0.59, 0.70]}
	bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}

	interp_bins =  {'H2O':[], 'CO2': [], 'O3': [], 'CH4':[], 'O2':[]}

	for key in bins:
		#truncate spectra
		beg, end = bins[key]
		wave_trunc, flux_trunc = truncate(beg, end, wave, flux)

		#interpolate spectra (this returns an interpolation function that's a function of wavelength)
		interp_bins[key] = interp1d(wave_trunc, flux_trunc, kind='nearest', fill_value="extrapolate")


	return interp_bins



def get_binned_Djs(earth_wave, earth_flux, exo_wave, exo_flux):


	#bins = {'H2O':[0.52, 0.72], 'CO2': [1.42, 1.58], 'O3': [0.94, 0.99], 'CH4':[0.73, 0.81], 'O2':[0.59, 0.70]}
	bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}

	binned_Djs = {'H2O': np.nan, 'CO2': np.nan, 'O3': np.nan, 'CH4':np.nan, 'O2':np.nan}


	interp_bins_earth = get_binned_interp(earth_wave, earth_flux)
	interp_bins_exo = get_binned_interp(exo_wave, exo_flux)


	for key in bins:
		##truncate spectra
		beg, end = bins[key]


		int_earth = interp_bins_earth[key]
		int_exo = interp_bins_exo[key]

		#new, uniform array of wavelengths to use for all spectra
		wave_new = np.linspace(beg, end, num=50)


		binned_Djs[key] = djs(int_earth(wave_new),int_exo(wave_new))
	

	return binned_Djs


##where to truncate wavelength so all spectra have same number of bins
#start = max(min(earth_wave), min(mars_wave), min(jupiter_wave))
#stop = min(max(earth_wave), max(mars_wave), max(jupiter_wave))

















