
import numpy as np
import scipy as sp
import math
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from djs_spectra import *



def get_entropy_bias(waves, flux):

    nu = fft(waves)
    #bias = math.log((math.e),2)*(np.amax(nu)-1)) / (2*sum(flux))
    bias = math.log((math.e),2)*(np.amax(nu[1:])-np.amin(nu[1:])) / (2*sum(flux))

    return bias




#truncates wavelength and flux arrays from the beginning and ending wavelengths in a wavelength packet
def get_bin_width(beg_wave, end_wave, wave_arr, flux_arr):

    beg_idx, end_idx = (np.abs(wave_arr-beg_wave)).argmin(), (np.abs(wave_arr-end_wave)).argmin()
    wave_arr, flux_arr = wave_arr[beg_idx:end_idx], flux_arr[beg_idx:end_idx]


    return wave_arr, flux_arr

if __name__ == '__main__':



    #-------------------------------------------------------------
    #
    #               USER PARAMETERS
    #
    #-------------------------------------------------------------

    #what plot(s) do you want to make?

    plt_plain_spectra = False #just plot plain spectra
    plt_binned = True #plot both spectra with interpolation in bins near biomarkers

	#do you want to normalize the flux of all the spectra to the scale of the flux of Earth? (good for Djs calculations)
    normalize = True

    ###change these to the location of your Earth and exoplanet spectra files
    earth_infile = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/earth.txt'
    exo_infile = '/Users/saravannah/AstroSpec/Exo_Transmit_spectra/Jupiter_spectrum.txt'



    ##########You should not need to edit anything below this line
	#-------------------------------------------------------------
    #
    #                IMPORT SPECTRAL DATA
    #
    #-------------------------------------------------------------
    #W/m^2*Ansgstrom, microns
    earth_wave, earth_flux = np.loadtxt(earth_infile, unpack=True, skiprows=2)
    exo_wave, exo_flux = np.loadtxt(exo_infile, unpack=True) #np.loadtxt(exo_infile, unpack=True, skiprows=2)

    ###for exo-transmit only!! convert wavelength from m to um 
    earth_wave = earth_wave * 1e6
    exo_wave = exo_wave * 1e6

    if normalize:
        exo_flux = exo_flux*np.mean(earth_flux/exo_flux)

    #wavelength bins for biomarkers from HITRAN
    bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'CO2_high': [1.4,1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70], 'O2_high':[4.38, 4.85]}

    binned_Djs = get_binned_Djs(earth_wave, earth_flux, exo_wave, exo_flux, bins)
    print('Djs near H20 is ', binned_Djs['H2O'], ', near CO2 is ', binned_Djs['CO2'], ', near O3 is ', binned_Djs['O3'], ', near CH4 is ', binned_Djs['CH4'], ', near O2 is ', binned_Djs['O2'], ' near the high wavelength CO2 transition is ', binned_Djs['CO2_high'], ' near the high wavelength O2 transition is ', binned_Djs[
        'O2_high'])



    #-------------------------------------------------------------
    #
    #                PLOT RESULTS
    #
    #-------------------------------------------------------------

    if plt_plain_spectra:

        plt.scatter(earth_wave, earth_flux, s=5, label='Earth')
        plt.scatter(exo_wave, exo_flux, s=5, label='Exoplanet')

        plt.ylabel('Transit depth (percent)')
        plt.xlabel('Wavelength (um)')
        plt.legend()
        plt.show()

       





    if plt_binned:

        #interp_bins_earth, interp_bins_exo = get_binned_interp(earth_wave, earth_flux), get_binned_interp(exo_wave, exo_flux)
        iso_colors = {'H2O':'red', 'CO2': 'green', 'CO2_high': 'green', 'O3': 'blue', 'CH4':'purple', 'O2':'gray', 'O2_high': 'gray'}

        fig, axs = plt.subplots(nrows=2, figsize=(5,7), sharex=True)

        axs[0].scatter(earth_wave, earth_flux, s=5)
        axs[0].set_title('Earth')
        for key in bins:
            beg, end = bins[key]
            axs[0].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])


        axs[1].scatter(exo_wave, exo_flux, s=5)
        axs[1].set_title('Exoplanet')
        for key in bins:
            beg, end = bins[key]
            axs[1].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])

        for ax in axs:
            #ax.set_ylim([-1e-11, 6e-11])
            ax.set_xscale('log')
            ax.set_xlim([0.3,10])
            ax.legend(loc='upper left')
            ax.set_ylabel('Transit depth (percent)')


        plt.xlabel('Wavelength (um)')
        plt.tight_layout()
        plt.show()
