
import numpy as np
import scipy as sp
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from djs_spectra import *

if __name__ == '__main__':


    #-------------------------------------------------------------
    #
    #               USER PARAMETERS
    #
    #-------------------------------------------------------------

    #what plot(s) do you want to make?

    plt_plain_spectra = True #just plot plain spectra
    plt_binned_interp = False #plot both spectra with interpolation in bins near biomarkers

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

    binned_Djs = get_binned_Djs(earth_wave, earth_flux, exo_wave, exo_flux)
    print('Djs near H20 is ', binned_Djs['H2O'], ', near CO2 is ', binned_Djs['CO2'], ', near O3 is ', binned_Djs['O3'], ', near CH4 is ', binned_Djs['CH4'], ', near O2 is ', binned_Djs['O2'])



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

       





    if plt_binned_interp:

        interp_bins_earth, interp_bins_exo = get_binned_interp(earth_wave, earth_flux), get_binned_interp(exo_wave, exo_flux)
        iso_colors = {'H2O':'red', 'CO2': 'green', 'O3': 'blue', 'CH4':'purple', 'O2':'gray'}

        fig, axs = plt.subplots(nrows=2, figsize=(5,7), sharex=True)

        axs[0].scatter(earth_wave, earth_flux, s=5)
        axs[0].set_title('Earth')
        for key in interp_bins_earth:
            beg, end = bins[key]
            wave_new = np.linspace(beg, end, num=50)
            axs[0].plot(wave_new, interp_bins_earth[key](wave_new), linestyle='dashed', color='black')
            axs[0].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])


        axs[1].scatter(exo_wave, exo_flux, s=5)
        axs[1].set_title('Exoplanet')
        for key in interp_bins_exo:
            beg, end = bins[key]
            wave_new = np.linspace(beg, end, num=50)
            axs[1].plot(wave_new, interp_bins_exo[key](wave_new), linestyle='dashed', color='black')
            axs[1].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])

        for ax in axs:
            #ax.set_ylim([-1e-11, 6e-11])
            #ax.set_xscale('log')
            ax.set_xlim([0.1,2])
            ax.legend(loc='upper left')
            ax.set_ylabel('Transit depth (percent)')


        plt.xlabel('Wavelength (um)')
        plt.tight_layout()
        plt.show()
