import numpy as np
import scipy as sp
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
import inquirer
from djs_spectra import *


if __name__ == '__main__':

    #moreFiles = True
    #fileList = []

    #while moreFiles:
    #
    #    inFile = input("Enter file location for first spectrum: ")
    #    fileList.append(inFile)
    #    moreFiles = bool(int(input('Add more spectra? 1 for yes, 0 for no: ')))


    #normalize = bool(int(input('Would you like to normalize the spectral flux? 1 for yes, 0 for no: ')))

    #do you want to normalize the flux of all the spectra to the scale of the flux of Earth?
    normalize = False

    #is this Exo_Transmit data? (false = Sagan Institute data)
    isExoTrans = False

    #What do you want to plot?
    plt_binned_djs = False
    plt_plain_spectra = True
    plt_djs_dens = False




    #


    #-------------------------------------------------------------
    #
    #                IMPORT SPECTRAL DATA
    #
    #-------------------------------------------------------------
    #W/m^2*Ansgstrom, microns
    earth_wave, earth_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Earth_Lundock081121_Spec_Sun_HiRes.txt', unpack=True)
    mars_wave, mars_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Mars_McCord1971_Spec_Sun_HiRes.txt', unpack=True)
    jupiter_wave, jupiter_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Jupiter_Lundock080507_Spec_Sun_HiRes.txt', unpack=True)

    spectra = {'earth': [earth_wave, earth_flux], 'mars': [mars_wave, mars_flux], 'jupiter': [jupiter_wave, jupiter_flux]}

    earth_wave, earth_flux = spectra['earth'][0], spectra['earth'][1]


    if isExoTrans:
        for key in spectra:
            #for exo-transmit only!! convert wavelength from m to um 
            spectra[key][0] = spectra[key][0]*1e6
    else:
        #clean spectra
        for key in spectra:
            spectra[key][0], spectra[key][1] = clean(spectra[key][0], spectra[key][1])
            #spectra[key][0], spectra[key][1] = clean(spectra[key][0], spectra[key][1])
            #spectra = spectra[key][0], spectra[key][1]

    if normalize:
        avg_flux = np.mean(earth_flux)
        for key in spectra:
            #normalize flux to Earth flux
            flux = spectra[key][1]
            foo = flux * avg_flux/np.mean(flux)
            #print('before, spec[key][1] is ', spectra[key][1])
            spectra[key][1] = foo
            #print('after, spec[key][1] is ', spectra[key][1])
            #print('mean of foo minus mean of Earth is ', np.mean(spectra[key][1])-np.mean(earth_flux))


    bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}
    binned_Djs = spectra.copy() #force Python to make a copy of spectra so editing binned_Djs doesn't edit spectra elements 

    for key in spectra:
        wave = spectra[key][0]
        flux = spectra[key][1]
        binned_Djs[key] = get_binned_Djs(earth_wave, earth_flux, wave, flux)



    #binned_Djs_mars, binned_Djs_jupiter = get_binned_Djs(earth_wave, earth_flux, mars_wave, mars_flux, jupiter_wave, jupiter_flux)

    #most common isotopes of CO and H20 from HITRAN
    nu_spec, S_spec, ID_spec, _, _, _, _ , _ = np.loadtxt('6140cfc7.out.txt', unpack = True)


    for key in spectra:
        print('Djs for ', key ,' near H2O is ', binned_Djs[key]['H2O'], ', near CO2 is ', binned_Djs[key]['CO2'], ', near O3 is ', binned_Djs[key]['O3'], ', near CH4 is ', binned_Djs[key]['CH4'], ', near O2 is ', binned_Djs[key]['O2'])
    #print('Djs for Mars near H2O is ', binned_Djs_mars['H2O'], ', near CO2 is ', binned_Djs_mars['CO2'], ', near O3 is ', binned_Djs_mars['O3'], ', near CH4 is ', binned_Djs_mars['CH4'], ', near O2 is ', binned_Djs_mars['O2'])
    #print('Djs for Jupiter near H2O is ', binned_Djs_jupiter['H2O'], ', near CO2 is ', binned_Djs_jupiter['CO2'], ', near O3 is ', binned_Djs_jupiter['O3'], ', near CH4 is ', binned_Djs_jupiter['CH4'], ', near O2 is ', binned_Djs_jupiter['O2'])

    #-------------------------------------------------------------
    #
    #                GET HITRAN SPECTRA
    #
    #-------------------------------------------------------------

    nu_CO, S_CO = [], []
    nu_H2O, S_H2O = [], []

    for i in range(len(nu_spec)):
        if ID_spec[i] == 1: 
            nu_CO.append(nu_spec[i]/1e3)
            S_CO.append(S_spec[i])
        elif ID_spec[i] == 5:
            nu_H2O.append(nu_spec[i]/1e3)
            S_H2O.append(S_spec[i])
        else:
            print(ID_spec[i])


    #-------------------------------------------------------------
    #
    #                PLOT RESULTS
    #
    #-------------------------------------------------------------


    if plt_plain_spectra:


        #plt.scatter(spectra['earth'][0], spectra['earth'][1], label='Earth', s=5)#earth_wave, earth_flux, s=5, label='Earth')
        #plt.scatter(spectra['mars'][0], spectra['mars'][1], s=5, label='Mars')
        #plt.scatter(spectra['jupiter'][0], spectra['jupiter'][1], s=5, label='Jupiter')

        plt.plot(spectra['earth'][0], spectra['earth'][1], label='Earth')#earth_wave, earth_flux, s=5, label='Earth')
        plt.plot(spectra['mars'][0], spectra['mars'][1], label='Mars')
        plt.plot(spectra['jupiter'][0], spectra['jupiter'][1], label='Jupiter')

        #TO-DO: change based on dataset
        plt.ylim([0, 5e-11])
        #plt.ylabel('Transit depth (percent)')
        plt.ylabel('Flux')
        plt.xlabel('Wavelength (um)')
        plt.legend()
        plt.show()



    if plt_binned_djs:

        interp_bins_earth, interp_bins_mars, interp_bins_jupiter = get_binned_interp(spectra['earth'][0], spectra['earth'][1]), get_binned_interp(spectra['mars'][0], spectra['mars'][1]), get_binned_interp(spectra['jupiter'][0], spectra['jupiter'][1])
        iso_colors = {'H2O':'red', 'CO2': 'green', 'O3': 'blue', 'CH4':'purple', 'O2':'gray'}

        fig, axs = plt.subplots(nrows=len(spectra.keys()), figsize=(5,10), sharex=True)

        axs[0].scatter(spectra['earth'][0], spectra['earth'][1], s=5)
        axs[0].set_title('Earth')
        for key in interp_bins_earth:
            beg, end = bins[key]
            wave_new = np.linspace(beg, end, num=50)
            axs[0].plot(wave_new, interp_bins_earth[key](wave_new), linestyle='dashed', color='black')
            axs[0].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])

        axs[1].scatter(spectra['mars'][0], spectra['mars'][1], s=5)
        axs[1].set_title('Mars')
        for key in interp_bins_mars:
                beg, end = bins[key]
                wave_new = np.linspace(beg, end, num=50)
                axs[1].plot(wave_new, interp_bins_mars[key](wave_new), linestyle='dashed', color='black')
                axs[1].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])


        axs[2].scatter(spectra['jupiter'][0], spectra['jupiter'][1], s=5)
        axs[2].set_title('Jupiter')
        for key in interp_bins_jupiter:
            beg, end = bins[key]
            wave_new = np.linspace(beg, end, num=50)
            axs[2].plot(wave_new, interp_bins_jupiter[key](wave_new), linestyle='dashed', color='black')
            axs[2].axvspan(beg, end, alpha=0.5, label=key, color=iso_colors[key])
        

        for ax in axs:
            #ax.set_ylim([-1e-11, 6e-11])
            #ax.set_xscale('log')
            ax.set_xlim([0.1,2])
            ax.legend(loc='upper right')
            ax.set_ylabel('Flux')
            ax.set_ylim([0, 5e-11])


        plt.xlabel('Wavelength (um)')
        plt.tight_layout()
        plt.show()




