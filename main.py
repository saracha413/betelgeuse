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
    normalize = True

    #is this Exo_Transmit data? (false = Sagan Institute data)
    isExoTrans = False

    #What do you want to plot?
    plt_binned_djs = True
    plt_spectra = False
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

    if normalize:
        avg_flux = np.mean(spectra['earth'][1])
        for key in spectra:
            #normalize flux to Earth flux
            spectra[key][1] = spectra[key][1] * avg_flux/np.mean(spectra[key][1])


    bins = {'H2O':[0.5, 0.7], 'CO2': [1.4, 1.6], 'O3': [0.9, 1.0], 'CH4':[0.7, 0.8], 'O2':[0.6, 0.70]}
    
    binned_Djs = spectra
    #for key in spectra:
    #    binned_Djs[key] = bins 
    #    for bin_key in bins:
    #        binned_Djs[key][bin_key] =
     
    for key in spectra:
        wave = spectra[key][0]
        flux = spectra[key][1]
        binned_Djs[key] = get_binned_Djs(earth_wave, earth_flux, wave, flux)

    #binned_Djs_mars, binned_Djs_jupiter = get_binned_Djs(earth_wave, earth_flux, mars_wave, mars_flux, jupiter_wave, jupiter_flux)

    #most common isotopes of CO and H20 from HITRAN
    nu_spec, S_spec, ID_spec, _, _, _, _ , _ = np.loadtxt('6140cfc7.out.txt', unpack = True)

    print(binned_Djs['earth'])
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





