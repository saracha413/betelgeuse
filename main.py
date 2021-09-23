import numpy as np
import scipy as sp
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from djs_spectra import *


if __name__ == '__main__':


    #do you want to normalize the flux of all the spectra to the scale of the flux of Earth?
    normalize = True


    #What do you want to plot?
    plt_spectra = False
    plt_djs_dens = True



    #-------------------------------------------------------------
    #
    #                IMPORT SPECTRAL DATA
    #
    #-------------------------------------------------------------
    #W/m^2*Ansgstrom, microns
    earth_wave, earth_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Earth_Lundock081121_Spec_Sun_HiRes.txt', unpack=True)
    mars_wave, mars_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Mars_McCord1971_Spec_Sun_HiRes.txt', unpack=True)
    jupiter_wave, jupiter_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Jupiter_Lundock080507_Spec_Sun_HiRes.txt', unpack=True)

    if normalize:
        avg_flux = np.mean(earth_flux)
        mars_flux = mars_flux * avg_flux/np.mean(mars_flux)
        jupiter_flux = jupiter_flux *avg_flux/np.mean(jupiter_flux)

    binned_Djs_mars, binned_Djs_jupiter = get_binned_Djs(earth_wave, earth_flux, mars_wave, mars_flux, jupiter_wave, jupiter_flux)

    #most common isotopes of CO and H20 from HITRAN
    nu_spec, S_spec, ID_spec, _, _, _, _ , _ = np.loadtxt('6140cfc7.out.txt', unpack = True)


    print('Djs for Mars near H20 is ', binned_Djs_mars['H2O'], ', near CO2 is ', binned_Djs_mars['CO2'], ', near O3 is ', binned_Djs_mars['O3'], ', near CH4 is ', binned_Djs_mars['CH4'], ', near O2 is ', binned_Djs_mars['O2'])
    print('Djs for Jupiter near H20 is ', binned_Djs_jupiter['H2O'], ', near CO2 is ', binned_Djs_jupiter['CO2'], ', near O3 is ', binned_Djs_jupiter['O3'], ', near CH4 is ', binned_Djs_jupiter['CH4'], ', near O2 is ', binned_Djs_jupiter['O2'])


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




    #------------- plot spectra --------------
    if plt_spectra:
        fig, axs = plt.subplots(nrows=3, figsize=(5,10))

        axs[0].plot(earth_wave, earth_flux, label='Earth')
        axs[0].plot(wave_new, int_earth(wave_new), label='interp Earth', linestyle='--')
        axs[1].plot(mars_wave, mars_flux, label='Mars')
        axs[1].plot(wave_new, int_mars(wave_new), label='interp Mars', linestyle='--')
        axs[2].plot(jupiter_wave, jupiter_flux, label='Jupiter')
        axs[2].plot(wave_new, int_jupiter(wave_new), label='interp Jupiter', linestyle='--')

        for ax in axs:
            ax.set_yscale('log')
            ax.set_xlim([0.5,2.5])
            ax.set_ylim([0, 1e-10])
            #ax.set_ylim([-34, -22])
            ax.set_xlabel('Wavelength (um)')
            ax.set_ylabel('Flux (W/(m^2*Angstroms)')
            ax.legend()

        plt.tight_layout()
        #plt.savefig('cleaned_n=1000', format='png')
        plt.show()
        
        
    #------------- plot djs density --------------
    if plt_djs_dens:
        fig, axs = plt.subplots(nrows=3, figsize=(10,10), sharex=True)


        #clean spectra
        earth_wave, earth_flux = clean(earth_wave, earth_flux)
        mars_wave, mars_flux = clean(mars_wave, mars_flux)
        jupiter_wave, jupiter_flux = clean(jupiter_wave, jupiter_flux)

        #interpolate full spectra
        int_earth = interp1d(earth_wave, earth_flux, kind='nearest')
        int_mars = interp1d(mars_wave, mars_flux, kind='nearest')
        int_jupiter = interp1d(jupiter_wave, jupiter_flux, kind='nearest')

        #where to truncate wavelength so all spectra have same number of bins
        start = max(min(earth_wave), min(mars_wave), min(jupiter_wave))
        stop = min(max(earth_wave), max(mars_wave), max(jupiter_wave))

        #new, uniform array of wavelengths to use for all spectra
        wave_new = np.linspace(start, stop, num=1000)



        djs_dens_jup = djs_density(int_earth(wave_new), int_jupiter(wave_new))
        djs_dens_mars = djs_density(int_earth(wave_new), int_mars(wave_new))
        
        
        #angular wavenumnber
        k = 2*np.pi/wave_new
        #linear wavenumber
        nu = 1/wave_new
        
        matplotlib.rcParams['text.usetex'] = True

        axs[0].axvline(x=0.6, ymin=0, ymax=1, label=r'H$_2$O', color='red', linestyle='dashed') #H2O
        axs[0].axvline(x=1.47, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
        axs[0].axvline(x=0.96, ymin=0, ymax=1, label=r'O$_3$', color='blue', linestyle='dashed')
        axs[0].axvline(x=0.75, ymin=0, ymax=1, label=r'CH$_4$', color='purple', linestyle='dashed')
        axs[0].axvline(x=0.64, ymin=0, ymax=1, label=r'O$_2$', color='black', linestyle='dashed')

        axs[0].plot(wave_new, int_mars(wave_new), color='black')

        axs[0].scatter(mars_wave, mars_flux, s=5)  
        #axs[0].set_ylim([0, 0.4e-11])
        axs[0].set_ylabel(r'Flux (W/(m$^2 \AA$)')
        axs[0].legend()
        axs[0].set_title('Mars spectrum')

        
        #axs[1].vlines(x=[0.6], ymin=3.7e-11, ymax=4e-11, label='H2O', color='green') #H2O
        axs[1].axvline(x=0.6, ymin=0, ymax=1, label=r'H$_2$O', color='red', linestyle='dashed')
        #axs[1].text(0.6, 3.23e-11, 'H2O' )
        #axs[1].vlines(x=[0.42, 1.47], ymin=3.7e-11, ymax=4e-11, label='CO2', color='red') #CO2
        #axs[1].axvline(x=0.42, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
        axs[1].axvline(x=1.47, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
        #axs[1].text(0.42, 2e-11, 'CO2' )
        #axs[1].text(1.47, 2.9e-12, 'CO2' )
        #axs[1].vlines(x=[0.96], ymin=1.7e-12, ymax=4e-11, label='O3', color='blue') #O3
        axs[1].axvline(x=0.96, ymin=0, ymax=1, label=r'O$_3$', color='blue', linestyle='dashed')
        #axs[1].text(0.96, 1.93e-11, 'O3')
        #axs[1].vlines(x=[0.75], ymin=3.7e-11, ymax=4e-11, label='CH4', color='purple') #CH4
        axs[1].axvline(x=0.75, ymin=0, ymax=1, label=r'CH$_4$', color='purple', linestyle='dashed')
        #axs[0].text(0.75, 1.25e-13, 'CH4')
        #axs[1].text(0.75, 1.41e-11, 'CH4')
        #axs[1].vlines(x=[0.64], ymin=3.7e-11, ymax=4e-11, label='O2', color='black') #O2
        axs[1].axvline(x=0.64, ymin=0, ymax=1, label=r'O$_2$', color='black', linestyle='dashed')
        #axs[1].text(0.64, 2.5e-11, 'O2')



        axs[1].plot(wave_new, int_earth(wave_new),color='black')
        axs[1].scatter(earth_wave, earth_flux, s=5)
        #axs[1].set_ylim([0, 5e-11])
        axs[1].set_ylabel(r'Flux (W/(m$^2 \AA$)')
        axs[1].set_title('Earth spectrum')
        #axs[1].legend()


        
        axs[2].axvline(x=0.6, ymin=0, ymax=1, label=r'H$_2$O', color='red', linestyle='dashed')
        axs[2].axvline(x=1.47, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
        axs[2].axvline(x=0.96, ymin=0, ymax=1, label=r'O$_3$', color='blue', linestyle='dashed')
        axs[2].axvline(x=0.75, ymin=0, ymax=1, label=r'CH$_4$', color='purple', linestyle='dashed')
        axs[2].axvline(x=0.64, ymin=0, ymax=1, label=r'O$_2$', color='black', linestyle='dashed')

        axs[2].scatter(wave_new, djs_dens_mars, s=10)
        axs[2].set_ylim([1e-14, 0.3])
        axs[2].set_yscale('log')
        axs[2].set_xlabel(r'Wavelength ($\mu$ m)')
        axs[2].set_ylabel(r'$\mathcal{D}_{JS}$(Earth|Mars)')
        axs[2].set_title('Information difference between Earth and Mars spectra')


        #print JSD of interpolation function for Earth vs. interpolation function for Mars
        print('Djs of Earth vs. Mars is ', djs(int_earth(wave_new),int_mars(wave_new)))
        print('Djs of Earth vs. Jupiter is ', djs(int_earth(wave_new),int_jupiter(wave_new)))

        print('Djs density at Mars H2O is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.6)[1])
        #print('Djs density at first CO2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.42)[1])
        print('Djs density at Mars CO2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 1.47)[1])
        print('Djs density at Mars O3 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.96)[1])
        print('Djs density at Mars CH4 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.75)[1])
        print('Djs density at Mars O2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.64)[1])


        print('Djs density at Jupiter H2O is ', find_nearest_Djs(wave_new, djs_dens_jup, 0.6)[1])
        #print('Djs density at first CO2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.42)[1])
        print('Djs density at Jupiter CO2 is ', find_nearest_Djs(wave_new, djs_dens_jup, 1.47)[1])
        print('Djs density at Jupiter O3 is ', find_nearest_Djs(wave_new, djs_dens_jup, 0.96)[1])
        print('Djs density at Jupiter CH4 is ', find_nearest_Djs(wave_new, djs_dens_jup, 0.75)[1])
        print('Djs density at Jupiter O2 is ', find_nearest_Djs(wave_new, djs_dens_jup, 0.64)[1])

        
        plt.tight_layout()
        #plt.savefig('Djs_density_vs_lambda.png', format='png')
        plt.show()

