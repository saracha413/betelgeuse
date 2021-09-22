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
	Djs = 1/2 * sum(p*np.log(p/r)) + 1/2 * sum(q*np.log(q/r))
 
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




#adapted from https://www.codegrepper.com/code-examples/python/python+find+nearest+point+in+array
#get Djs density to match a wavelength value
def find_nearest_Djs(wav_arr, djs_arr, wav_value):

    idx = (np.abs(wav_arr-wav_value)).argmin()

    return wav_arr[idx], djs_arr[idx]
    
#-------------------------------------------------------------
#
#                IMPORT SPECTRAL DATA
#
#-------------------------------------------------------------
#W/m^2*Ansgstrom, microns
earth_wave, earth_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Earth_Lundock081121_Spec_Sun_HiRes.txt', unpack=True)
mars_wave, mars_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Mars_McCord1971_Spec_Sun_HiRes.txt', unpack=True)
jupiter_wave, jupiter_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/NativeResolution/Sun/Jupiter_Lundock080507_Spec_Sun_HiRes.txt', unpack=True)

print('Mars has ', len(mars_wave), ' points while Jupiter has ', len(jupiter_wave),'.' )

#earth_wave, earth_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/R8resolution/Sun/Earth_Lundock081121_Spec_Sun_LoRes.txt', unpack=True)
#mars_wave, mars_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/R8resolution/Sun/Mars_McCord1971_Spec_Sun_LoRes.txt', unpack=True)
#jupiter_wave, jupiter_flux = np.loadtxt('CatalogofSolarSystemObjects/Spectra/R8resolution/Sun/Jupiter_Lundock080507_Spec_Sun_LoRes.txt', unpack=True)


#clean spectra
earth_wave, earth_flux = clean(earth_wave, earth_flux)
mars_wave, mars_flux = clean(mars_wave, mars_flux)
jupiter_wave, jupiter_flux = clean(jupiter_wave, jupiter_flux)

#interpolate spectra (this returns an interpolation function that's a function of wavelength)
int_earth = interp1d(earth_wave, earth_flux, kind='nearest')
int_mars = interp1d(mars_wave, mars_flux, kind='nearest')
int_jupiter = interp1d(jupiter_wave, jupiter_flux, kind='nearest')


#where to truncate wavelength so all spectra have same number of bins
start = max(min(earth_wave), min(mars_wave), min(jupiter_wave))
stop = min(max(earth_wave), max(mars_wave), max(jupiter_wave))

#new, uniform array of wavelengths to use for all spectra
wave_new = np.linspace(start, stop, num=1000)



#most common isotopes of CO and H20 from HITRAN
nu_spec, S_spec, ID_spec, _, _, _, _ , _ = np.loadtxt('6140cfc7.out.txt', unpack = True)

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

#What do you want to plot?
plt_spectra = False
plt_djs_dens = True


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
    plt.savefig('cleaned_n=1000', format='png')
    plt.show()
    
    
#------------- plot djs density --------------
if plt_djs_dens:
    fig, axs = plt.subplots(nrows=3, figsize=(10,10), sharex=True)
    
    djs_dens_jup = djs_density(int_earth(wave_new), int_jupiter(wave_new))
    djs_dens_mars = djs_density(int_earth(wave_new), int_mars(wave_new))



    
    
    #angular wavenumnber
    k = 2*np.pi/wave_new
    #linear wavenumber
    nu = 1/wave_new
    
    matplotlib.rcParams['text.usetex'] = True

    
    #axs[0].vlines(x=[0.6], ymin=1.7e-12, ymax=1.8e-12, label='H2O', color='green') #H2O
    axs[0].axvline(x=0.6, ymin=0, ymax=1, label=r'H$_2$O', color='red', linestyle='dashed') #H2O
    #axs[0].text(0.6, 1.4e-12, 'H2O' )
    #axs[0].vlines(x=[0.42, 1.47], ymin=1.8e-12, ymax=2e-12, label='CO2', color='red') #CO2
    #axs[0].axvline(x=0.42, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
    axs[0].axvline(x=1.47, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
    #axs[0].text(0.52, 5e-13, 'CO2' )
    #axs[0].text(1.47, 2.8e-13, 'CO2' )
    #axs[0].vlines(x=[0.96], ymin=1.7e-12, ymax=1.8e-12, label='O3', color='blue') #O3
    axs[0].axvline(x=0.96, ymin=0, ymax=1, label=r'O$_3$', color='blue', linestyle='dashed')
    #axs[0].text(0.96, 2.72e-13, 'O3')
    #axs[0].vlines(x=[0.75], ymin=1.7e-12, ymax=1.8e-12, label='CH4', color='purple') #CH4
    axs[0].axvline(x=0.75, ymin=0, ymax=1, label=r'CH$_4$', color='purple', linestyle='dashed')
    #axs[0].text(0.75, 1.25e-13, 'CH4')
    #axs[0].vlines(x=[0.64], ymin=1.7e-12, ymax=1.8e-12, label='O2', color='black') #O2
    axs[0].axvline(x=0.64, ymin=0, ymax=1, label=r'O$_2$', color='black', linestyle='dashed')
    #axs[0].text(0.64, 1.49e-12, 'O2')

    axs[0].plot(wave_new, int_mars(wave_new), color='black')

    axs[0].scatter(mars_wave, mars_flux, s=5)  
    axs[0].set_ylim([0, 0.4e-11])
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
    axs[1].set_ylim([0, 5e-11])
    axs[1].set_ylabel(r'Flux (W/(m$^2 \AA$)')
    axs[1].set_title('Earth spectrum')
    #axs[1].legend()


    
    #axs[2].vlines(x=[0.6], ymin=7e-3, ymax=1e-2, label='H2O', color='green') #H2O
    axs[2].axvline(x=0.6, ymin=0, ymax=1, label=r'H$_2$O', color='red', linestyle='dashed')
    #axs[2].text(0.6, 2e-5, 'H2O' )
    #axs[2].vlines(x=[0.42, 1.47, 1.5, 1.53], ymin=7e-3, ymax=1e-2, label='CO2', color='red') #CO2
    #axs[2].axvline(x=0.42, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
    axs[2].axvline(x=1.47, ymin=0, ymax=1, label=r'CO$_2$', color='green', linestyle='dashed')
    #axs[2].vlines(x=[0.96], ymin=7e-3, ymax=1e-2, label='O3', color='blue') #O3
    axs[2].axvline(x=0.96, ymin=0, ymax=1, label=r'O$_3$', color='blue', linestyle='dashed')
    #axs[1].text(0.96, 2.996e-6, 'O3')
    #axs[2].vlines(x=[0.75], ymin=7e-3, ymax=1e-2, label='CH4', color='purple') #CH4
    axs[2].axvline(x=0.75, ymin=0, ymax=1, label=r'CH$_4$', color='purple', linestyle='dashed')
    axs[2].axvline(x=0.64, ymin=0, ymax=1, label=r'O$_2$', color='black', linestyle='dashed')

    axs[2].scatter(wave_new, djs_dens_mars, s=10)
    axs[2].set_ylim([1e-14, 0.3])
    axs[2].set_yscale('log')
    axs[2].set_xlabel(r'Wavelength ($\mu$ m)')
    axs[2].set_ylabel(r'$\mathcal{D}_{JS}$(Earth|Mars)')
    axs[2].set_title('Information difference between Earth and Mars spectra')


    
    #for ax in axs:
    #    ax.set_xscale('log')
    #    ax.set_yscale('log')
    #    #ax.set_xlim([0.3,2.1])
    #    #ax.set_ylim([1e-13, 1e-2])
    #    ax.set_xlabel('wavenumber (1/um)')
    #    ax.set_ylabel('Djs density')
        
        
    #axs[2].set_title('Djs density Earth vs. Mars')

    #plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.tight_layout()
    #plt.savefig('Djs_density_vs_lambda.png', format='png')
    plt.show()

#print JSD of interpolation function for Earth vs. interpolation function for Mars
print('Djs of Earth vs. Mars is ', djs(int_earth(wave_new),int_mars(wave_new)))
print('Djs of Earth vs. Jupiter is ', djs(int_earth(wave_new),int_jupiter(wave_new)))


print('Djs density at H2O is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.6)[1])
#print('Djs density at first CO2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.42)[1])
print('Djs density at CO2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 1.47)[1])
print('Djs density at O3 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.96)[1])
print('Djs density at CH4 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.75)[1])
print('Djs density at O2 is ', find_nearest_Djs(wave_new, djs_dens_mars, 0.64)[1])


Djs_min =  wave_new[(djs_dens_mars).argmin()]

print('Wavelength at Djs minimum is ',Djs_min)
print('The wavelengths near there are ', earth_wave[(np.abs(earth_wave-Djs_min)).argmin()], ' for Earth and ', mars_wave[(np.abs(mars_wave-Djs_min)).argmin()], ' for Mars.')









