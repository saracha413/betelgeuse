import numpy as np
import scipy as sp
from scipy.fft import fft
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
 
