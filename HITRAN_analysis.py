import numpy as np
import matplotlib
import matplotlib.pyplot as plt


nu_spec, S_spec, ID_spec, _, _, _, _, _ = np.loadtxt('/Users/saravannah/AstroSpec/10-21-astropec_HITRAN/61718276.out.txt', unpack = True)



nu_H2O, S_H2O = [], []
nu_CO2, S_CO2 = [], []
nu_O3, S_O3 = [], [] #ozone
nu_CO, S_CO = [], []
nu_CH4, S_CH4 = [], [] #methane
nu_O2, S_O2 = [], []

for i in range(len(nu_spec)):
	if ID_spec[i] == 1:
		nu_H2O.append(nu_spec[i])
		S_H2O.append(S_spec[i])
	elif ID_spec[i] == 2:
		nu_CO2.append(nu_spec[i])
		S_CO2.append(S_spec[i])
	elif ID_spec[i] == 3:
		nu_O3.append(nu_spec[i])
		S_O3.append(S_spec[i])
	elif ID_spec[i] == 5:
		nu_CO.append(nu_spec[i])
		S_CO.append(S_spec[i])
	elif ID_spec[i] == 6:
		nu_CH4.append(nu_spec[i])
		S_CH4.append(S_spec[i])
	elif ID_spec[i] == 7:
		nu_O2.append(nu_spec[i])
		S_O2.append(S_spec[i])
	else:
		print('ID ', ID_spec[i], ' does not match metatadata.')




plt.scatter(1e3*np.reciprocal(nu_H2O), S_H2O, label='H2O', s=0.5)
plt.scatter(1e3*np.reciprocal(nu_CO2), S_CO2, label='CO2', s=0.5)
plt.scatter(1e3*np.reciprocal(nu_O3), S_O3, label='O3', s=0.5)
plt.scatter(1e3*np.reciprocal(nu_CO), S_CO, label='CO', s=0.5)
plt.scatter(1e3*np.reciprocal(nu_CH4), S_CH4, label='CH4', s=0.5)
plt.scatter(1e3*np.reciprocal(nu_O2), S_O2, label='O2', s=0.5)
plt.xlabel('Wavelength (um)')
plt.ylabel('Line intensity')

#plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()