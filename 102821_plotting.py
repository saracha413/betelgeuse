#just a quick plotting routine

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def order_by_wavelength(wave, Djs_vals):
	wave_new = wave
	wave_new.sort()
	Djs_vals_new = np.zeros(len(wave))
	for wavelength in wave:
		original_idx = wave.index(wavelength)
		new_idx = wave_new.index(wavelength)
		Djs_vals_new[new_idx] = Djs_vals[original_idx]

	return wave_new, Djs_vals_new 


#name, H2O, CO2, O3, CH4, O2, CO2high, O2high = np.loadtxt('djs_hot_jup.csv',skiprows=1, unpack=True)

data = {}
with open('Djs-hotJup-vs-Earth.txt') as f:
	for line in f:
	   key, H2O, CO2, O3, CH4, O2, CO2High, O2High = line.split(',')
	   if not key=='Planet':
	      data[key] = [float(H2O), float(CO2), float(O3), float(CH4), float(O2), float(CO2High), float(O2High)]
data_jup = {}
with open('Djs-hotJup-vs-Jupiter_v2.txt') as f:
	for line in f:
	   key, H2O, CO2, O3, CH4, O2, CO2High, O2High = line.split(',')
	   if not key=='Planet':
	      data_jup[key] = [float(H2O), float(CO2), float(O3), float(CH4), float(O2), float(CO2High), float(O2High)]

wavs = {'H2O': 0.596, 'CO2':0.43, 'O3': 0.955, 'CH4': 0.752, 'O2': 0.643, 'CO2high': 1.471, 'O2high':4.5}

wavs_list = list(wavs.values())


#HAT-P-1b and Sim1 have 5x graphite, 
#Sim4 has solar metallicity bu 0.8 C-to-O ratio; Sim3 ahs same but 0.2 C-toO
X = {'WASP-39b':1, 'WASP-6b':1, 'HD 189733b':1, 'Sim2':5,'HAT-P-1b':5,'Sim4':0.8,'Sim3':0.8,'Sim1':1,'HAT-P-12b':1,'HD 209458b':0.1}

colors = {1:'blue', 5:'gray', 0.8:'orange', 0.2:'red', 0.1:'lightblue'}

for key in data:
	if key!='Jupiter':
		plt.scatter(wavs_list, data[key], color=colors[X[key]]) #color = 'orange')
		plt.scatter(wavs_list, data_jup[key], marker='*', color=colors[X[key]]) #color='blue')
		#test_wavs, test_Djs = order_by_wavelength(wavs_list, data[key])

		#plt.plot(test_wavs, test_Djs, label=key)


for key in wavs:
	plt.annotate(key, (wavs[key],2e-5))

plt.legend()
plt.xlabel('wavelength (um)')
plt.ylabel('Djs')
plt.xscale('log')
plt.yscale('log')




plt.show()