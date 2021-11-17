#just a quick plotting routine

import numpy as np
import matplotlib
import matplotlib.pyplot as plt



data_jup = {}
with open('Djs-density-hotJup-vs-Jupiter-v2.txt') as f:
	for line in f:
		print(line)
		key, H2O, CO2, O3, CH4, O2, CO2High, O2High = line.split(',')
		if not key=='Planet':
			data_jup[key] = [H2O, CO2, O3, CH4, O2, CO2High, O2High]
			#[float(H2O), float(CO2), float(O3), float(CH4), float(O2), float(CO2High), float(O2High)]

wavs = {'H2O': 0.596, 'CO2':0.43, 'O3': 0.955, 'CH4': 0.752, 'O2': 0.643, 'CO2high': 1.471, 'O2high':4.5}

wavs_list = list(wavs.values())