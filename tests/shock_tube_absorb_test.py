import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.simulations.absorbance.curve_superimpose as csp  
import MSI.simulations.yaml_parser as yp
import cantera as ct
import pandas as pd 
import matplotlib.pyplot as plt

test_p = pr.Processor('MSI/data/test_data/Hong.cti')
test_tube = st.shockTube(pressure=1.635,
                         temperature=1283,
                         observables=['OH','H2O'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'H2O2':0.003094,'H2O':0.001113,
                                     'O2':0.000556,'Ar':0.999623},
                         initialTime=0,
                         finalTime=0.001,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

#test_p = pr.Processor('MSI/data/test_data/optimized_burke.cti')
#test_tube = st.shockTube(pressure=3.44187,
#                         temperature=1079,
#                         observables=['OH','H2O'],
#                         kineticSens=1,
#                         physicalSens=0,
#                         conditions={'H2O2':0.00195373,'Ar':0.99804627},
#                         initialTime=0,
#                         finalTime=0.0014,
#                         thermalBoundary='Adiabatic',
#                         mechanicalBoundary='constant pressure',
#                         processor=test_p,
#                         save_timeHistories=1,
#                         save_physSensHistories=1)



OH = pd.read_csv('MSI/data/test_data/hong_oh_4.csv')
OH_time = OH['Time']
OH_ppm = OH['OH_ppm']



H2O = pd.read_csv('MSI/data/test_data/hong_h2o_4.csv')
H2O_time = H2O['Time']
H2O_ppm = H2O['H2O_ppm']


abs_csv = pd.read_csv('MSI/data/test_data/hong_abs_4.csv')
abs_time = abs_csv['time']
abs_values = abs_csv['Absorbance_227']

test_tube.run()
time_history = test_tube.timeHistory
parser = yp.Parser()
abs_instance = csp.Absorb()
#exp_loaded = parser.load_to_obj('MSI/data/test_data/Hong_4.yaml')
abs_loaded = parser.load_to_obj('MSI/data/test_data/Hong_6_high_temp_abs.yaml')
abs_data = abs_instance.superimpose_shock_tube(test_tube,abs_loaded,15.2,kinetic_sens=0)
abs_data=abs_data





plt.plot(test_tube.timeHistories[0]['time']*1000,abs_data[227])
time = test_tube.timeHistories[0]['time']
abs_data = abs_data[227]
plt.plot(abs_time*1000,abs_values,color='r')
#plt.scatter(np.array(abs_time)*1000,absorb,color='r')
#plt.axis([.01,1.4,0,.15])
#plt.plot(test_tube.timeHistories[0]['time'],test_tube.timeHistories[0]['time'])
#loaded_tube = parser.parse_shock_tube_obj(loaded_exp=exp, loaded_absorption=absp)
#uneeded for just testing absorbance



plt.figure()
plt.plot(time_history['time'],time_history['OH']*1e6)
plt.plot(OH_time,OH_ppm,color='r')

plt.figure()
plt.plot(time_history['time'],time_history['H2O']*1e6)
plt.plot(H2O_time,H2O_ppm,color='r')



