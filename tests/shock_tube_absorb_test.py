
import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')
 #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.simulations.absorbance.curve_superimpose as csp  
import MSI.simulations.yaml_parser as yp
import cantera as ct
import pandas as pd 
import matplotlib.pyplot as plt

#test_p = pr.Processor('MSI/data/hong_H2O2_fake_data/Hong_new_half_temperature_dep_new_approach.cti')
test_p = pr.Processor('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/cl2_mechanism_extra_reactions.cti')


test_tube = st.shockTube(pressure=.19,
                         temperature=240,
                         observables=['OH','H2O'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'Cl(2)':0.0027300000000000002 ,
                                     'CH4O(1)':0.000496976,
                                     'O2(7)': 0.045556129,
                                     'Cl2(5)':0.000745033,
                                     'Ar':0.95},
                         initialTime=0,
                         finalTime=0.02,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)


# OH = pd.read_csv('MSI/data/test_data/hong_oh_4.csv')
# OH_time = OH['Time']
# OH_ppm = OH['OH_ppm']



# H2O = pd.read_csv('MSI/data/test_data/hong_h2o_4.csv')
# H2O_time = H2O['Time']
# H2O_ppm = H2O['H2O_ppm']


abs_csv = pd.read_csv('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/kircher_abs_0.csv')
abs_time = abs_csv['time']
abs_values = abs_csv['Absorbance_227']

test_tube.run()
#test_tube2.run()
#test_tube3.run()

time_history = test_tube.timeHistory
# #time_history2 = test_tube2.timeHistory
# #time_history3 = test_tube3.timeHistory

parser = yp.Parser()
# #parser2 = yp.Parser()
# #parser3 = yp.Parser()

abs_instance = csp.Absorb()
# #abs_instance2 = csp.Absorb()
# #abs_instance3 = csp.Absorb()


# #exp_loaded = parser.load_to_obj('MSI/data/test_data/Hong_4.yaml')
abs_loaded = parser.load_to_obj('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/Kircher_0_abs.yaml')
abs_data = abs_instance.superimpose_shock_tube(test_tube,abs_loaded,670,kinetic_sens=0)

# #abs_loaded2 = parser2.load_to_obj('MSI/data/hong_H2O2_fake_data/Hong_HO2_fake_data_3_abs.yaml')
# #abs_data2 = abs_instance2.superimpose_shock_tube(test_tube2,abs_loaded2,15.2,kinetic_sens=0)
# #
# #abs_loaded3 = parser.load_to_obj('MSI/data/hong_H2O2_fake_data/Hong_HO2_fake_data_3_abs.yaml')
# #abs_data3 = abs_instance.superimpose_shock_tube(test_tube3,abs_loaded3,15.2,kinetic_sens=0)


plt.xlabel('Time (ms)')
plt.ylabel('Absorbance 227nm')
plt.plot(test_tube.timeHistories[0]['time']*1e3,abs_data[227],label='constant value rate constants')
# #plt.plot(test_tube2.timeHistories[0]['time']*1e3,abs_data2[227],label='k1,k2=constant value  k3,k4=expression from paper')
# #plt.plot(test_tube3.timeHistories[0]['time']*1e3,abs_data3[227],label='k1,k2=derived A value  k3,k4=expression from paper')

# time = test_tube.timeHistories[0]['time'].values
# abs_data = abs_data[227]
plt.plot(abs_time*1e3,abs_values,color='r',label = 'Experimental Data')
# plt.legend(loc='lower left')

# OH_data = time_history['OH'].values*1e6
# H2O_data = time_history['H2O'].values*1e6



# plt.figure()
# plt.xlabel('Time (ms)')
# plt.ylabel('OH ppm')
# plt.plot(time_history['time']*1e3,time_history['OH']*1e6,label='constant value rate constants')
# #plt.plot(time_history2['time']*1e3,time_history2['OH']*1e6,label='k1,k2=constant value  k3,k4=expression from paper')
# #plt.plot(time_history3['time']*1e3,time_history3['OH']*1e6,label='k1,k2=derived A value  k3,k4=expression from paper')

# plt.plot(OH_time*1e3,OH_ppm,color='r',label = 'Experimental Data')
# plt.legend()
# plt.figure()
# plt.xlabel('Time (ms)')
# plt.ylabel(r'H$_2$O ppm')
# plt.plot(time_history['time']*1e3,time_history['H2O']*1e6,label='constant value rate constants')
# #plt.plot(time_history2['time']*1e3,time_history2['H2O']*1e6,label='k1,k2=constant value  k3,k4=expression from paper')
# #plt.plot(time_history3['time']*1e3,time_history3['H2O']*1e6,label='k1,k2=derived A value  k3,k4=expression from paper')


# plt.plot(H2O_time*1e3,H2O_ppm,color='r',label = 'Experimental Data')
# plt.legend()



