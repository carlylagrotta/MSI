import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# test_p = pr.Processor('MSI/data/test_data/FFCM1_custom.cti')
# #test_p = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb.cti')
# #test_p = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb_test.cti')

# test_tube = st.shockTube(pressure=1.909,
#                          temperature=1398,
#                          observables=['H2O'],
#                          kineticSens=0,
#                          physicalSens=0,
#                          conditions={'O2':0.00062 ,
#                                      'H2O':0.001234,
#                                      'H2O2':0.00254 ,
#                                      'Ar': 0.9974968858350951},
#                          initialTime=0,
#                          finalTime=.07,
#                          thermalBoundary='Adiabatic',
#                          mechanicalBoundary='constant volume',
#                          processor=test_p,
#                          save_timeHistories=1,
#                          save_physSensHistories=1,
#                          volumeTrace='MSI/data/RCM/volume_trace_file.csv')

# test_tube.run()
# #test_tube.printVars()
# time_History = test_tube.timeHistory

test_p = pr.Processor('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/branching_reaction_study/FFCM1_custom_cheb_extra_zeros_new_extra_reactions_precursor.cti')

test_tube = st.shockTube(pressure=0.6008089141308347 ,
                         temperature=1925, 
                         observables=['OH'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'CH3OH': 0.000023898,
                                     'Ar': 0.999976102},
                         initialTime=0,
                         finalTime=.002,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1
                         )

test_tube.run()
#test_tube.printVars()
time_History = test_tube.timeHistory
# plt.plot(time_History['time'],time_History['pressure']/100000)
# plt.ylim(0,140)
# plt.figure()
# plt.plot(time_History['time'],time_History['temperature'])
plt.figure() 

concentration = np.true_divide(1,time_History['temperature'].to_numpy())*time_History['pressure'].to_numpy()
                           
concentration *= (1/(8.314e6))*time_History['OH'].dropna().to_numpy()
plt.plot(time_History['time'],concentration)
plt.savefig('/Users/carlylagrotta/Desktop/OH.pdf')