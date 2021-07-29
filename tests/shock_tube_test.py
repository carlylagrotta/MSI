import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

test_p = pr.Processor('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/branching_reaction_study/ffcm1_branching_reactions_CH3OH.cti')
test_p2 = pr.Processor('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/branching_reaction_study/ffcm1_branching_reactions_theory.cti')


test_tube = st.shockTube(pressure=0.31232922337880797414 ,
                         temperature=1226, 
                         observables=['OH'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'C4H10O2': 0.000013636,
                                     'Ar':0.999986364},
                         initialTime=0,
                         finalTime=0.0003,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1
                         )

test_tube.run()
#test_tube.printVars()
time_History = test_tube.timeHistory


test_tube2 = st.shockTube(pressure=0.31232922337880797414 ,
                         temperature=1226, 
                         observables=['OH'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'C4H10O2': 0.000013636,
                                     'Ar':0.999986364},
                         initialTime=0,
                         finalTime=0.0003,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p2,
                         save_timeHistories=1,
                         save_physSensHistories=1
                         )
test_tube2.run()
time_History2 = test_tube2.timeHistory


plt.figure() 
df_experiment = pd.read_csv('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/branching_reaction_study/Hong_H2O_raw_data_fig_1_2_OH_O.csv')
concentration = np.true_divide(1,time_History['temperature'].to_numpy())*time_History['pressure'].to_numpy()
                           
concentration *= (1/(8.314e6))*time_History['H2O'].dropna().to_numpy()
plt.plot(time_History['time'],time_History['OH']*1e6,label='no theory')
plt.plot(time_History2['time'],time_History2['OH']*1e6,label=' theory')
#plt.plot(df_experiment['Time'],df_experiment['H2O_ppm'])

plt.legend( )

plt.savefig('/Users/carlylagrotta/Desktop/OH_new.pdf')