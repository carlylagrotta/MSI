import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt

test_p = pr.Processor('MSI/data/test_data/FFCM1_custom.cti')
test_p2 = pr.Processor('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/rmg_chlorine_custom_calories.cti')
#test_p = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb.cti')
#test_p = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb_test.cti')

test_tube = st.shockTube(pressure=1.909,
                          temperature=1398,
                          observables=['H2O'],
                          kineticSens=0,
                          physicalSens=0,
                          conditions={'O2':0.00062 ,
                                      'H2O':0.001234,
                                      'H2O2':0.00254 ,
                                      'Ar': 0.9974968858350951},
                          initialTime=0,
                          finalTime=.0002,
                          thermalBoundary='Adiabatic',
                          mechanicalBoundary='constant volume',
                          processor=test_p,
                          save_timeHistories=1,
                          save_physSensHistories=1)

test_tube.run()
#test_tube.printVars()
time_History = test_tube.timeHistory



test_tube2 = st.shockTube(pressure=1.909,
                          temperature=1398,
                          observables=['H2O(10)'],
                          kineticSens=0,
                          physicalSens=0,
                          conditions={'O2(6)':0.00062 ,
                                      'H2O(10)':0.001234,
                                      'H2O2(13)':0.00254 ,
                                      'Ar': 0.9974968858350951},
                          initialTime=0,
                          finalTime=.0002,
                          thermalBoundary='Adiabatic',
                          mechanicalBoundary='constant volume',
                          processor=test_p2,
                          save_timeHistories=1,
                          save_physSensHistories=1)

test_tube2.run()
#test_tube.printVars()
time_History2 = test_tube2.timeHistory


plt.plot(time_History['time'],time_History['H2O'])
plt.plot(time_History2['time'],time_History2['H2O(10)'])
