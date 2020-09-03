import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt

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

test_p = pr.Processor('MSI/data/RCM/mech.cti')

test_tube = st.shockTube(pressure=30.11e5/101325,
                         temperature=682, 
                         observables=['H2O'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'C3H8': 0.024950099800399205,
                                     'CH3OCH3': 0.024950099800399205,
                                     'N2': 0.7504990019960079,
                                     'O2': 0.1996007984031936},
                         initialTime=0,
                         finalTime=.1,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1
                         )

test_tube.run()
#test_tube.printVars()
time_History = test_tube.timeHistory
plt.plot(time_History['time'],time_History['pressure']/100000)
plt.ylim(0,140)
plt.figure()
plt.plot(time_History['time'],time_History['temperature'])

