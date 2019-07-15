import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt

test_p = pr.Processor('MSI/data/H_O2/one_reaction.cti')
test_tube = st.shockTube(pressure=0.35757249605664937,
                         temperature=1387,
                         observables=['OH'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'O2':0.001238 ,
                                     'H2O':0.001263,
                                     'H':0.0000010570824524312896,
                                     'OH':0.0000010570824524312896 ,
                                     'Ar': 0.9974968858350951},
                         initialTime=0,
                         finalTime=0.02,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

#test_tube = st.shockTube(pressure=0.35757249605664937,
#                         temperature=1387,
#                         observables=['OH'],
#                         kineticSens=0,
#                         physicalSens=0,
#                         conditions={'O2':0.0025 ,
#                                     'NH3':.0001501,
#                                     'H': 0.0000004347826086956522 ,
#                                     'NH2':0.0000004347826086956522 ,
#                                     'Ar': 0.9973490304347826},
#                         initialTime=0,
#                         finalTime=0.002,
#                         thermalBoundary='Adiabatic',
#                         mechanicalBoundary='constant volume',
#                         processor=test_p,
#                         save_timeHistories=1,
#                         save_physSensHistories=1)
test_tube.run()
#test_tube.printVars()
time_History = test_tube.timeHistory
test = time_History['H'][0]/time_History['H']

#plt.plot(time_History['time']*1e3,time_History['co']*1e6)
#time_OH = time_History['time']
#OH_ppm = time_History['co']*1e6
#plt.figure()
#plt.plot(time_History['time']*1e3,time_History['po[ome]3']*1e6)
#time_H2O = time_History['time']
#H2O_ppm = time_History['po[ome]3']*1e6
#
#plt.figure()
#plt.plot(time_History['time'],time_History['pressure']/101325)
#
#plt.figure()

