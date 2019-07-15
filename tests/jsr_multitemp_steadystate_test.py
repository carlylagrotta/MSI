import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.jsr_steadystate as jsr
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

#test_p = pr.Processor('C:\\Users\\HP USER\\Google Drive\\Burke #Group\\Codes\\Mechanisms\\CH4_DME\\chem.cti')
test_p=pr.Processor('C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\Mechanisms\\grimech3\\gri30.cti')
jsr1 = jsr.JSR_multiTemp_steadystate(volume=8.5e-5,pressure=1.0125239,
                         temperatures=np.linspace(1050,1200,2),
                         observables=['OH','H2O'],
                         kineticSens=1,
						 physicalSens=0,
                         conditions={'CH4':0.0294,'O2':0.0296,'N2':0.941},
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_physSensHistories=0,save_timeHistories=1)
						 
solution,ksens=jsr1.run()
methane_profile=[]
#for i in jsr1.JSR_objects:
#	print(i.pressure,i.temperature,i.conditions)
#	print(i.solution['ch4'],i.reactorPressure)
#for i in range(len(jsr1.JSR_objects)):
#	methane_profile.append(jsr1.JSR_objects[i].solution['ch4'])

	
#plt.plot(np.linspace(858,1258,25),methane_profile)
#plt.savefig('C:\\Users\\HP USER\\Google Drive\\Burke Group\\Mark\\MSI\\data\\jsr_test\\methane.pdf',
#				dpi=1200, bbox_inches='tight')
print(solution)
print(ksens)
measT=[800,820,840,860,880,900]
measT=np.array(measT)+273.15
plt.figure()
plt.plot(solution['temperature'],solution['CH4'],'r')
measured=[0.032636,0.031043,0.028097,0.01566,0.01716,0.009131]
plt.plot(measT,measured,'ko')
plt.savefig('C:\\Users\\HP USER\\Desktop\\methane.pdf',dpi=1200,bbox_inches='tight')
jsr1.sensitivity_adjustment(temp_del=0.01)
jsr1.sensitivity_adjustment(pres_del=0.01)
jsr1.species_adjustment(spec_del=0.01)
print(jsr1.timeHistories)
print(len(jsr1.timeHistories))				 
#jsr1.set_geometry(volume=9.19523225755e-5)
#a,b=jsr1.run()
#print(a)
#test_tube.printVars()
#time_History = test_tube.timeHistory


#plt.plot(time_History['time']*1e3,time_History['OH']*1e6)
#time_OH = time_History['time']
#OH_ppm = time_History['OH']*1e6
#plt.figure()
#plt.plot(time_History['time']*1e3,time_History['H2O']*1e6)
#time_H2O = time_History['time']
#H2O_ppm = time_History['H2O']*1e6
