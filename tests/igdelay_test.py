import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built
import os
os.chdir('C:\\Users\\Skoron\\Desktop')
import MSI.simulations.instruments.flames as f
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import MSI.simulations.instruments.ignition_delay as ig

#test_p = pr.Processor('C:\\Users\\HP USER\\Google Drive\\Burke #Group\\Codes\\Mechanisms\\heptane_lowT\\Mech.cti')


test_p=pr.Processor('C:\\Users\\Skoron\\Desktop\\MSI\\data\\igdelay_test_H2_O2\\chem_mixed_rxns.cti')
s = ig.ignition_delay(pressure=1.909,
                         temperature=875.0,
                         observables=['OH','H2O'],
                         kineticSens=1,
						 physicalSens=0,
                         conditions={'H2':0.09506,'O2':0.19011,'N2':0.7405413023911442},
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         finalTime=1.0,
                         target='temperature',
                         target_type='max derivative',
                         save_physSensHistories=0,save_timeHistories=1,n_processors=2)		

				 
solution,sens=s.run_single()
print(solution,sens)
#methane_profile=[]
#for i in jsr1.JSR_objects:
#	print(i.pressure,i.temperature,i.conditions)
#	print(i.solution['ch4'],i.reactorPressure)
#for i in range(len(jsr1.JSR_objects)):
#	methane_profile.append(jsr1.JSR_objects[i].solution['ch4'])

	
#plt.plot(np.linspace(858,1258,25),methane_profile)
#plt.savefig('C:\\Users\\HP USER\\Google Drive\\Burke Group\\Mark\\MSI\\data\\jsr_test\\methane.pdf',
#				dpi=1200, bbox_inches='tight')
#print(solution['ch4'])
#print(ksens)



