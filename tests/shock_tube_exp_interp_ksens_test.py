import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built
import numpy as np
import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct

test_p = pr.Processor('MSI/data/test_data/glarborg_custom.cti')
test_tube = st.shockTube(pressure=0.986923,
                         temperature=315,
                         observables=['CH3','OH'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'CH3':0.0000012094157562676408,'HO2':0.0000007614839946870331,
                                     'OH':0.000003041863871824672,'H2O2':0.0001531112203220986,'N2O':0.0007737003154574132},
                         initialTime=0,
                         finalTime=0.003,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

csv_paths = ['MSI/data/test_data/sangwan_ch3_0.csv','MSI/data/test_data/sangwan_oh_0.csv']
exp_data = test_tube.importExperimentalData(csv_paths)
test_tube.run()
data = test_tube.interpolate_experimental_kinetic()
int_ksens_exp_mapped= test_tube.map_and_interp_ksens()
#print(data)


