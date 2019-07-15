import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built
import numpy as np
import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import pandas as pd
import matplotlib.pyplot as plt
test_p = pr.Processor('MSI/data/test_data/FFCM1_custom_updated.cti')
#test_tube = st.shockTube(pressure=1.672,
#                         temperature=1182,
#                         observables=['OH','H2O'],
#                         kineticSens=1,
#                         physicalSens=0,
#                         conditions={'H2O':0.001113,
#                                     'H2O2':0.002046,
#                                     'O2':0.000556,
#                                     'Ar':0.996285},
#                         initialTime=0,
#                         finalTime=0.003,
#                         thermalBoundary='Adiabatic',
#                         mechanicalBoundary='constant pressure',
#                         processor=test_p,
#                         save_timeHistories=1,
#                         save_physSensHistories=1)
#test_tube = st.shockTube(pressure=1.635,
#                         temperature=1283,
#                         observables=['OH','H2O'],
#                         kineticSens=1,
#                         physicalSens=0,
#                         conditions={'H2O':0.001113,
#                                     'H2O2':0.003094 ,
#                                     'O2':0.000556,
#                                     'Ar':0.996285},
#                         initialTime=0,
#                         finalTime=0.003,
#                         thermalBoundary='Adiabatic',
#                         mechanicalBoundary='constant pressure',
#                         processor=test_p,
#                         save_timeHistories=1,
#                         save_physSensHistories=1)

test_p = pr.Processor('MSI/data/test_data/FFCM1.cti')

test_tube = st.shockTube(pressure=0.986923,
                         temperature=295,
                         observables=['CH3','OH','HO2'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'CH3': 1.2094157562676408e-06,
                                     'HO2': 7.614839946870331e-07,
                                     'OH': 3.041863871824672e-06,
                                     'H2O2': 0.0001531112203220986,
                                     'CH4': 0.0007737003154574132,
                                     'H2O': 0.9990681757005978},
                         initialTime=0,
                         finalTime=0.003,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

test_tube.run()
k_sens = test_tube.kineticSensitivities
reactions = test_tube.processor.solution.reaction_equations()
def sort_top_uncertainty_weighted_sens(observables,S_matrix,top_sensitivity=10):
    sensitivitys =[[] for x in range(len(observables))]
    topSensitivities = [[] for x in range(len(observables))]  
    for i,observable in enumerate(observables):
        
        temp = S_matrix[:,:,i] 
        sort_s= pd.DataFrame(temp).reindex(pd.DataFrame(temp).abs().max().sort_values(ascending=False).index, axis=1)
        cc=pd.DataFrame(sort_s).iloc[:,:top_sensitivity]
        top_five_reactions=cc.columns.values.tolist()
        topSensitivities[i].append(top_five_reactions)
        ccn=pd.DataFrame(cc).as_matrix()
        sensitivitys[i].append(ccn) 
    return sensitivitys,topSensitivities
    
    

sensitivitys,topSensitivities = sort_top_uncertainty_weighted_sens(test_tube.observables,k_sens)

def plotting_senstivities(sensitivitys,
                          topSensitivities,
                          reactions,observables,
                          time,
                          working_directory):
    for i,observable in enumerate(observables):
        plt.figure()
        for column in range(np.shape(sensitivitys[i][0])[1]):
            plt.plot(time,sensitivitys[i][0][:,column],label=reactions[topSensitivities[i][0][column]])
            plt.xlabel('time (s)')
            plt.ylabel(observable)
            #plt.legend(loc='best',ncol=2)
            plt.legend(ncol=5, loc='upper left',bbox_to_anchor=(-.5,-.15))
            plt.savefig(working_directory+'/'+observable+'_ksens'+'.pdf', bbox_inches='tight')
            plt.title(str(test_tube.temperature)+'_'+str(test_tube.conditions))
plotting_senstivities(sensitivitys,topSensitivities,reactions,
                      test_tube.observables,
                      test_tube.timeHistories[0]['time'],
                      '/home/carly/Dropbox/Columbia')