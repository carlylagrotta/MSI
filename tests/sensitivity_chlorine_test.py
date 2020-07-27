import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')#get rid of this at some point with central test script or when package is built
import numpy as np
import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import pandas as pd
import matplotlib.pyplot as plt


#test_p = pr.Processor('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/chem_calories_extra_reaction.cti')
#test_p = pr.Processor('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/chem_calories_extra_reaction_editing.cti')
test_p = pr.Processor('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/cl2_mechanism_extra_reactions.cti')

#test_p = pr.Processor('MSI/data/klip_optimization_with_raw_hong_data_chlorine_mechanism/chem_calories_extra_reaction_editing_AnEa.cti')

# test_tube = st.shockTube(pressure=1,
#                          temperature=777,
#                          observables=['HO2(6)','H2O2(15)'],
#                          kineticSens=1,
#                          physicalSens=0,
#                          conditions={'Cl2(5)': 0.0037048 ,'CH4O(1)':0.000211705,
#                                      'HO2(6)':.00001852,'O2(7)':0.1961539751,'N2(8)':0.7846159003},
#                          initialTime=0,
#                          finalTime=0.01,
#                          thermalBoundary='Adiabatic',
#                          mechanicalBoundary='constant volume',
#                          processor=test_p,
#                          save_timeHistories=1,
#                          save_physSensHistories=1)

test_tube = st.shockTube(pressure=.19,
                          temperature=240,
                          observables=['HO2(6)','H2O2(15)'],
                          kineticSens=1,
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
# test_tube = st.shockTube(pressure=.92,
#                          temperature=298,
#                          observables=['HO2(6)','H2O2(15)'],
#                          kineticSens=1,
#                          physicalSens=0,
#                          conditions={'OCCO(401)':.00525 ,
#                                      'O2(7)': .129,
#                                      'Ar':0.886},
#                          initialTime=0,
#                          finalTime=0.02,
#                          thermalBoundary='Adiabatic',
#                          mechanicalBoundary='constant pressure',
#                          processor=test_p,
#                          save_timeHistories=1,
#                          save_physSensHistories=1)
test_tube.run()
k_sens_new = test_tube.kineticSensitivities
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
        ccn=pd.DataFrame(cc).values
        sensitivitys[i].append(ccn) 
    return sensitivitys,topSensitivities
    
    

sensitivitys,topSensitivities = sort_top_uncertainty_weighted_sens(test_tube.observables,k_sens_new,top_sensitivity=10)

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
            plt.legend(ncol=2, loc='upper left',bbox_to_anchor=(-.5,-.15))
            #plt.savefig(working_directory+'/'+observable+'_ksens'+'.pdf', bbox_inches='tight')
            plt.title(str(test_tube.temperature)+'_'+str(test_tube.conditions))

plotting_senstivities(sensitivitys,topSensitivities,reactions,
                      test_tube.observables,
                      test_tube.timeHistories[0]['time'],'/home/carly/Dropbox/Columbia')

