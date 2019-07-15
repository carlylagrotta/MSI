
import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.optimization.opt_runner as opt
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import MSI.optimization.shock_tube_optimization_shell as stMSI
import MSI.master_equation.master_equation as meq
import cantera as ct
import numpy as np 



                 
MSI_st_instance = stMSI.MSI_shocktube_optimization('FFCM1.cti',.01,1,1,'MSI/data/test_data',
                                 [['Hong_4.yaml','Hong_4_abs.yaml'],['Hong_1.yaml'],
                  ['Troe_6.yaml','Troe_6_abs.yaml']],                 
                   'uncertainty_test.csv','FFCM1_target_values.csv' )
MSI_st_instance.one_run_shock_tube_optimization()

X = MSI_st_instance.X
experiment_dictonaries = MSI_st_instance.experiment_dictonaries
list_of_parsed_yamls = MSI_st_instance.list_of_parsed_yamls

#just use this array as a test 
a1 = np.array([-4.71844785e-14,  -1.11022302e-14,   4.16333634e-15,   5.55111512e-15,
  5.55111512e-15,  -1.64246079e-02,   5.42951118e-03,   1.81168588e-03,
 -2.88407459e-03,   4.68262347e-03,  -8.32667268e-15,   2.22044605e-14,
 -1.11022302e-14,  -0.00000000e+00,  -1.11022302e-14,  -1.64246079e-02,
  5.42951118e-03,   1.81168588e-03,  -2.88407459e-03,   4.68262347e-03,
 -1.94289029e-14,  -2.49800181e-14,   2.77555756e-15,   2.77555756e-15,
 -2.77555756e-15,  -1.64246079e-02,   5.42951118e-03,   1.81168588e-03,
 -2.88407459e-03,   4.68262347e-03,  -4.16333634e-14,  -2.49800181e-14,
  2.77555756e-15,   2.77555756e-15,  -2.77555756e-15,  -1.64246079e-02,
  5.42951118e-03,   1.81168588e-03,  -2.88407459e-03,   4.68262347e-03,
 -4.16333634e-14,  -2.49800181e-14,   2.77555756e-15,   2.77555756e-15,
 -2.77555756e-15,  -1.64246079e-02,   5.42951118e-03,   1.81168588e-03,
 -2.88407459e-03,   4.68262347e-03,  -4.16333634e-14,  -2.49800181e-14,
  2.77555756e-15,   2.77555756e-15,  -2.77555756e-15,  -1.64246079e-02,
  5.42951118e-03,   1.81168588e-03,  -2.88407459e-03,   4.68262347e-03,
 -4.16333634e-14,  -2.49800181e-14,   2.77555756e-15,   2.77555756e-15,
 -2.77555756e-15,  -1.64246079e-02,   5.42951118e-03,   1.81168588e-03,
 -2.88407459e-03,   4.68262347e-03,  -4.16333634e-14,  -2.49800181e-14,
  2.77555756e-15,   2.77555756e-15,  -2.77555756e-15,  -1.64246079e-02,
  5.42951118e-03,   1.81168588e-03,  -2.88407459e-03,   4.68262347e-03,
 -4.16333634e-14,  -2.49800181e-14,   2.77555756e-15,   2.77555756e-15,
 -2.77555756e-15,  -1.64246079e-02,   5.42951118e-03,   1.81168588e-03,
 -2.88407459e-03,   4.68262347e-03,  -4.16333634e-14,  -2.49800181e-14,
  2.77555756e-15,   2.77555756e-15,  -2.77555756e-15,  -1.64246079e-02,
  5.42951118e-03,   1.81168588e-03,  -2.88407459e-03,   4.68262347e-03])
a1 = a1.reshape((25,4))
  
sensitivity_dict = {'H + HCO <=> CO + H2':[a1,a1,a1]}


master_equation_instance = meq.Master_Equation()
mapped_to_alpha_full_simulation,nested_list = master_equation_instance.map_to_alpha(sensitivity_dict,
                                      experiment_dictonaries,
                                      list_of_parsed_yamls,
                                      ['H + HCO <=> CO + H2'])

    
S_matrix = master_equation_instance.map_to_S(nested_list,sensitivity_dict, ['H + HCO <=> CO + H2'])

matrix_instance = ml.Adding_Target_Values(MSI_st_instance.S_matrix,MSI_st_instance.Y_matrix,MSI_st_instance.sigma)
matrix_instance.target_values_for_S(MSI_st_instance.data_directory+'/'+MSI_st_instance.k_target_values_csv,MSI_st_instance.experiment_dictonaries,
                                   ['H + HCO <=> CO + H2'],sensitivity_dict)
