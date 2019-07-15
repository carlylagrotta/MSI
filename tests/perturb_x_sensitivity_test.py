import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.optimization.opt_runner as opt
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import MSI.optimization.shock_tube_optimization_shell as stMSI
import cantera as ct
import MSI.utilities.plotting_script as plotter
import MSI.utilities.perturb_X_for_testing_shock_tube_shell as perturbX
import matplotlib.pyplot as plt

import numpy as np
 #burke_target_value_test.csv                
files_to_include = [['Hong_1.yaml']]                                  
numer_of_iterations = 3
cti_file = 'chem_original_burke.cti'
working_directory = 'MSI/data/test_data'
reaction_uncertainty_csv = 'burke_uncertainty_test.csv'

rate_constant_target_value_data = ''
#this would be an empty string '' if you do not want to include it 

#this could be 'On'




MSI_st_instance_one = stMSI.MSI_shocktube_optimization(cti_file,
                                                   .01,
                                                   1,
                                                   1,
                                                   working_directory,
                                                   files_to_include,                 
                                                   reaction_uncertainty_csv,'' )
MSI_st_instance_one.one_run_shock_tube_optimization()

S_matrix_original = MSI_st_instance_one.S_matrix
exp_dict_list_original = MSI_st_instance_one.experiment_dictonaries
X_one_itteration = MSI_st_instance_one.X
experimental_dict_uncertainty_original = MSI_st_instance_one.experiment_dict_uncertainty_original
Y_matrix_original = MSI_st_instance_one.Y_matrix
OH =  exp_dict_list_original[0]['simulation'].timeHistories[0]['OH']
time =  exp_dict_list_original[0]['simulation'].timeHistories[0]['time']

Sij_list = []
Y_difference_list = []
S_new_list= []
S_percent_difference_list = []
print('_________________________perturbing______________________________')
for row_in_X in range(np.shape(X_one_itteration)[0]):
#for row_in_X in range(5):
    MSI_st_perturb_instance = perturbX.perturb_X_Shell(cti_file,
                                                          .01,
                                                           1,
                                                           1,
                                                           working_directory,
                                                           files_to_include,
                                                           reaction_uncertainty_csv,'',
                                                           shape_of_X = np.shape(X_one_itteration),
                                                           shape_of_X_counter= row_in_X,
                                                           S_matrix_original = S_matrix_original,
                                                           Y_matrix_original = Y_matrix_original,
                                                           experimental_dict_uncertainty_original = experimental_dict_uncertainty_original,
                                                           original_experimental_dicts = exp_dict_list_original)
    MSI_st_perturb_instance.multiple_shock_tube_runs(2)
    
    experimental_dict_perturbed = MSI_st_perturb_instance.experiment_dictonaries
    plt.figure()
    plt.title(str(row_in_X))
    plt.plot(time,OH)
    plt.plot(experimental_dict_perturbed[0]['simulation'].timeHistories[0]['time'], experimental_dict_perturbed[0]['simulation'].timeHistories[0]['OH'])
    y_new = MSI_st_perturb_instance.Y_matrix
    Sij_list.append(MSI_st_perturb_instance.Sij)
    Y_difference_list.append(MSI_st_perturb_instance.Y_difference)
    S_new_list.append(MSI_st_perturb_instance.S_new)
    S_percent_difference_list.append(MSI_st_perturb_instance.S_percent_difference)
                                                       
    
Sij = sum(Sij_list)
Y_difference = sum(Y_difference_list)
S_new = sum(S_new_list)
Percent_difference = sum(S_percent_difference_list)




#run the new testing class
#make a new instnace of the testing class and run the residuals calculator 
#store all the resudals and then add the matrix together 


