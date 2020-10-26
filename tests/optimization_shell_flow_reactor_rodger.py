import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')

#os.chdir('../../') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.optimization.opt_runner as opt
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import MSI.optimization.optimization_shell as stMSI
import cantera as ct
import MSI.utilities.plotting_script as plotter


files_to_include = [['YAML_ignition_delay_Mathieu_Tab6.yaml']]
                    

                    # ['YAML_batch_reactor_Mueller_Fig2.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig4a.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig5a.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig5b.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig5c.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_1.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_3.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_5.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_65.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_8.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_10.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig12_025equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig12_05equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig12_1equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13a_025equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13a_05equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13a_1equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13b_025equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13b_05equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13b_1equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_12ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_24ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_52ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_205ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_313ppmNO.yaml']]
                    # ['YAML_batch_reactor_Mueller_Fig15.yaml']]
    
                    # ['YAML_flow_reactor_Glar_Fig34_N2O_O2_01.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_N2O_O2_04.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_N2O_O2_10.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_N2O_O2_20.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_N2O_O2_50.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_NO2_O2_20.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_NO2_O2_50.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_NO2_O2_01.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_NO2_O2_04.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig34_NO2_O2_10.yaml'],
    
                    # ['YAML_batch_reactor_Yetter_Fig8.yaml'],
                    # ['YAML_batch_reactor_Yetter_Fig10_1atm.yaml'],
                    # ['YAML_batch_reactor_Yetter_Fig10_3atm.yaml'],
                    # ['YAML_batch_reactor_Yetter_Fig10_6atm.yaml'],
                    # ['YAML_batch_reactor_Yetter_Fig10_9atm.yaml'],
                    
                    # ['YAML_batch_reactor_Mueller_Fig2.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig4a.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig5a.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig5b.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig5c.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_1.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_3.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_5.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_65.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_8.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig8_10.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig12_025equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig12_05equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig12_1equ.yaml']]
                    # ['YAML_batch_reactor_Mueller_Fig13a_025equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13a_05equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13a_1equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13b_025equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13b_05equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig13b_1equ.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_12ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_24ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_52ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_205ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig14_313ppmNO.yaml'],
                    # ['YAML_batch_reactor_Mueller_Fig15.yaml']]

                    # ['YAML_flow_reactor_Glar_Fig28_25ppm.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig28_50ppm.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig28_100ppm.yaml'],
                    
                    # ['YAML_flow_reactor_Glar_Fig30_NO_H2O_00.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig30_NO_H2O_10.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig30_NO_H2O_28.yaml'],
                    
                    # ['YAML_flow_reactor_Glar_Fig33_NO_O2_01.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig33_NO_O2_04.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig33_NO_O2_10.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig33_NO_O2_20.yaml'],
                    # ['YAML_flow_reactor_Glar_Fig33_NO_O2_50.yaml'],
                    
                    # ['YAML_flow_reactor_Zhang_S1_1_H2.yaml'],
                    # ['YAML_flow_reactor_Zhang_S1_1_O2.yaml'],
                    # ['YAML_flow_reactor_Zhang_S1_2_H2.yaml'],
                    # ['YAML_flow_reactor_Zhang_S1_2_O2.yaml'],
                    # ['YAML_flow_reactor_Zhang_S1_3_H2.yaml'],
                    # ['YAML_flow_reactor_Zhang_S1_4_H2.yaml']]


# Currently not being used for optimization...
#['YAML_flow_reactor_Zhang_S2_NO_8_O2.yaml'],
#['YAML_flow_reactor_Zhang_S2_NO_18_O2.yaml'],
#['YAML_flow_reactor_Zhang_S2_NO_209_O2.yaml'],
#['YAML_flow_reactor_Zhang_S2_NO_045_O2_short.yaml'],
                   

#files_to_include = [['Masten_0.yaml'],
#                    ['Masten_1.yaml'],
#                    ['Masten_2.yaml'],
#                    ['Masten_3.yaml'],
#                    ['Masten_4.yaml'],
#                    ['Masten_5.yaml'],
#                    ['Masten_6.yaml'],
#                    ['Masten_7.yaml'],
#                    ['Masten_8.yaml'],
#                    ['Masten_9.yaml']]
#------------------------------------------------------------------
# Reactions of Importance
#------------------------------------------------------------------
#
#------------------------------------------------------------------
numer_of_iterations = 2           
working_directory = 'MSI\All_Optimization_Data\EnergeticMaterials'
cti_file = 'Glarborg_HNO.cti'
            # 'Glarborg_HNO.cti'
reaction_uncertainty_csv = 'Glarborg_HNO_uncertanties.csv'

#rate_constant_target_value_data = 'burke_target_value_single_reactions.csv'
#rate_constant_target_value_data = 'target_reactions_test.csv'
rate_constant_target_value_data = ''
#this would be an empty string '' if you do not want to include it 
run_with_k_target_values = 'Off'
#this could be 'On'

rate_constant_target_value_data_for_plotting = 'target_reactions_test.csv'



MSI_st_instance_one = stMSI.MSI_optimization(cti_file,
                                                   .01,
                                                   1,
                                                   1,
                                                   working_directory,
                                                   files_to_include,                 
                                                   reaction_uncertainty_csv,rate_constant_target_value_data )
MSI_st_instance_one.one_run_optimization()

S_matrix_original = MSI_st_instance_one.S_matrix
exp_dict_list_original = MSI_st_instance_one.experiment_dictonaries
original_covariance = MSI_st_instance_one.covarience
X_one_itteration = MSI_st_instance_one.X
MSI_st_instance_one.deltaXAsNsEas
z_df_original = MSI_st_instance_one.z_data_frame
    



#need to fix this and return _s_matrix and y_matrix



MSI_st_instance_two = stMSI.MSI_optimization(cti_file,
                                                   .01,
                                                   1,
                                                   1,
                                                   working_directory,
                                                   files_to_include,                 
                                                   reaction_uncertainty_csv,rate_constant_target_value_data )






X_list = MSI_st_instance_two.multiple_runs(numer_of_iterations)


deltaXAsNsEas = MSI_st_instance_two.deltaXAsNsEas
physical_obervable_updates_list = MSI_st_instance_two.physical_obervable_updates_list
absorbance_observables_updates_list = MSI_st_instance_two.absorbance_coef_update_dict
Ydf = MSI_st_instance_two.Y_data_frame
Zdf = MSI_st_instance_two.z_data_frame
experimental_dicts = MSI_st_instance_two.experiment_dictonaries
z_matrix = MSI_st_instance_two.z_matrix
s_matrix = MSI_st_instance_two.s_matrix
y = MSI_st_instance_two.y_matrix
Y_matrix = MSI_st_instance_two.Y_matrix
S_matrix = MSI_st_instance_two.S_matrix

X = MSI_st_instance_two.X
Xdf = MSI_st_instance_two.X_data_frame
covarience = MSI_st_instance_two.covarience
exp_dict_list_optimized = MSI_st_instance_two.experiment_dictonaries
parsed_yaml_list = MSI_st_instance_two.list_of_parsed_yamls
sigma = MSI_st_instance_two.sigma
X = MSI_st_instance_two.X
delta_X = MSI_st_instance_two.delta_X
#target_value_rate_constant_csv = 'MSI/data/test_data/FFCM1_custom_target_value_test.csv'
original_cti_file = MSI_st_instance_two.data_directory +'/'+ MSI_st_instance_two.cti_file_name

experiment_dict_uncertainty = MSI_st_instance_two.experiment_dict_uncertainty_original
target_value_csv = MSI_st_instance_two.data_directory +'/'+ MSI_st_instance_two.k_target_values_csv

if run_with_k_target_values == 'On' or run_with_k_target_values == 'on':
    k_target_value_S_matrix = MSI_st_instance_two.k_target_values_for_s
else:
    k_target_value_S_matrix = None


##########################################################################################################################
#PLOTTING##
##########################################################################################################################


plotting_instance = plotter.Plotting(S_matrix,
                                      s_matrix,
                                      Y_matrix,
                                      Y_matrix,
                                      z_matrix,
                                      X,
                                      sigma,
                                      covarience,
                                      original_covariance,
                                      S_matrix_original,
                                      exp_dict_list_optimized,
                                      exp_dict_list_original,
                                      parsed_yaml_list,
                                      Ydf,
                                      target_value_rate_constant_csv= '',
                                      k_target_value_S_matrix =k_target_value_S_matrix,
                                      k_target_values=run_with_k_target_values,
                                      working_directory=working_directory,
                                      shock_tube_instance = MSI_st_instance_two,
                                      optimized_cti_file=MSI_st_instance_two.new_cti_file,
                                      original_cti_file=original_cti_file)

observable_counter_and_absorbance_wl,length_of_experimental_data = plotting_instance.lengths_of_experimental_data()
sigmas_optimized,test = plotting_instance.calculating_sigmas(S_matrix,covarience)
sigmas_original,test2 = plotting_instance.calculating_sigmas(S_matrix_original,original_covariance)
plotting_instance.plotting_observables(sigmas_original = sigmas_original,sigmas_optimized= sigmas_optimized)
diag = plotting_instance.getting_matrix_diag(covarience)

#plotting_instance.Y_matrix_plotter(Y_matrix,exp_dict_list_optimized,y,sigma)



plotting_instance.plotting_rate_constants(optimized_cti_file=MSI_st_instance_two.new_cti_file,
                                original_cti_file=original_cti_file,
                                initial_temperature=250,
                                final_temperature=2500)
                                


sensitivity, top_sensitivity = plotting_instance.sort_top_uncertainty_weighted_sens()
obs = plotting_instance.plotting_uncertainty_weighted_sens()




#plotting_instance.difference_plotter(['A_1','n_1','Ea_1','A_0'],
#                                                optimized_cti_file=MSI_st_instance_two.new_cti_file,
#                                                pdf_distribution_file='MSI/data/automating_HO2_masten/graph_read_pdf.csv')


