# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:33:39 2020

@author: Skoron
"""

import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
#sys.path.append('C:\\Users\\Skoron\\Desktop')
os.chdir('../../')
#print(os.getcwd())
#import sys
import fire
import MSI.optimization.optimization_shell as shell
import MSI.utilities.plotting_script as plotter
import MSI.optimization.optimization_shell_chebyshev as MSIcheb
import re
import pandas as pd
import numpy as np
import ast

class multiscale_informatics:
    
    
    def __init__(self,input_options):
        self.input_options=input_options
        self.get_parameters()
        
        
    def run_msi(self):
        print(self.wdir)
        os.chdir(self.wdir)
        for m,mech in enumerate(self.models):
            if self.run_individual_yaml:
                #
                for y,yamlfile in enumerate(self.yaml_files):
                    if not self.master_equation_models:
                        self.msi_no_master(m,yaml_input=[self.yaml_files[y]])
                        self.plotting_no_master(m,file_identifier=self.yaml_files[y][0].strip('.yaml'))
                    elif self.master_equation_models:
                        self.msi_master(m,yaml_input=[self.yaml_files[y]])
                        self.plotting_master(m, file_identifier=self.yaml_files[y][0].strip('.yaml'))
                
            elif not self.run_individual_yaml:
                
                if not self.master_equation_models:
                    self.msi_no_master(m,yaml_input=self.yaml_files)
                    self.plotting_no_master(m,file_identifier='all_yamls')
                    
                elif self.master_equation_models:
                    self.msi_master(m,yaml_input=self.yaml_files)
                    self.plotting_master(m,file_identifier='all_yamls')
                
                
    def write_convergence(self):
        #print(type(self.X_list))
        #print(len(self.X_list))
        convergence=pd.DataFrame(columns=['Parameter'])
        for i in range(len(self.X_list)):
            #print(np.shape(self.X_list[i]))
            convergence[str(i)]=np.array(self.X_list[i]).flatten()
        print(self.Xdf['value'])
        convergence['Parameter']=self.Xdf['value']
        #print(self.Xdf.columns)
        convergence.to_csv(self.wdir+'\\convergence_data.csv',index=False)
        
    def msi_master(self,iteration,yaml_input=[]):
        current_master_uncertainties=pd.read_csv(self.master_uncertainties[iteration])
        self.get_chebyshev_coefficients(iteration)
        #print(self.cheb_coeffs)
        MSI_instance_one=MSIcheb.MSI_optimization_chebyshev(self.models[iteration],
                                                            0.01,
                                                            1,
                                                            1,
                                                            self.wdir,
                                                            yaml_input,
                                                            self.reaction_uncertainties[iteration],
                                                            k_target_values_csv=self.targets[iteration],
                                                            master_equation_reactions=self.master_equation_reactions[iteration],
                                                            chebyshev_sensitivities=self.cheb_coeffs[iteration],
                                                            master_reaction_equation_cti_name=self.master_equation_models[iteration],
                                                            master_index=self.indices[iteration],
                                                            master_equation_uncertainty_df=current_master_uncertainties,
                                                            chebyshev_fit_nominal_parameters_dict= None)
        MSI_instance_one.one_run_optimization()
        self.S_matrix_original=MSI_instance_one.S_matrix
        self.exp_dict_list_original=MSI_instance_one.experiment_dictonaries
        self.original_covariance=MSI_instance_one.covarience
        self.X_one_iteration=MSI_instance_one.X
        self.z_df_original=MSI_instance_one.z_data_frame
        self.y_df_original=MSI_instance_one.Y_data_frame
        
        
        self.MSI_instance_two=MSIcheb.MSI_optimization_chebyshev(self.models[iteration],
                                                            0.01,
                                                            1,
                                                            1,
                                                            self.wdir,
                                                            yaml_input,
                                                            self.reaction_uncertainties[iteration],
                                                            k_target_values_csv=self.targets[iteration],
                                                            master_equation_reactions=self.master_equation_reactions[iteration],
                                                            chebyshev_sensitivities=self.cheb_coeffs[iteration],
                                                            master_reaction_equation_cti_name=self.master_equation_models[iteration],
                                                            master_index=self.indices[iteration],
                                                            master_equation_uncertainty_df=current_master_uncertainties,
                                                            chebyshev_fit_nominal_parameters_dict= None)
        
        self.X_list=self.MSI_instance_two.multiple_runs(self.iterations)
        self.deltaXAsNsEas=self.MSI_instance_two.deltaXAsNsEas
        self.physical_obervable_updates_list = self.MSI_instance_two.physical_obervable_updates_list
        self.absorbance_observables_updates_list = self.MSI_instance_two.absorbance_coef_update_dict
        self.Ydf = self.MSI_instance_two.Y_data_frame
        self.Zdf = self.MSI_instance_two.z_data_frame
        self.experimental_dicts = self.MSI_instance_two.experiment_dictonaries
        self.z_matrix = self.MSI_instance_two.z_matrix
        self.s_matrix = self.MSI_instance_two.s_matrix
        self.y = self.MSI_instance_two.y_matrix
        self.Y_matrix = self.MSI_instance_two.Y_matrix
        self.S_matrix = self.MSI_instance_two.S_matrix
        
        self.X = self.MSI_instance_two.X
        self.Xdf = self.MSI_instance_two.X_data_frame
        self.covarience = self.MSI_instance_two.covarience
        self.exp_dict_list_optimized = self.MSI_instance_two.experiment_dictonaries
        self.parsed_yaml_list = self.MSI_instance_two.list_of_parsed_yamls
        self.sigma = self.MSI_instance_two.sigma
        self.X = self.MSI_instance_two.X
        self.delta_X = self.MSI_instance_two.delta_X
        self.molecular_parameter_updates=self.MSI_instance_two.delta_x_molecular_params_by_reaction_dict
        self.original_diag=np.diag(self.original_covariance)
        self.original_cti_file=self.MSI_instance_two.data_directory+'/'+self.MSI_instance_two.cti_file_name
        self.experiment_dict_uncertainty=self.MSI_instance_two.experiment_dict_uncertainty_original
        self.target_value_csv=self.MSI_instance_two.data_directory+'/'+self.MSI_instance_two.k_target_values_csv
        if self.targets[iteration]:
            self.k_target_value_S_matrix = self.MSI_instance_two.k_target_values_for_S
            self.run_with_k_target_values='On'
        else:
            self.k_target_value_S_matrix = None
            self.run_with_k_target_values='Off'
        self.csv_file_sigma=''
    def plotting_master(self,iteration,file_identifier=''):
        plotting_instance = plotter.Plotting(self.S_matrix,
                                     self.s_matrix,
                                     self.Y_matrix,
                                     self.y,
                                     self.z_matrix,
                                     self.X,
                                     self.sigma,
                                     self.covarience,
                                     self.original_covariance,
                                     self.S_matrix_original,
                                     self.exp_dict_list_optimized,
                                     self.exp_dict_list_original,
                                     self.parsed_yaml_list,
                                     self.Ydf,
                                     target_value_rate_constant_csv= self.optional_targets[iteration],
                                     target_value_rate_constant_csv_extra_values = self.optional_targets[iteration],
                                     k_target_value_S_matrix =self.k_target_value_S_matrix,
                                     k_target_values=self.run_with_k_target_values,
                                     working_directory = self.wdir,
                                     sigma_uncertainty_weighted_sensitivity_csv=self.csv_file_sigma,
                                     cheby_sensitivity_dict = self.cheb_coeffs[iteration],
                                     mapped_to_alpha_full_simulation=self.MSI_instance_two.mapped_to_alpha_full_simulation)

        #csv_file_sigma = MSI_st_instance_two.data_directory +'/'+'sigma_for_uncertainty_weighted_sensitivity_updated.csv'
        self.observable_counter_and_absorbance_wl,self.length_of_experimental_data = plotting_instance.lengths_of_experimental_data()
        self.sigmas_optimized,self.test = plotting_instance.calculating_sigmas(self.S_matrix,self.covarience)
        self.sigmas_original,self.test2 = plotting_instance.calculating_sigmas(self.S_matrix_original,self.original_covariance)
        plotting_instance.plotting_observables(sigmas_original = self.sigmas_original,sigmas_optimized= self.sigmas_optimized)
        self.diag = plotting_instance.getting_matrix_diag(self.covarience)
        
        plotting_instance.plotting_rate_constants(optimized_cti_file=self.MSI_instance_two.new_cti_file,
                                original_cti_file=self.original_cti_file,
                                initial_temperature=250,
                                final_temperature=2500,
                                master_equation_reactions = self.master_equation_reactions[iteration])



        self.sensitivity, self.top_sensitivity = plotting_instance.sort_top_uncertainty_weighted_sens()
        self.obs = plotting_instance.plotting_uncertainty_weighted_sens()
        
        
        
    def get_chebyshev_coefficients(self,iteration):
        self.cheb_coeffs=[]
        self.cheb_sensitivity_dict=None
        import importlib.util
        cheb=importlib.util.spec_from_file_location(self.master_sens[iteration].strip('.py'),self.wdir+'\\'+self.master_sens[iteration])
        #print(self.master_sens[iteration].strip('.py'))
        module=importlib.util.module_from_spec(cheb)
        cheb.loader.exec_module(module)
        #print(module.cheb_sensitivity_dict)
        #import self.master_sens[iteration].strip('.py') as c
        # file=open(self.master_sens[iteration],'r')
        
        # with open(self.master_sens[iteration],'r') as f:
        #     contents=f.read()
        #     dictionary=ast.literal_eval(contents)
        
        # print(dictionary)
        self.cheb_coeffs.append(module.cheb_sensitivity_dict)
        
        
        
        
        
    def msi_no_master(self,iteration,yaml_input=[]):
        MSI_instance_one=shell.MSI_optimization(self.models[iteration],
                                                0.01,
                                                1,
                                                1,
                                                self.wdir,
                                                yaml_input,
                                                self.reaction_uncertainties[iteration],
                                                k_target_values_csv=self.targets[iteration])
        MSI_instance_one.one_run_optimization()
        self.S_matrix_original=MSI_instance_one.S_matrix
        self.exp_dict_list_original=MSI_instance_one.experiment_dictonaries
        self.original_covariance=MSI_instance_one.covarience
        self.X_one_iteration=MSI_instance_one.X
        self.z_df_original=MSI_instance_one.z_data_frame
        
        self.MSI_instance_two=shell.MSI_optimization(self.models[iteration],
                                                0.01,
                                                1,
                                                1,
                                                self.wdir,
                                                yaml_input,
                                                self.reaction_uncertainties[iteration],
                                                k_target_values_csv=self.targets[iteration])
        self.X_list=self.MSI_instance_two.multiple_runs(self.iterations)
        
        self.deltaXAsNsEas=self.MSI_instance_two.deltaXAsNsEas
        self.physical_obervable_updates_list = self.MSI_instance_two.physical_obervable_updates_list
        self.absorbance_observables_updates_list = self.MSI_instance_two.absorbance_coef_update_dict
        self.Ydf = self.MSI_instance_two.Y_data_frame
        self.Zdf = self.MSI_instance_two.z_data_frame
        self.experimental_dicts = self.MSI_instance_two.experiment_dictonaries
        self.z_matrix = self.MSI_instance_two.z_matrix
        self.s_matrix = self.MSI_instance_two.s_matrix
        self.y = self.MSI_instance_two.y_matrix
        self.Y_matrix = self.MSI_instance_two.Y_matrix
        self.S_matrix = self.MSI_instance_two.S_matrix
        
        self.X = self.MSI_instance_two.X
        self.Xdf = self.MSI_instance_two.X_data_frame
        self.covarience = self.MSI_instance_two.covarience
        self.exp_dict_list_optimized = self.MSI_instance_two.experiment_dictonaries
        self.parsed_yaml_list = self.MSI_instance_two.list_of_parsed_yamls
        self.sigma = self.MSI_instance_two.sigma
        self.X = self.MSI_instance_two.X
        self.delta_X = self.MSI_instance_two.delta_X
        #target_value_rate_constant_csv = 'MSI/data/test_data/FFCM1_custom_target_value_test.csv'
        self.original_cti_file = self.MSI_instance_two.data_directory +'/'+ self.MSI_instance_two.cti_file_name
        self.experiment_dict_uncertainty = self.MSI_instance_two.experiment_dict_uncertainty_original
        self.target_value_csv = self.MSI_instance_two.data_directory +'/'+ self.MSI_instance_two.k_target_values_csv
        
        if self.targets[iteration]:
            self.k_target_value_S_matrix = self.MSI_instance_two.k_target_values_for_s
            self.run_with_k_target_values='On'
        else:
            self.k_target_value_S_matrix = None
            self.run_with_k_target_values='Off'
            
    def plotting_no_master(self, iteration,file_identifier=''):
        plotting_instance = plotter.Plotting(self.S_matrix,
                                          self.s_matrix,
                                          self.Y_matrix,
                                          self.Y_matrix,
                                          self.z_matrix,
                                          self.X,
                                          self.sigma,
                                          self.covarience,
                                          self.original_covariance,
                                          self.S_matrix_original,
                                          self.exp_dict_list_optimized,
                                          self.exp_dict_list_original,
                                          self.parsed_yaml_list,
                                          self.Ydf,
                                          target_value_rate_constant_csv= self.optional_targets[iteration],
                                          k_target_value_S_matrix =self.k_target_value_S_matrix,
                                          k_target_values=self.run_with_k_target_values,
                                          working_directory=self.wdir,
                                          shock_tube_instance = self.MSI_instance_two,
                                          optimized_cti_file=self.MSI_instance_two.new_cti_file,
                                          original_cti_file=self.original_cti_file)
    
        self.observable_counter_and_absorbance_wl,self.length_of_experimental_data = plotting_instance.lengths_of_experimental_data()
        self.sigmas_optimized,test = plotting_instance.calculating_sigmas(self.S_matrix,self.covarience)
        self.sigmas_original,self.test2 = plotting_instance.calculating_sigmas(self.S_matrix_original,self.original_covariance)
        plotting_instance.plotting_observables(sigmas_original = self.sigmas_original,sigmas_optimized= self.sigmas_optimized,
                                               file_identifier=self.models[iteration].rstrip('.cti')+'_'+file_identifier,
                                               filetype='.pdf')
        self.diag = plotting_instance.getting_matrix_diag(self.covarience)
        
        #plotting_instance.Y_matrix_plotter(Y_matrix,exp_dict_list_optimized,y,sigma)
        
        
        
        plotting_instance.plotting_rate_constants(optimized_cti_file=self.MSI_instance_two.new_cti_file,
                                        original_cti_file=self.original_cti_file,
                                        initial_temperature=250,
                                        final_temperature=2500)
                                        
        
        
        self.sensitivity, self.top_sensitivity = plotting_instance.sort_top_uncertainty_weighted_sens()
        self.obs = plotting_instance.plotting_uncertainty_weighted_sens()
        
        
        
    def get_parameters(self):
        self.get_working_directory()
        self.get_iterations()
        self.get_yaml_option()
        self.get_yaml_files()
        self.get_kinetic_models()
        self.get_reaction_uncertainties_list()
        self.get_master_equation_models()
        self.get_master_equation_uncertainties()
        self.get_targets()
        self.get_optional_plotting_targets()
        self.get_master_equation_sens()
        self.get_master_reactions()
        self.get_master_index()
        
    def get_working_directory(self):
        self.wdir=self.input_options[0].split('=')[1].rstrip('\'').lstrip('\'')
    
    def get_iterations(self):
        self.iterations=int(self.input_options[1].split('=')[1])
        
    def get_yaml_option(self):
        if re.match('[Ff][Aa][Ll][Ss][Ee]',self.input_options[2].split('=')[1]):
            self.run_individual_yaml=False
        elif re.match('[Tt][Rr][Uu][Ee]',self.input_options[2].split('=')[1]):
            self.run_individual_yaml=True
        else:
            raise Exception('Please set run_individual_yaml=True/False in the input file to run.')
        
    def get_yaml_files(self):
        yaml_bool=False
        self.yaml_files=[]
        for i,string in enumerate(self.input_options):
            if 'end_yaml' in string:
                break
            if yaml_bool:
                temp_entry=string.lstrip('[')
                temp_entry=temp_entry.rstrip(']')
                entry=temp_entry.split(',')
                #print(entry)
                if entry[1]=='':
                    entry=[entry[0]]
                self.yaml_files.append(entry)               
                
                
            if 'begin_yaml_list' in string:
                yaml_bool=True
            
            
    def get_kinetic_models(self):
        model_bool=False
        self.models=[]
        for i,string in enumerate(self.input_options):
            if 'end_model_list' in string:
                break
            if model_bool:
                self.models.append(string)
            if 'begin_model_list' in string:
                model_bool=True
            
    def get_master_index(self):
        index_bool=False
        self.indices=[]
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_index' in string:
                break
            if index_bool:
                self.indices.append(list(map(int, string.lstrip('[').rstrip(']').split(','))))
            if 'begin_master_equation_index' in string:
                index_bool=True
                
        
    def get_reaction_uncertainties_list(self):
        uncertainties_bool=False
        self.reaction_uncertainties=[]
        for i,string in enumerate(self.input_options):
            if 'end_reaction_uncertainty_list' in string:
                break
            if uncertainties_bool:
                self.reaction_uncertainties.append(string)
            if 'begin_reaction_uncertainty_list' in string:
                uncertainties_bool=True
    def get_master_equation_models(self):
        master_bool=False
        self.master_equation_models=[]
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_model_list' in string:
                break
            if master_bool:
                self.master_equation_models.append(string)
            if 'begin_master_equation_model_list' in string:
                master_bool=True
                
                
    def get_master_equation_uncertainties(self):
        master_bool=False
        self.master_uncertainties=[]
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_uncertainties' in string:
                break
            if master_bool:
                self.master_uncertainties.append(string)
            if 'begin_master_equation_uncertainties' in string:
                master_bool=True
                
    def get_targets(self):
        targets_bool=False
        self.targets=[]
        for i,string in enumerate(self.input_options):  
            if 'end_rate_constant_targets' in string:
                break
            if targets_bool:
                self.targets.append(string)
            if 'begin_rate_constant_targets' in string:
                targets_bool=True
        if not self.targets:
            self.targets=['']*len(self.models)
        
    def get_optional_plotting_targets(self):
        targets_bool=False
        self.optional_targets=[]
        for i,string in enumerate(self.input_options):  
            if 'end_optional_plotting_targets' in string:
                break
            if targets_bool:
                self.optional_targets.append(string)
            if 'begin_optional_plotting_targets' in string:
                targets_bool=True
        if not self.optional_targets and all('' == s or s.isspace() for s in self.targets):
    
            self.optional_targets=['']*len(self.models)
        elif not self.optional_targets and not all('' == s or s.isspace() for s in self.targets):
            self.optional_targets=self.targets
                
                
                
    
    def get_master_equation_sens(self):
        sens_bool=False
        self.master_sens=[]
        for i,string in enumerate(self.input_options):  
            if 'end_master_equation_sensitivities' in string:
                break
            if sens_bool:
                self.master_sens.append(string)
            if 'begin_master_equation_sensitivities' in string:
                sens_bool=True
            
            
    def get_master_reactions(self):
        rxn_list_len=len(self.models)
        set_names=[]
        reactions_bool=False
        set_bools=[False]*rxn_list_len
        self.master_equation_reactions=[]
        for i in range(rxn_list_len):
            set_names.append('set'+str(i+1))
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_reactions' in string:
                break
            if reactions_bool:
                for j,string1 in enumerate(set_bools):
                    #print(j)
                    index_top=self.input_options.index('set'+str(j+1))
                    index_bottom=self.input_options.index('end set'+str(j+1))
                    self.master_equation_reactions.append(self.input_options[index_top+1:index_bottom])
                break
                
                self.master_sens.append(string)
            if 'begin_master_equation_reactions' in string:
                reactions_bool=True
        


def parser(input_file):
    #print(#input_file)
    
    
    with open(input_file,'r') as f:
        input_lines=f.readlines()
        
    for i in range(len(input_lines)):
        #rint('char='+str(input_lines[i]))
        input_lines[i]=input_lines[i].lstrip()
        input_lines[i]=input_lines[i].rstrip('\n')
        input_lines[i]=input_lines[i].rstrip()
        if input_lines[i]:
            if input_lines[i][0]=='#':
               input_lines[i]=''
        
        
    while '' in input_lines:
        #print('a')
        input_lines.remove('')
    
      
    return input_lines
    
        



def main(input_file=''):

    if input_file=='':
        print('Please run program with defined input file using --input_file=FILEPATH')

    elif input_file !='':
        input_options=parser(input_file)
        simulation=multiscale_informatics(input_options)
        simulation.run_msi()
        simulation.write_convergence()
        
        return simulation

if __name__ == '__main__':
    a=fire.Fire(main)