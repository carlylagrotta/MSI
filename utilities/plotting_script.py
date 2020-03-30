import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cantera as ct
import copy
from textwrap import wrap
import scipy.stats as stats
import math
from scipy.stats import multivariate_normal
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import MSI.master_equation.master_equation as meq 
import re
import os






class Plotting(object):
    def __init__(self,S_matrix,
                 s_matrix,
                 Y_matrix,
                 y_matrix,
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
                 target_value_rate_constant_csv='',
                 target_value_rate_constant_csv_extra_values = '',
                 k_target_value_S_matrix = None,
                 k_target_values='Off',
                 working_directory='',
                 sigma_uncertainty_weighted_sensitivity_csv='',
                 simulation_run=None,
                 shock_tube_instance = None,
                 cheby_sensitivity_dict = None,
                 mapped_to_alpha_full_simulation=None):
        self.S_matrix = S_matrix
        self.s_matrix = s_matrix
        self.Y_matrix = Y_matrix
        self.y_matrix = y_matrix
        self.z_matrix = z_matrix
        self.X = X
        self.sigma = sigma
        #self.sigma = sigma
        self.covarience=covarience
        self.original_covariance=original_covariance
        #original
        self.S_matrix_original=S_matrix_original
        self.exp_dict_list_optimized = exp_dict_list_optimized
        self.exp_dict_list_original = exp_dict_list_original
        self.parsed_yaml_list = parsed_yaml_list
        self.target_value_rate_constant_csv = target_value_rate_constant_csv
        self.k_target_value_S_matrix = k_target_value_S_matrix
        self.Ydf = Ydf
        self.k_target_values=k_target_values
        self.target_value_rate_constant_csv_extra_values = target_value_rate_constant_csv_extra_values
        self.working_directory = working_directory
        self.sigma_uncertainty_weighted_sensitivity_csv  = sigma_uncertainty_weighted_sensitivity_csv
        self.simulation_run = simulation_run
        self.shock_tube_instance = shock_tube_instance
        self.cheby_sensitivity_dict=cheby_sensitivity_dict
        self.mapped_to_alpha_full_simulation = mapped_to_alpha_full_simulation
        
 #fix all the indexing to have a captial or lowercase time situation or add the module that lets you do either to all the scripts  

    def lengths_of_experimental_data(self):
        simulation_lengths_of_experimental_data = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            length_of_experimental_data=[]
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if observable in exp['mole_fraction_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Time'].shape[0])
                        observable_counter+=1
                    if re.match('[Jj][Ss][Rr]',exp['simulation_type']):
                        length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Temperature'].shape[0])
                        observable_counter+=1
                if observable in exp['concentration_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Time'].shape[0])
                        observable_counter+=1
                    if re.match('[Jj][Ss][Rr]',exp['simulation_type']):
                        length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Temperature'].shape[0])
                        observable_counter+=1
            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                absorbance_wl=0
                for k,wl in enumerate(wavelengths):
                    length_of_experimental_data.append(exp['absorbance_experimental_data'][k]['time'].shape[0])
                    absorbance_wl+=1
            else:
                absorbance_wl=0
                    
            simulation_lengths_of_experimental_data.append(length_of_experimental_data)
            
                    
        self.simulation_lengths_of_experimental_data=simulation_lengths_of_experimental_data
        
        return observable_counter+absorbance_wl,length_of_experimental_data
                    

        
    def calculating_sigmas(self,S_matrix,covarience):  
        sigmas =[[] for x in range(len(self.simulation_lengths_of_experimental_data))]
                 
        counter=0
        for x in range(len(self.simulation_lengths_of_experimental_data)):
            for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                temp=[]
                for z in np.arange(counter,(self.simulation_lengths_of_experimental_data[x][y]+counter)):       
                    SC = np.dot(S_matrix[z,:],covarience)
                    sigma = np.dot(SC,np.transpose(S_matrix[z,:]))
                    test = sigma
                    sigma = np.sqrt(sigma)
                    temp.append(sigma)
                temp = np.array(temp)            
                sigmas[x].append(temp)
        
                
                counter = counter + self.simulation_lengths_of_experimental_data[x][y]
        
        return sigmas, test
    
    
    
    def plotting_observables(self,sigmas_original=[],sigmas_optimized=[]):
        
        
        
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                plt.figure()
                
                if observable in exp['mole_fraction_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable],'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original['simulation'].timeHistories[0][observable],'r',label= "$\it{A priori}$ model")
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable],'o',color='black',label='Experimental Data')
                        plt.xlabel('Time (ms)')
                        plt.ylabel('Mole Fraction '+''+str(observable))
                        plt.title('Experiment_'+str(i+1))
                        
                        
                        
                        
    
                        
                        if bool(sigmas_optimized) == True:
                            
                            high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                            high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                            low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                            low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,'b--')
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,'b--')
                            
                            
                            
                            high_error_original = np.exp(sigmas_original[i][observable_counter])
                            high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                            low_error_original = np.exp(sigmas_original[i][observable_counter]*-1)
                            low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                            
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'r--')
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'r--')
                        
                        plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight',dpi=1000)
                        
    
                        observable_counter+=1
                    elif re.match('[Jj][Ss][Rr]',exp['simulation_type']):
                        plt.plot(exp['simulation'].timeHistories[0]['temperature'],exp['simulation'].timeHistories[0][observable],'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['temperature'],self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable],'r',label= "$\it{A priori}$ model")
                        plt.plot(exp['experimental_data'][observable_counter]['Temperature'],exp['experimental_data'][observable_counter][observable],'o',color='black',label='Experimental Data')
                        plt.xlabel('Temperature (K)')
                        plt.ylabel('Mole Fraction '+''+str(observable))
                        plt.title('Experiment_'+str(i+1))
                        
                        if bool(sigmas_optimized) == True:
                            
                            high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                            high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistories[0][observable].dropna().values)
                            low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                            low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistories[0][observable].dropna().values)
                            #plt.figure()
                            plt.plot(exp['experimental_data'][observable_counter]['Temperature'],  high_error_optimized,'b--')
                            plt.plot(exp['experimental_data'][observable_counter]['Temperature'],low_error_optimized,'b--')
                            
                            
                            
                            #high_error_original = np.exp(sigmas_original[i][observable_counter])
                           # high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                            #low_error_original = np.exp(sigmas_original[i][observable_counter]*-1)
                            #low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                            #plt.figure()
                           # plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'r--')
                            #plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'r--')
                        
                        plt.savefig(os.path.join(self.working_directory,'Experiment_'+str(i+1)+'_'+str(observable)+'.pdf'), bbox_inches='tight',dpi=1000)
                        observable_counter+=1
                if observable in exp['concentration_observables']:
                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable]*1e6,'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable]*1e6,'r',label= "$\it{A priori}$ model")
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Experimental Data') 
                    plt.xlabel('Time (ms)')
                    plt.ylabel('ppm'+''+str(observable))
                    plt.title('Experiment_'+str(i+1))
                    
                    if bool(sigmas_optimized)==True:
                        high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                        high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        low_error_optimized = np.exp(np.array(sigmas_optimized[i][observable_counter])*-1)
                        low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,'b--')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,'b--')                    
                        
    
    
                        high_error_original = np.exp(sigmas_original[i][observable_counter])
                        high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        low_error_original = np.exp(np.array(sigmas_original[i][observable_counter])*-1)
                        low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        
                        #plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'r--')
                        #plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'r--')
                    
                    plt.plot([],'w' ,label= 'T:'+ str(self.exp_dict_list_original[i]['simulation'].temperature))
                    plt.plot([],'w', label= 'P:'+ str(self.exp_dict_list_original[i]['simulation'].pressure))
                    key_list = []
                    for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():
                        
                        plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                        key_list.append(key)
                   
                    #plt.legend(handlelength=3)
                    plt.legend(ncol=2)
                    sp = '_'.join(key_list)
                    #print(sp)
                    #plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K'+'_'+str(self.exp_dict_list_original[i]['simulation'].pressure)+'_'+sp+'_'+'.pdf', bbox_inches='tight')
                    
                    #stub
                    plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                    plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)
                    


                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                plt.figure()
                for k,wl in enumerate(wavelengths):
                    plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')

                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['absorbance_calculated_from_model'][wl],'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['absorbance_calculated_from_model'][wl],'r',label= "$\it{A priori}$ model")
                    #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')
                    plt.xlabel('Time (ms)')
                    plt.ylabel('Absorbance'+''+str(wl))
                    plt.title('Experiment_'+str(i+1))
                    
                    if bool(sigmas_optimized)==True:
                        high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])
                        high_error_optimized = np.multiply(high_error_optimized,exp['absorbance_model_data'][wl])
                        low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                        low_error_optimized = np.multiply(low_error_optimized,exp['absorbance_model_data'][wl])
                        
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,high_error_optimized,'b--')
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,low_error_optimized,'b--')
                        
                        high_error_original = np.exp(sigmas_original[i][observable_counter])
                        high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['absorbance_model_data'][wl])
                        low_error_original =  np.exp(sigmas_original[i][observable_counter]*-1)
                        low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['absorbance_model_data'][wl])
                        
                        
                        #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,high_error_original,'r--')
                        #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,low_error_original,'r--')

#                    if bool(sigmas_optimized)==True and  i+1 == 11:    
#                        plt.ylim(top=.35)
                    
                    #start here
                    
                    plt.plot([],'w' ,label= 'T:'+ str(self.exp_dict_list_original[i]['simulation'].temperature))
                    plt.plot([],'w', label= 'P:'+ str(self.exp_dict_list_original[i]['simulation'].pressure))
                    for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():                        
                        plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                        

                    #plt.legend(handlelength=3)
                    plt.legend(ncol=2)
                    #plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')

                    plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                    plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)

                    
                    

# make function to plot rate constants 
                    
    def plotting_rate_constants(self,optimized_cti_file='',
                                original_cti_file='',
                                initial_temperature=250,
                                final_temperature=2500,
                                master_equation_reactions=[]):
        
        gas_optimized = ct.Solution(optimized_cti_file)
        gas_original = ct.Solution(original_cti_file)

        def unique_list(seq):
            checked = []
            for e in seq:
                if e not in checked:
                    checked.append(e)
            return checked
        
        def target_values_for_S(target_value_csv,
                                exp_dict_list,
                                S_matrix,
                                master_equation_reaction_list = [],
                                master_equation_sensitivites = {}):
                    
                    
                    
                target_value_csv = pd.read_csv(target_value_csv)
                target_reactions = target_value_csv['Reaction']
                target_temp = target_value_csv['temperature']
                target_press = target_value_csv['pressure']
                target_k = target_value_csv['k']
                reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
                number_of_reactions_in_cti = len(reactions_in_cti_file)
                As = []
                Ns =  []
                Eas = []
                    
    
                
                def create_empty_nested_reaction_list():
                    
                    
                    nested_reaction_list = [[] for x in range(len(master_equation_reaction_list))]
                    for reaction in master_equation_reaction_list:
                        for i,MP in enumerate(master_equation_sensitivites[reaction]):
                            nested_reaction_list[master_equation_reaction_list.index(reaction)].append(0)
                    return nested_reaction_list      
                
                
                def create_tuple_list(array_of_sensitivities):
                    tuple_list = []
                    for ix,iy in np.ndindex(array_of_sensitivities.shape):
                        tuple_list.append((ix,iy))
                    return tuple_list
                    
                MP_stack = []
                target_values_to_stack =  []
                for i,reaction in enumerate(target_reactions):
                    if reaction in master_equation_reaction_list:
                        nested_reaction_list = create_empty_nested_reaction_list()
                        for j, MP_array in enumerate(master_equation_sensitivites[reaction]):
                            tuple_list = create_tuple_list(MP_array)
                            temp = []
                            counter = 0    
                            for sensitivity in np.nditer(MP_array,order='C'):
                                k = tuple_list[counter][0]
                                l= tuple_list[counter][1]
                                counter +=1
                                   #need to add reduced p and t, and check these units were using to map
                                    
                                #these might not work
                                
                                t_alpha= meq.Master_Equation.chebyshev_specific_poly(self,k,meq.Master_Equation.calc_reduced_T(self,target_temp[i]))
                                
                                if target_press[i] ==0:
                                    target_press_new = 1e-9
                                else:
                                    target_press_new=target_press[i]
                                p_alpha = meq.Master_Equation.chebyshev_specific_poly(self,l,meq.Master_Equation.calc_reduced_P(self,target_press_new*101325))
                                #these might nowt work 
                                single_alpha_map = t_alpha*p_alpha*sensitivity
                                temp.append(single_alpha_map)
                            temp =sum(temp)
                            #should there be an = temp here 
                            #nested_reaction_list[master_equation_reaction_list.index(reaction)][j]=temp
                            nested_reaction_list[master_equation_reaction_list.index(reaction)][j]=temp
    
                        temp2  = nested_reaction_list
                        flat_list = [item for sublist in temp2 for item in sublist]
                        #print(flat_list)
                        MP_stack.append(nested_reaction_list)
                        flat_list = np.array(flat_list)
                        flat_list = flat_list.reshape((1,flat_list.shape[0])) 
                        target_values_to_stack.append(flat_list)
                    
                    else:
                        #this will need to get fixed if we want to handle all reactions as chevy
                        A_temp = np.zeros((1,number_of_reactions_in_cti-len(master_equation_reaction_list)))
        
                        N_temp = np.zeros((1,number_of_reactions_in_cti-len(master_equation_reaction_list)))
                        Ea_temp = np.zeros((1,number_of_reactions_in_cti-len(master_equation_reaction_list)))
                            #decide if this mapping is correct             
                        A_temp[0,reactions_in_cti_file.index(reaction)] = 1
                        N_temp [0,reactions_in_cti_file.index(reaction)] = np.log(target_temp[i])
                        Ea_temp[0,reactions_in_cti_file.index(reaction)] = (-1/target_temp[i])
                        
                        As.append(A_temp)
                        Ns.append(N_temp)
                        Eas.append(Ea_temp)
                        A_temp = A_temp.reshape((1,A_temp.shape[1]))
                        N_temp = N_temp.reshape((1,N_temp.shape[1]))
                        Ea_temp = Ea_temp.reshape((1,Ea_temp.shape[1]))
                        target_values_to_stack.append(np.hstack((A_temp,N_temp,Ea_temp)))
                        
                        
                   # might need to edit this to pass in s? and  
                S_matrix = S_matrix
                shape_s = S_matrix.shape
                S_target_values = []
                for i,row in enumerate(target_values_to_stack):
                    if target_reactions[i] in master_equation_reaction_list:
                        zero_to_append_infront = np.zeros((1,((number_of_reactions_in_cti-len(master_equation_reaction_list))*3)))
                        
                        zero_to_append_behind = np.zeros((1, shape_s[1] - ((number_of_reactions_in_cti-len(master_equation_reaction_list))*3) - np.shape(row)[1] ))                
                        temp_array = np.hstack((zero_to_append_infront,row,zero_to_append_behind))
                        S_target_values.append(temp_array)
                    else:
                        zero_to_append_behind = np.zeros((1,shape_s[1]-np.shape(row)[1]))
                        temp_array = np.hstack((row,zero_to_append_behind))
                        S_target_values.append(temp_array)
    
    
                S_target_values = np.vstack((S_target_values))
                return S_target_values      
            
        def sort_rate_constant_target_values(parsed_csv,unique_reactions,gas):
            reaction_list_from_mechanism = gas.reaction_equations()
            target_value_ks = [[] for reaction in range(len(unique_reactions))]
            target_value_temps = [[] for reaction in range(len(unique_reactions))]
            reaction_list_from_mechanism = gas.reaction_equations()
            
            for i,reaction in enumerate(parsed_csv['Reaction']):
                idx = reaction_list_from_mechanism.index(reaction)
                target_value_ks[unique_reactions.index(idx)].append(parsed_csv['k'][i])
                target_value_temps[unique_reactions.index(idx)].append(parsed_csv['temperature'][i])
                
            return target_value_temps,target_value_ks
        
        def rate_constant_over_temperature_range_from_cantera(reaction_number,
                                                              gas,
                                                              initial_temperature=250,
                                                              final_temperature=2500,
                                                              pressure=1,
                                                              conditions = {'H2':2,'O2':1,'N2':4}):
            Temp = []
            k = []
            
            for temperature in np.arange(initial_temperature,final_temperature,1):
                gas.TPX = temperature,pressure*101325,conditions
                Temp.append(temperature)
                k.append(gas.forward_rate_constants[reaction_number]*1000)
            return Temp,k
#
#        def calculate_sigmas_for_rate_constants(k_target_value_S_matrix,k_target_values_parsed_csv,unique_reactions,gas,covarience):
#
#            
#            reaction_list_from_mechanism = gas.reaction_equations()
#            sigma_list_for_target_ks = [[] for reaction in range(len(unique_reactions))]
#            shape = k_target_value_S_matrix.shape
#            for row in range(shape[0]):
#                SC = np.dot(k_target_value_S_matrix[row,:],covarience)
#                sigma_k = np.dot(SC,np.transpose(k_target_value_S_matrix[row,:]))
#                sigma_k = np.sqrt(sigma_k)
#                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
#                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
#                
#            return sigma_list_for_target_ks
        
        def calculate_sigmas_for_rate_constants(k_target_value_S_matrix,k_target_values_parsed_csv,unique_reactions,gas,covarience):

            
            reaction_list_from_mechanism = gas.reaction_equations()
            sigma_list_for_target_ks = [[] for reaction in range(len(unique_reactions))]
            shape = k_target_value_S_matrix.shape
            for row in range(shape[0]):
                #print(row)
                SC = np.dot(k_target_value_S_matrix[row,:],covarience)
                sigma_k = np.dot(SC,np.transpose(k_target_value_S_matrix[row,:]))
                sigma_k = np.sqrt(sigma_k)
                #print(row)
                #print(k_target_values_parsed_csv['Reaction'][row])
                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
                
            return sigma_list_for_target_ks
        
        def calculating_target_value_ks_from_cantera_for_sigmas(k_target_values_parsed_csv,gas,unique_reactions):
            target_value_ks = [[] for reaction in range(len(unique_reactions))]
            
            
            target_reactions = k_target_values_parsed_csv['Reaction']
            target_temp = k_target_values_parsed_csv['temperature']
            target_press = k_target_values_parsed_csv['pressure']
            reactions_in_cti_file = gas.reaction_equations()
            #print(reactions_in_cti_file)
            
            for i,reaction in enumerate(target_reactions): 

                if target_press[i] == 0:
                        pressure = 1e-9
                else:
                    pressure = target_press[i]
                        
                gas.TPX = target_temp[i],pressure*101325,{'H2O2':0.003094,'O2':0.000556,'H2O':0.001113,'Ar':0.995237}
                reaction_number_in_cti = reactions_in_cti_file.index(reaction)
                k = gas.forward_rate_constants[reaction_number_in_cti]
                indx = reactions_in_cti_file.index(reaction)
                target_value_ks[unique_reactions.index(indx)].append(k*1000)


            return target_value_ks
    
    
        if bool(self.target_value_rate_constant_csv) and self.k_target_values=='On':

                                
            S_matrix_k_target_values_extra = target_values_for_S(self.target_value_rate_constant_csv_extra_values,
                                                                 self.exp_dict_list_optimized,
                                                                 self.S_matrix,
                                                                 master_equation_reaction_list = master_equation_reactions,
                                                                 master_equation_sensitivites=self.cheby_sensitivity_dict)
            

#paste here
            
            unique_reactions_optimized=[]
            unique_reactions_original = []
            
            reaction_list_from_mechanism_original = gas_original.reaction_equations()
            reaction_list_from_mechanism = gas_optimized.reaction_equations()
            k_target_value_csv_extra = pd.read_csv(self.target_value_rate_constant_csv_extra_values)     
            k_target_value_csv = pd.read_csv(self.target_value_rate_constant_csv)
        
            for row in range(k_target_value_csv_extra.shape[0]):
                unique_reactions_optimized.append(reaction_list_from_mechanism.index(k_target_value_csv_extra['Reaction'][row]))
                unique_reactions_original.append(reaction_list_from_mechanism_original.index(k_target_value_csv_extra['Reaction'][row]))
            unique_reactions_optimized = unique_list(unique_reactions_optimized)
            unique_reactions_original = unique_list(unique_reactions_original)

            
            sigma_list_for_target_ks_optimized = calculate_sigmas_for_rate_constants(S_matrix_k_target_values_extra,k_target_value_csv_extra,unique_reactions_optimized,gas_optimized,self.covarience)
            self.sigma_list_for_target_ks_optimized = sigma_list_for_target_ks_optimized
            sigma_list_for_target_ks_original = calculate_sigmas_for_rate_constants(S_matrix_k_target_values_extra,k_target_value_csv_extra,unique_reactions_original,gas_original,self.original_covariance)
            self.sigma_list_for_target_ks_original = sigma_list_for_target_ks_original
            ######################  
            
            target_value_temps_optimized,target_value_ks_optimized = sort_rate_constant_target_values(k_target_value_csv_extra,unique_reactions_optimized,gas_optimized)
            target_value_temps_original,target_value_ks_original = sort_rate_constant_target_values(k_target_value_csv_extra,unique_reactions_original,gas_original)
           
            
            
            ############################################# 
            unique_reactions_optimized_for_plotting=[]
            unique_reactions_original_for_plotting = []
            
            for row in range(k_target_value_csv.shape[0]):
                unique_reactions_optimized_for_plotting.append(reaction_list_from_mechanism.index(k_target_value_csv['Reaction'][row]))
                unique_reactions_original_for_plotting.append(reaction_list_from_mechanism_original.index(k_target_value_csv['Reaction'][row]))
            
            unique_reactions_optimized_for_plotting = unique_list(unique_reactions_optimized)
            unique_reactions_original_for_plotting = unique_list(unique_reactions_original)            
            
            target_value_temps_optimized_for_plotting,target_value_ks_optimized_for_plotting = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_optimized_for_plotting,gas_optimized)
            target_value_temps_original_for_plotting,target_value_ks_original_for_plotting = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_original_for_plotting,gas_original)
           #############################################
           
           
            target_value_ks_calculated_with_cantera_optimized = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv_extra,gas_optimized,unique_reactions_optimized)
            target_value_ks_calculated_with_cantera_original = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv_extra,gas_original,unique_reactions_original)    
            
            
            for i,reaction in enumerate(unique_reactions_optimized):
                plt.figure()
                Temp_optimized,k_optimized = rate_constant_over_temperature_range_from_cantera(reaction,
                                                                  gas_optimized,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_optimized,k_optimized,'b')
                #calculate sigmas 
                high_error_optimized = np.exp(np.array(sigma_list_for_target_ks_optimized[i]))
                high_error_optimized = np.multiply(high_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])
                
                
                low_error_optimized = np.exp(np.array(sigma_list_for_target_ks_optimized[i])*-1)
                low_error_optimized = np.multiply(low_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])    
                
                #plt.semilogy(target_value_temps_optimized[i],high_error_optimized,'b--')   
                
                a, b = zip(*sorted(zip(target_value_temps_optimized[i],high_error_optimized)))
                plt.semilogy(a,b,'b--')
               # print(a,b)
                
                
                a, b = zip(*sorted(zip(target_value_temps_optimized[i],low_error_optimized)))  
                plt.semilogy(a,b,'b--')

                Temp_original,k_original = rate_constant_over_temperature_range_from_cantera(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]),
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_original,k_original,'r')

                high_error_original = np.exp(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])


                high_error_original = np.multiply(high_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])
                
                
                low_error_original = np.exp(np.array(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])*-1)
                low_error_original = np.multiply(low_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])  
                
                a, b = zip(*sorted(zip(target_value_temps_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))],high_error_original)))  
                plt.semilogy(a,b,'r--')
                
                a, b = zip(*sorted(zip(target_value_temps_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))],low_error_original)))  
                plt.semilogy(a,b,'r--')
                       
                #plt.semilogy(target_value_temps_optimized[i],target_value_ks_optimized[i],'o',color='black')
                plt.semilogy(target_value_temps_optimized_for_plotting[i],target_value_ks_optimized_for_plotting[i],'o',color='black')
                
                plt.xlabel('Temperature (K)')
                plt.ylabel('Kmol/m^3-s')
                plt.title(reaction_list_from_mechanism[reaction])
                #plt.savefig(self.working_directory+'/'+reaction_list_from_mechanism[reaction]+'.pdf', bbox_inches='tight')
                #plt.savefig(self.working_directory+'/'+reaction_list_from_mechanism[reaction]+'.svg', bbox_inches='tight')
                
        elif bool(self.target_value_rate_constant_csv) and self.k_target_values=='Off':
            
            unique_reactions_optimized=[]
            unique_reactions_original = []
            reaction_list_from_mechanism_original = gas_original.reaction_equations()
            reaction_list_from_mechanism = gas_optimized.reaction_equations()
            
            k_target_value_csv = pd.read_csv(self.target_value_rate_constant_csv)     
            for row in range(k_target_value_csv.shape[0]):
                unique_reactions_optimized.append(reaction_list_from_mechanism.index(k_target_value_csv['Reaction'][row]))
                unique_reactions_original.append(reaction_list_from_mechanism_original.index(k_target_value_csv['Reaction'][row]))
            unique_reactions_optimized = unique_list(unique_reactions_optimized)
            unique_reactions_original = unique_list(unique_reactions_original)
            
            
          ######################  
            target_value_temps_optimized,target_value_ks_optimized = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_optimized,gas_optimized)
            target_value_temps_original,target_value_ks_original = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_original,gas_original)
           ############################################# 
            target_value_ks_calculated_with_cantera_optimized = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv,gas_optimized,unique_reactions_optimized)
            target_value_ks_calculated_with_cantera_original = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv,gas_original,unique_reactions_original)
           
            for i,reaction in enumerate(unique_reactions_optimized):
                plt.figure()
                Temp_optimized,k_optimized = rate_constant_over_temperature_range_from_cantera(reaction,
                                                                  gas_optimized,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_optimized,k_optimized,'b')

                Temp_original,k_original = rate_constant_over_temperature_range_from_cantera(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]),
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_original,k_original,'r')

                
                plt.semilogy(target_value_temps_optimized[i],target_value_ks_optimized[i],'o',color='black')
                
                plt.xlabel('Temperature (K)')
                plt.ylabel('Kmol/m^3-s')
                plt.title(reaction_list_from_mechanism[reaction])
                plt.savefig(self.working_directory+'/'+reaction_list_from_mechanism[reaction]+'.pdf', bbox_inches='tight')
                plt.savefig(self.working_directory+'/'+reaction_list_from_mechanism[reaction]+'.svg', bbox_inches='tight')
                          
                
                
                
    def plotting_X_itterations(self,list_of_X_values_to_plot = [], list_of_X_array=[],number_of_iterations=None):
        for value in list_of_X_values_to_plot:
            temp = []
            for array in list_of_X_array:
                temp.append(array[value][0])
            
            plt.figure()
            plt.plot(np.arange(0,number_of_iterations,1),temp)
        
        return
        
                
                
        


    
    def getting_matrix_diag(self,cov_matrix):
        diag = cov_matrix.diagonal()
        return diag
            
            

                
    def Y_matrix_plotter(self,Y_matrix,exp_dict_list_optimized,y_matrix,sigma):  
        #sigmas =[[] for x in range(len(self.simulation_lengths_of_experimental_data))]
             
        counter=0
        for x in range(len(self.simulation_lengths_of_experimental_data)):
            observable_counter = 0 
            for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                #for z in np.arange(counter,(self.simulation_lengths_of_experimental_data[x][y]+counter)):       
                
                    
                plt.figure()   
                Y_values_to_plot = list(Y_matrix[counter:self.simulation_lengths_of_experimental_data[x][y]+counter,:])
                y_values_to_plot = list(y_matrix[counter:self.simulation_lengths_of_experimental_data[x][y]+counter,:])  
                sigmas_to_plot = list(sigma[counter:self.simulation_lengths_of_experimental_data[x][y]+counter,:])  
                if 'perturbed_coef' in exp_dict_list_optimized[x].keys():
                    wavelengths = self.parsed_yaml_list[x]['absorbanceCsvWavelengths'][0]
                    time = exp_dict_list_optimized[x]['absorbance_experimental_data'][0]['time']                     
                    plt.subplot(4, 1, 1)    
                    plt.title('Experiment_'+str(x+1)+'_Wavelength_'+str(wavelengths))
                    plt.plot(time*1e3,Y_values_to_plot)
                    plt.tick_params(labelbottom=False)
                    plt.ylabel('Y_matrix')
                    plt.subplot(plt.subplot(4, 1, 2))
                    plt.plot(time*1e3,y_values_to_plot)
                    plt.tick_params(labelbottom=False)
                    plt.ylabel('y_matrix')
                    plt.subplot(plt.subplot(4, 1, 3))
                    plt.plot(time*1e3,sigmas_to_plot)
                    plt.tick_params(labelbottom=False)
                    plt.ylabel('sigma')
                    plt.subplot(plt.subplot(4, 1, 4))
                    plt.plot(time*1e3,np.array(Y_values_to_plot)/np.array(sigmas_to_plot))
                    plt.ylabel('Y/sigma')
                    plt.xlabel('time')
                    
                    plt.savefig(self.working_directory+'/'+'Experiment_'+str(x+1)+' '+'Absorbance at'+'_'+str(wavelengths)+'.pdf', bbox_inches='tight')
                else:
                      
                    time = exp_dict_list_optimized[x]['experimental_data'][y]['Time']
                    plt.subplot(4, 1, 1)                  
                    plt.plot(time*1e3,Y_values_to_plot)
                    plt.tick_params(labelbottom=False)
                    plt.title('Experiment_'+str(x+1)+'_observable_'+exp_dict_list_optimized[0]['observables'][observable_counter])
                    plt.ylabel('Y_matrix')
                    plt.subplot(plt.subplot(4, 1, 2))
                    plt.plot(time*1e3,y_values_to_plot)
                    plt.tick_params(labelbottom=False)
                    plt.ylabel('y_matrix')
                    plt.subplot(plt.subplot(4, 1, 3))
                    plt.plot(time*1e3,sigmas_to_plot)
                    plt.tick_params(labelbottom=False)
                    plt.ylabel('sigma')
                    plt.subplot(plt.subplot(4, 1, 4))
                    plt.plot(time*1e3,np.array(Y_values_to_plot)/np.array(sigmas_to_plot))
                    plt.ylabel('Y/sigma')
                    plt.xlabel('time')
                    
                    plt.savefig('Experiment_'+str(x+1)+'_observable_'+exp_dict_list_optimized[0]['observables'][observable_counter]+'.pdf', bbox_inches='tight')
                    observable_counter+=1
        
                
                counter = counter + self.simulation_lengths_of_experimental_data[x][y]
        
        return   


    def shorten_sigma(self):
        flat_list = [item for sublist in self.simulation_lengths_of_experimental_data for item in sublist]
        length = sum(flat_list)
        observables_list = self.Ydf['value'].tolist()[length:]
        short_sigma = list(self.sigma)[length:]
        #print(flat_list)
        if bool(self.target_value_rate_constant_csv) and self.k_target_values=='On':
           
           k_target_value_csv = pd.read_csv(self.target_value_rate_constant_csv) 
           shape = k_target_value_csv.shape[0]
           slc = len(observables_list) - shape
           observables_list = observables_list[:slc]
           short_sigma = short_sigma[:slc]
           short_sigma = np.array(short_sigma)
        self.short_sigma =  short_sigma
           
        
        return 
            
    def sort_top_uncertainty_weighted_sens(self,top_sensitivity=10):
        S_matrix_copy = copy.deepcopy(self.S_matrix)
        self.shorten_sigma()
        sigma_csv = self.sigma_uncertainty_weighted_sensitivity_csv
        if bool(sigma_csv):
            df = pd.read_csv(sigma_csv)
            Sig = np.array(df['Sigma'])
            Sig = Sig.reshape((Sig.shape[0],1))
            
        else:
            Sig = self.short_sigma
            
        #Sig = self.sigma
        for pp  in range(np.shape(S_matrix_copy)[1]):
            S_matrix_copy[:,pp] *=Sig[pp]

        sensitivitys =[[] for x in range(len(self.simulation_lengths_of_experimental_data))]
        topSensitivities = [[] for x in range(len(self.simulation_lengths_of_experimental_data))]   
        start=0
        stop = 0
        for x in range(len(self.simulation_lengths_of_experimental_data)):
            for y in range(len(self.simulation_lengths_of_experimental_data[x])):           
    
                stop = self.simulation_lengths_of_experimental_data[x][y] + start
                temp = S_matrix_copy[start:stop,:]
                sort_s= pd.DataFrame(temp).reindex(pd.DataFrame(temp).abs().max().sort_values(ascending=False).index, axis=1)
                cc=pd.DataFrame(sort_s).iloc[:,:top_sensitivity]
                top_five_reactions=cc.columns.values.tolist()
                topSensitivities[x].append(top_five_reactions)
                #ccn=pd.DataFrame(cc).as_matrix()
                ccn=pd.DataFrame(cc).to_numpy()

                sensitivitys[x].append(ccn)           
                start = start + self.simulation_lengths_of_experimental_data[x][y]
                
               
        return sensitivitys,topSensitivities
    
    
    def getting_time_profiles_for_experiments(self, exp_dict_list_optimized):
        time_profiles =[[] for x in range(len(self.simulation_lengths_of_experimental_data))]
        observables = [[] for x in range(len(self.simulation_lengths_of_experimental_data))]
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue                                
                if observable in exp['mole_fraction_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        observables[i].append(observable)
                        observable_counter+=1
                    elif re.match('[Jj][Ss][Rr]',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Temperature'])
                        observables[i].append(observable)
                        observable_counter+=1
                if observable in exp['concentration_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Time']*1e3)        
                        observables[i].append(observable)                                
                        observable_counter+=1
                    elif re.match('[Jj][Ss][Rr]',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Temperature'])        
                        observables[i].append(observable)                                
                        observable_counter+=1
            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    time_profiles[i].append(exp['absorbance_experimental_data'][k]['time']*1e3)       
                    observables[i].append('Absorbance_'+str(wl))                                 
        self.time_profiles = time_profiles
        self.observable_list = observables
        return time_profiles
    
    
    def get_observables_list(self):
        #use this function to return observable list and uncertainty  pass in csv and get unc and csv
        sigma_csv = self.sigma_uncertainty_weighted_sensitivity_csv
        if bool(sigma_csv):
            df = pd.read_csv(sigma_csv)
            Sig = df['Sigma'].values
            Sig = np.array(Sig)
            Sig = Sig.reshape((Sig.shape[0],1))
            observable_list  = df['Observable'].tolist()
            self.sigma_list = Sig
            
            #print(self.sigma_list)
            return observable_list
            
        else:    
            flat_list = [item for sublist in self.simulation_lengths_of_experimental_data for item in sublist]
            length = sum(flat_list)
            observables_list = self.Ydf['value'].tolist()[length:]
            
            if bool(self.target_value_rate_constant_csv) and self.k_target_values=='On':
               
               k_target_value_csv = pd.read_csv(self.target_value_rate_constant_csv) 
               shape = k_target_value_csv.shape[0]
               slc = len(observables_list) - shape
               observables_list = observables_list[:slc]
        
           
        
            return observables_list
    
    
    def plotting_uncertainty_weighted_sens(self):
        sensitivities,top_sensitivities = self.sort_top_uncertainty_weighted_sens()
        observables_list = self.get_observables_list()
        if bool(self.sigma_uncertainty_weighted_sensitivity_csv):
            
            sigma_list = self.sigma_list
        else:
            sigma_list = list(self.short_sigma)
        #start here
        time_profiles = self.getting_time_profiles_for_experiments(self.exp_dict_list_optimized)
        list_of_experiment_observables = self.observable_list
        def subplot_function(number_of_observables_in_simulation,time_profiles,sensitivities,top_sensitivity_single_exp,observables_list,list_of_experiment_observables,experiment_number):
            #plt.figure(figsize=(2,6))
            #stub
            plt.figure()
            for plot_number in range(number_of_observables_in_simulation):
                for c,top_columns in enumerate(top_sensitivity_single_exp[plot_number]):
                    plt.subplot(number_of_observables_in_simulation,1,plot_number+1)
                    if plot_number==0:
                        plt.title('Experiment_'+str(experiment_number+1))
                    plt.plot(time_profiles[plot_number],sensitivities[plot_number][:,c],label = observables_list[top_columns] +'_'+str(sigma_list[top_columns])) 
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.5)
                    plt.ylabel(list_of_experiment_observables[plot_number])
                    top,bottom = plt.ylim()
                    left,right = plt.xlim()
                    plt.legend(ncol=5, loc='upper left',bbox_to_anchor=(-.5,-.3))
                    #plt.legend(ncol=3, loc='upper left',bbox_to_anchor=(1.2,2),fontsize=2)

            if self.simulation_run==None:
                
                plt.savefig(self.working_directory+'/'+'Experiment_'+str(experiment_number+1)+'.pdf', bbox_inches='tight')
            else:
                
                plt.title('Experiment_'+str(self.simulation_run))
                plt.savefig(self.working_directory+'/'+'Experiment_'+str(self.simulation_run)+'.pdf', bbox_inches='tight')

               
        for x in range(len(sensitivities)):            
            number_of_observables_in_simulation = len(sensitivities[x])
            subplot_function(number_of_observables_in_simulation,time_profiles[x],sensitivities[x],top_sensitivities[x],observables_list,list_of_experiment_observables[x],x)
            
        return 
            
            
         
            
    def plotting_rate_constants_six_paramter_fit(self,optimized_cti_file='',
                                original_cti_file='',
                                initial_temperature=250,
                                final_temperature=2500,
                                master_equation_reactions = [],
                                six_parameter_fit_dict_optimized = {},
                                six_parameter_fit_dict_nominal = {},
                                six_parameter_fit_sensitivity_dict = {}):
        
       
        gas_optimized = ct.Solution(optimized_cti_file)
        gas_original = ct.Solution(original_cti_file)
        
        def unique_list(seq):
            checked = []
            for e in seq:
                if e not in checked:
                    checked.append(e)
            return checked
        
################################################################################        

        def target_values_for_S_six_parameter_fit(target_value_csv,
                                    exp_dict_list,
                                    S_matrix,
                                    master_equation_reaction_list = [],
                                    six_parameter_fit_sensitivity_dict = {}):
                    
                    
                    
                  
                target_value_csv = pd.read_csv(target_value_csv)
                target_reactions = target_value_csv['Reaction']
                target_temp = target_value_csv['temperature']
                target_press = target_value_csv['pressure']
                target_k = target_value_csv['k']
                reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
                number_of_reactions_in_cti = len(reactions_in_cti_file)
                As = []
                Ns =  []
                Eas = []
                    
                Number_of_MP = []
                #nested_reaction_list = [[] for x in range(len(master_equation_reaction_list))]
                #print(six_parameter_fit_sensitivity_dict.keys())
                
                def create_empty_nested_reaction_list():
                    nested_reaction_list = [[] for x in range(len(master_equation_reaction_list))]
                    
                    for reaction in master_equation_reaction_list:
                        for i,MP in enumerate(six_parameter_fit_sensitivity_dict[reaction]['A']):
                            nested_reaction_list[master_equation_reaction_list.index(reaction)].append(0)
                    #copy.deepcopy(nested_reaction_list) 
                    #don't think i need this 
                    return nested_reaction_list          
                      
                MP_stack = []
                target_values_to_stack =  []
                for i,reaction in enumerate(target_reactions):
                    #temp_array = np.zeros((1,Number_of_MP))
                    if reaction in master_equation_reaction_list:
                        nested_reaction_list = create_empty_nested_reaction_list()
                        for s,sensitivity in enumerate(six_parameter_fit_sensitivity_dict[reaction]['A']):
                            #stub
                            #start here tomorrow 
                            nested_reaction_list[master_equation_reaction_list.index(reaction)][s] = 1*six_parameter_fit_sensitivity_dict[reaction]['A'][s] + np.log(target_temp[i])*six_parameter_fit_sensitivity_dict[reaction]['n'][s] + (-1000/target_temp[i])*six_parameter_fit_sensitivity_dict[reaction]['Ea'][s] + (-(1000/target_temp[i])**3)*six_parameter_fit_sensitivity_dict[reaction]['c'][s]+ (-(1000/target_temp[i])**-1)*six_parameter_fit_sensitivity_dict[reaction]['d'][s] + (-(1000/target_temp[i])**-3)*six_parameter_fit_sensitivity_dict[reaction]['f'][s]
                            #nested_reaction_list[master_equation_reaction_list.index(reaction)][s] = 1*six_parameter_fit_sensitivity_dict[reaction]['A'][s] + np.log(target_temp[i])*six_parameter_fit_sensitivity_dict[reaction]['n'][s] + (-1/target_temp[i])*six_parameter_fit_sensitivity_dict[reaction]['Ea'][s] + (-(1000/target_temp[i])**3)*six_parameter_fit_sensitivity_dict[reaction]['c'][s]+ (-(1000/target_temp[i])**-1)*six_parameter_fit_sensitivity_dict[reaction]['d'][s]*(1000*4.184)**-1 + (-(1/target_temp[i])**-3)*six_parameter_fit_sensitivity_dict[reaction]['f'][s]*(1000*4.184)**-3
                        temp  = nested_reaction_list
                        flat_list = [item for sublist in temp for item in sublist]
                        MP_stack.append(nested_reaction_list)
                        flat_list = np.array(flat_list)
                        flat_list = flat_list.reshape((1,flat_list.shape[0])) 
                        target_values_to_stack.append(flat_list)
                        
                                
                                          
                            
                    else:
                        A_temp = np.zeros((1,number_of_reactions_in_cti-len(master_equation_reaction_list)))
        
                        N_temp = np.zeros((1,number_of_reactions_in_cti-len(master_equation_reaction_list)))
                        Ea_temp = np.zeros((1,number_of_reactions_in_cti-len(master_equation_reaction_list)))
                            #decide if this mapping is correct             
                        A_temp[0,reactions_in_cti_file.index(reaction)] = 1
                        N_temp [0,reactions_in_cti_file.index(reaction)] = np.log(target_temp[i]) 
                        Ea_temp[0,reactions_in_cti_file.index(reaction)] = (-1/target_temp[i])
                        
                        As.append(A_temp)
                        Ns.append(N_temp)
                        Eas.append(Ea_temp)
                        A_temp = A_temp.reshape((1,A_temp.shape[1]))
                        N_temp = N_temp.reshape((1,N_temp.shape[1]))
                        Ea_temp = Ea_temp.reshape((1,Ea_temp.shape[1]))
                        target_values_to_stack.append(np.hstack((A_temp,N_temp,Ea_temp)))
                        
                        
                   # might need to edit this to pass in s? and  
                S_matrix = S_matrix
                shape_s = S_matrix.shape
                S_target_values = []
                for i,row in enumerate(target_values_to_stack):
                    if target_reactions[i] in master_equation_reaction_list:
                        zero_to_append_infront = np.zeros((1,((number_of_reactions_in_cti-len(master_equation_reaction_list))*3)))
                        
                        zero_to_append_behind = np.zeros((1, shape_s[1] - ((number_of_reactions_in_cti-len(master_equation_reaction_list))*3) - np.shape(row)[1] ))                
                        temp_array = np.hstack((zero_to_append_infront,row,zero_to_append_behind))
                        S_target_values.append(temp_array)
                    else:
                        zero_to_append_behind = np.zeros((1,shape_s[1]-np.shape(row)[1]))
                        temp_array = np.hstack((row,zero_to_append_behind))
                        S_target_values.append(temp_array)
    
    
                S_target_values = np.vstack((S_target_values))
                return S_target_values
        
################################################################################        
        
        
        def calculate_six_parameter_fit(reaction,dictonary,temperature):
        #finish editing this 
        #calc Ea,c,d,F seprately 
            A = dictonary[reaction]['A']
            n = dictonary[reaction]['n']
            Ea_temp = dictonary[reaction]['Ea']/(1.987*temperature)
            c_temp = dictonary[reaction]['c']/((1.987*temperature)**3)
            d_temp = dictonary[reaction]['d']*(1.987*temperature)
            f_temp = dictonary[reaction]['f']* ((1.987*temperature)**3)

            k = A*(temperature**n)*np.exp(-Ea_temp-c_temp-d_temp-f_temp)
            return k 
        def sort_rate_constant_target_values(parsed_csv,unique_reactions,gas):
            reaction_list_from_mechanism = gas.reaction_equations()
            target_value_ks = [[] for reaction in range(len(unique_reactions))]
            target_value_temps = [[] for reaction in range(len(unique_reactions))]
            reaction_list_from_mechanism = gas.reaction_equations()
            
            for i,reaction in enumerate(parsed_csv['Reaction']):
                idx = reaction_list_from_mechanism.index(reaction)
                target_value_ks[unique_reactions.index(idx)].append(parsed_csv['k'][i])
                target_value_temps[unique_reactions.index(idx)].append(parsed_csv['temperature'][i])
                
            return target_value_temps,target_value_ks
        def rate_constant_over_temperature_range_from_cantera(reaction_number,
                                                              gas,
                                                              initial_temperature=250,
                                                              final_temperature=2500,
                                                              pressure=1,
                                                              conditions = {'H2':2,'O2':1,'N2':4},
                                                              dictonary={},
                                                              master_equation_reactions=[]):
            Temp = []
            k = []
            
            
            reaction_string = gas.reaction_equations()[reaction_number] 
            for temperature in np.arange(initial_temperature,final_temperature,1):
                
                if reaction_string in master_equation_reactions:
                    k.append(calculate_six_parameter_fit(reaction_string,dictonary,temperature))
                    Temp.append(temperature)
                #start editing here
                else:
                
                    gas.TPX = temperature,pressure*101325,conditions
                    Temp.append(temperature)
                    k.append(gas.forward_rate_constants[reaction_number]*1000)
            return Temp,k

        def calculate_sigmas_for_rate_constants(k_target_value_S_matrix,k_target_values_parsed_csv,unique_reactions,gas,covarience):

            
            reaction_list_from_mechanism = gas.reaction_equations()
            sigma_list_for_target_ks = [[] for reaction in range(len(unique_reactions))]
            shape = k_target_value_S_matrix.shape
            for row in range(shape[0]):
                #print(row)
                SC = np.dot(k_target_value_S_matrix[row,:],covarience)
                sigma_k = np.dot(SC,np.transpose(k_target_value_S_matrix[row,:]))
                sigma_k = np.sqrt(sigma_k)
                #print(row)
                #print(k_target_values_parsed_csv['Reaction'][row])
                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
                
            return sigma_list_for_target_ks
        
        def calculating_target_value_ks_from_cantera_for_sigmas(k_target_values_parsed_csv,gas,unique_reactions,six_parameter_fit_dictonary,master_equation_reactions):
            target_value_ks = [[] for reaction in range(len(unique_reactions))]
            
            
            target_reactions = k_target_values_parsed_csv['Reaction']
            target_temp = k_target_values_parsed_csv['temperature']
            target_press = k_target_values_parsed_csv['pressure']
            reactions_in_cti_file = gas.reaction_equations()
            #print(reactions_in_cti_file)
            
            for i,reaction in enumerate(target_reactions): 
                if reaction in master_equation_reactions:
                    k = calculate_six_parameter_fit(reaction,six_parameter_fit_dictonary,target_temp[i])
                    indx = reactions_in_cti_file.index(reaction)
                    target_value_ks[unique_reactions.index(indx)].append(k)

                else:
                    if target_press[i] == 0:
                        pressure = 1e-9
                    else:
                        pressure = target_press[i]
                        
                    gas.TPX = target_temp[i],pressure*101325,{'H2O2':0.003094,'O2':0.000556,'H2O':0.001113,'Ar':0.995237}
                    reaction_number_in_cti = reactions_in_cti_file.index(reaction)
                    k = gas.forward_rate_constants[reaction_number_in_cti]
                    indx = reactions_in_cti_file.index(reaction)

                    target_value_ks[unique_reactions.index(indx)].append(k*1000)


            return target_value_ks
    
    
        if bool(self.target_value_rate_constant_csv) and self.k_target_values=='On':
            
            ### make new s matrix with the new csv file, and make sure we are plotting the old one
           
            S_matrix_k_target_values_extra = target_values_for_S_six_parameter_fit(self.target_value_rate_constant_csv_extra_values,
                                                                                   self.exp_dict_list_optimized,
                                                                                   self.S_matrix,
                                                                                   master_equation_reaction_list = master_equation_reactions,
                                                                                   six_parameter_fit_sensitivity_dict = six_parameter_fit_sensitivity_dict)
            
            

            
            #make two unique
            unique_reactions_optimized=[]
            unique_reactions_original = []
            
            reaction_list_from_mechanism_original = gas_original.reaction_equations()
            reaction_list_from_mechanism = gas_optimized.reaction_equations()
            k_target_value_csv_extra = pd.read_csv(self.target_value_rate_constant_csv_extra_values)     
            k_target_value_csv = pd.read_csv(self.target_value_rate_constant_csv)
            for row in range(k_target_value_csv_extra.shape[0]):
                unique_reactions_optimized.append(reaction_list_from_mechanism.index(k_target_value_csv_extra['Reaction'][row]))
                unique_reactions_original.append(reaction_list_from_mechanism_original.index(k_target_value_csv_extra['Reaction'][row]))
            unique_reactions_optimized = unique_list(unique_reactions_optimized)
            unique_reactions_original = unique_list(unique_reactions_original)

            
            sigma_list_for_target_ks_optimized = calculate_sigmas_for_rate_constants(S_matrix_k_target_values_extra,k_target_value_csv_extra,unique_reactions_optimized,gas_optimized,self.covarience)
            self.sigma_list_for_target_ks_optimized = sigma_list_for_target_ks_optimized
            sigma_list_for_target_ks_original = calculate_sigmas_for_rate_constants(S_matrix_k_target_values_extra,k_target_value_csv_extra,unique_reactions_original,gas_original,self.original_covariance)
            self.sigma_list_for_target_ks_original = sigma_list_for_target_ks_original
            ######################  
            target_value_temps_optimized,target_value_ks_optimized = sort_rate_constant_target_values(k_target_value_csv_extra,unique_reactions_optimized,gas_optimized)
            target_value_temps_original,target_value_ks_original = sort_rate_constant_target_values(k_target_value_csv_extra,unique_reactions_original,gas_original)
           
            
            
            ############################################# 
            unique_reactions_optimized_for_plotting=[]
            unique_reactions_original_for_plotting = []
            
            for row in range(k_target_value_csv.shape[0]):
                unique_reactions_optimized_for_plotting.append(reaction_list_from_mechanism.index(k_target_value_csv['Reaction'][row]))
                unique_reactions_original_for_plotting.append(reaction_list_from_mechanism_original.index(k_target_value_csv['Reaction'][row]))
            unique_reactions_optimized_for_plotting = unique_list(unique_reactions_optimized)
            unique_reactions_original_for_plotting = unique_list(unique_reactions_original)            
            
            target_value_temps_optimized_for_plotting,target_value_ks_optimized_for_plotting = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_optimized_for_plotting,gas_optimized)
            target_value_temps_original_for_plotting,target_value_ks_original_for_plotting = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_original_for_plotting,gas_original)
           #############################################
           
           
           
            target_value_ks_calculated_with_cantera_optimized = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv_extra,gas_optimized,unique_reactions_optimized,six_parameter_fit_dict_optimized,master_equation_reactions)
            target_value_ks_calculated_with_cantera_original = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv_extra,gas_original,unique_reactions_original,six_parameter_fit_dict_nominal,master_equation_reactions)
            
            
            
            #print(target_value_ks_calculated_with_cantera_original)
            
            
            
            for i,reaction in enumerate(unique_reactions_optimized):
                plt.figure()
                Temp_optimized,k_optimized = rate_constant_over_temperature_range_from_cantera(reaction,
                                                                  gas_optimized,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1.635,
                                                                  conditions={'H2O2':0.003094,'O2':0.000556,'H2O':0.001113,'Ar':0.995237},
                                                                  dictonary = six_parameter_fit_dict_optimized,
                                                                  master_equation_reactions = master_equation_reactions)
                
                plt.semilogy(Temp_optimized,k_optimized,'b')
                #calculate sigmas 
                #print(sigma_list_for_target_ks_optimized[i])
                high_error_optimized = np.exp(np.array(sigma_list_for_target_ks_optimized[i]))
                #print(high_error_optimized)
                high_error_optimized = np.multiply(high_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])
                
                
                low_error_optimized = np.exp(np.array(sigma_list_for_target_ks_optimized[i])*-1)
                low_error_optimized = np.multiply(low_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])    
                
               # plt.semilogy(target_value_temps_optimized[i],high_error_optimized,'b--')   
                
                a, b = zip(*sorted(zip(target_value_temps_optimized[i],high_error_optimized)))

                #plt.scatter(a,b,color='blue')
                
                plt.semilogy(a,b,'b--')
                
                a, b = zip(*sorted(zip(target_value_temps_optimized[i],low_error_optimized)))  
                plt.semilogy(a,b,'b--')
                #plt.scatter(a,b,color='blue')
               # print(a,b)
                Temp_original,k_original = rate_constant_over_temperature_range_from_cantera(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]),
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1.635,
                                                                  conditions={'H2O2':0.003094,'O2':0.000556,'H2O':0.001113,'Ar':0.995237},
                                                                  dictonary = six_parameter_fit_dict_nominal,
                                                                  master_equation_reactions = master_equation_reactions)
                
                plt.semilogy(Temp_original,k_original,'r')
               # plt.xlim((0,3000))
                #plt.ylim((10**9,10**15))
                #print(unique_reactions_original)
               # print(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))
                #print(unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction])))
                high_error_original = np.exp(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])


                high_error_original = np.multiply(high_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])
                
                
                low_error_original = np.exp(np.array(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])*-1)
                low_error_original = np.multiply(low_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))])  
                
                a, b = zip(*sorted(zip(target_value_temps_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))],high_error_original)))  
                plt.semilogy(a,b,'r--')
                #plt.scatter(a,b,color='red')
                
                
                a, b = zip(*sorted(zip(target_value_temps_original[unique_reactions_original.index(reaction_list_from_mechanism_original.index(reaction_list_from_mechanism[reaction]))],low_error_original)))  
                plt.semilogy(a,b,'r--')
                #plt.scatter(a,b,color='red')
                
                plt.semilogy(target_value_temps_optimized_for_plotting[i],target_value_ks_optimized_for_plotting[i],'o',color='black')
                
                plt.xlabel('Temperature [K]')
                #plt.ylabel('Kmol/m^3-s')
                plt.ylabel(r'k [$\frac{cm^3}{{mol s}}$]')

                plt.title(reaction_list_from_mechanism[reaction])
                plt.tick_params(axis ='both', direction ='in') 
                plt.tick_params(axis ='both', direction ='in',which='minor') 

                #plt.savefig(os.path.join(self.working_directory,reaction_list_from_mechanism[reaction]+'.pdf'), bbox_inches='tight')
                #plt.savefig(os.path.join(self.working_directory,reaction_list_from_mechanism[reaction]+'.svg'), bbox_inches='tight')

        elif bool(self.target_value_rate_constant_csv) and self.k_target_values=='Off':
            
            unique_reactions_optimized=[]
            unique_reactions_original = []
            reaction_list_from_mechanism_original = gas_original.reaction_equations()
            reaction_list_from_mechanism = gas_optimized.reaction_equations()
            
            k_target_value_csv = pd.read_csv(self.target_value_rate_constant_csv)     
            for row in range(k_target_value_csv.shape[0]):
                unique_reactions_optimized.append(reaction_list_from_mechanism.index(k_target_value_csv['Reaction'][row]))
                unique_reactions_original.append(reaction_list_from_mechanism_original.index(k_target_value_csv['Reaction'][row]))
            unique_reactions_optimized = unique_list(unique_reactions_optimized)
            unique_reactions_original = unique_list(unique_reactions_original)

            
            
            
          ######################  
            target_value_temps_optimized,target_value_ks_optimized = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_optimized,gas_optimized)
            target_value_temps_original,target_value_ks_original = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_original,gas_original)
           ############################################# 
            target_value_ks_calculated_with_cantera_optimized = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv,gas_optimized,unique_reactions_optimized)
            target_value_ks_calculated_with_cantera_original = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv,gas_original,unique_reactions_original)
           
            for i,reaction in enumerate(unique_reactions_optimized):
                plt.figure()
                Temp_optimized,k_optimized = rate_constant_over_temperature_range_from_cantera(reaction,
                                                                  gas_optimized,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_optimized,k_optimized,'b')

                    
                Temp_original,k_original = rate_constant_over_temperature_range_from_cantera(unique_reactions_original[unique_reactions_original.index(reaction)],
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_original,k_original,'r')
                
                
                plt.semilogy(target_value_temps_optimized[i],target_value_ks_optimized[i],'o',color='black')
                plt.xlabel('Temperature (K)')
                plt.ylabel('Kmol/m^3-s')
                plt.title(reaction_list_from_mechanism[reaction])
        return S_matrix_k_target_values_extra
                
                
                
        
    def plotting_normal_distributions(self,
                                      paramter_list,
                                      optimized_cti_file='',
                                      pdf_distribution_file='',
                                      shock_tube_instance=None):
        
        all_parameters = shock_tube_instance.posterior_diag_df['parameter'].tolist()
        df = shock_tube_instance.posterior_diag_df
        gas_optimized = ct.Solution(optimized_cti_file)
        
        for parameter in paramter_list:
            indx = all_parameters.index(parameter)
            variance = df['value'][indx]
            if parameter[0]=='A' or parameter[0]=='n' or parameter[0]=='E':
                letter,number = parameter.split('_')
                number = int(number)
                if 'ElementaryReaction' in str(type(gas_optimized.reaction(number))):
                    A=gas_optimized.reaction(number).rate.pre_exponential_factor
                    n=gas_optimized.reaction(number).rate.temperature_exponent
                    Ea=gas_optimized.reaction(number).rate.activation_energy
                if 'FalloffReaction' in str(type(gas_optimized.reaction(number))):
                    A=gas_optimized.reaction(number).high_rate.pre_exponential_factor
                    n=gas_optimized.reaction(number).high_rate.temperature_exponent
                    Ea=gas_optimized.reaction(number).high_rate.activation_energy
                if 'ThreeBodyReaction' in   str(type(gas_optimized.reaction(number))):
                    A=gas_optimized.reaction(number).rate.pre_exponential_factor
                    n=gas_optimized.reaction(number).rate.temperature_exponent
                    Ea=gas_optimized.reaction(number).rate.activation_energy
            else: 
                letter = None
                
            if letter =='A':
                mu = np.log(A*1000)
                sigma = math.sqrt(variance)
                sigma = sigma
                
            elif letter == 'n':
                mu = n
                sigma = math.sqrt(variance)
                #sigma = sigma/2
            elif letter == 'Ea':
                mu=Ea/1000/4.184            
                sigma = math.sqrt(variance)
                sigma = sigma*ct.gas_constant/(1000*4.184)
                #sigma = sigma/2
            else:
                mu= 0 
                sigma = math.sqrt(variance)
                

            
            
            x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
            plt.figure()
            plt.plot(x, stats.norm.pdf(x, mu, sigma))
            plt.xlabel(parameter)
            plt.ylabel('pdf')
            plt.savefig(self.working_directory+'/'+parameter+'_distribution'+'_.pdf',bbox_inches='tight')

            if bool(pdf_distribution_file):
                df2 = pd.read_csv(pdf_distribution_file)
                #temp = np.log(np.exp(df2[parameter].values)/9.33e13)
                #plt.plot(temp,df2['pdf_'+parameter])
                plt.plot(df2[parameter],df2['pdf_'+parameter])
                plt.savefig(self.working_directory+'/'+parameter+'_distribution'+'_.pdf',bbox_inches='tight')

    

    
    def plotting_joint_normal_distributions(self,
                                            coupled_parameters,
                                            optimized_cti_file='',
                                            joint_data_csv=''):
                
        all_parameters = self.shock_tube_instance.posterior_diag_df['parameter'].tolist()
        df = self.shock_tube_instance.posterior_diag_df
        gas_optimized = ct.Solution(optimized_cti_file)
        for couple in coupled_parameters:
            indx1 = all_parameters.index(couple[0])
            indx2 = all_parameters.index(couple[1])
            variance1 = df['value'][indx1]
            variance2 = df['value'][indx2]
            if couple[0][0]=='A' or couple[0][0]=='n' or couple[0][0]=='E':
           
                letter1,number1 = couple[0].split('_')
                number1 = int(number1)
                number1_covariance = number1
                if letter1=='n':
                    number1_covariance = number1+len(gas_optimized.reaction_equations())
                if letter1=='Ea':
                    number1_covariance = number1+len(gas_optimized.reaction_equations())*2
                    

                    
                if 'ElementaryReaction' in str(type(gas_optimized.reaction(number1))):
                    A1=gas_optimized.reaction(number1).rate.pre_exponential_factor
                    n1=gas_optimized.reaction(number1).rate.temperature_exponent
                    Ea1=gas_optimized.reaction(number1).rate.activation_energy
                if 'FalloffReaction' in str(type(gas_optimized.reaction(number1))):
                    A1=gas_optimized.reaction(number1).high_rate.pre_exponential_factor
                    n1=gas_optimized.reaction(number1).high_rate.temperature_exponent
                    Ea1=gas_optimized.reaction(number1).high_rate.activation_energy
                if 'ThreeBodyReaction' in   str(type(gas_optimized.reaction(number1))):
                    A1=gas_optimized.reaction(number1).rate.pre_exponential_factor
                    n1=gas_optimized.reaction(number1).rate.temperature_exponent
                    Ea1=gas_optimized.reaction(number1).rate.activation_energy
            else:
                letter1 = None
                mu1=0
                mu_x=0
                sigma1= math.sqrt(variance1)
                number1_covariance = indx1
                variance_x = variance1
            if couple[1][0]=='A' or couple[1][0]=='n' or couple[1][0]=='E':
                letter2,number2 = couple[1].split('_')
                number2 = int(number2)
                number2_covariance = number2
                if letter2=='n':
                    number2_covariance = number2+len(gas_optimized.reaction_equations())
                if letter2 == 'Ea':
                    number2_covariance = number2+len(gas_optimized.reaction_equations())*2
                    
                if 'ElementaryReaction' in str(type(gas_optimized.reaction(number2))):   
                    A2=gas_optimized.reaction(number2).rate.pre_exponential_factor
                    n2=gas_optimized.reaction(number2).rate.temperature_exponent
                    Ea2=gas_optimized.reaction(number2).rate.activation_energy 
                if 'FalloffReaction' in str(type(gas_optimized.reaction(number2))):
                    A2=gas_optimized.reaction(number2).high_rate.pre_exponential_factor
                    n2=gas_optimized.reaction(number2).high_rate.temperature_exponent
                    Ea2=gas_optimized.reaction(number2).high_rate.activation_energy
                if 'ThreeBodyReaction' in   str(type(gas_optimized.reaction(number2))):
                    A2=gas_optimized.reaction(number2).rate.pre_exponential_factor
                    n2=gas_optimized.reaction(number2).rate.temperature_exponent
                    Ea2=gas_optimized.reaction(number2).rate.activation_energy
            else:
                mu_y=0
                mu2=0
                letter2=None
                variance_y = variance2
                sigma = math.sqrt(variance2)
                number2_covariance = indx2

            
            
            covariance_couple = self.covarience[number1_covariance,number2_covariance]
           # print(number1_covariance,number2_covariance)
            #covariance_couple = .00760122
            if letter1 =='A':
                mu1 = np.log(A1*1000)
                mu_x = mu1
                variance_x = variance1
                sigma = np.sqrt(variance_x)
                
                #sigma = np.exp(sigma)
                #sigma = sigma*1000
                #sigma = np.log(sigma)
                #sigma = sigma/2
                variance_x = sigma**2
                #convert to chemkin units
            if letter1 == 'n':
                mu1 = n1
                mu_x = mu1
                variance_x = variance1
                sigma = np.sqrt(variance_x)
                #sigma = sigma/2
                variance_x = sigma**2
            if letter1 == 'Ea':
                mu1=Ea1/1000/4.184 
                mu_x = mu1                
                variance_x = variance1
                sigma = math.sqrt(variance_x)
                sigma = sigma*ct.gas_constant/(1000*4.184)
                #sigma = sigma/2
                variance_x = sigma**2
   
            
            if letter2 =='A':
                mu2 = np.log(A2*1000)
                mu_y = mu2
                variance_y = variance2
                sigma = np.sqrt(variance_y)
                sigma = sigma
                #sigma = np.exp(sigma)
                #sigma = sigma*1000
                #sigma = np.log(sigma)
                #sigma = sigma/2
                variance_y = sigma**2
                #convert to chemkin units
            if letter2 == 'n':
                mu2 = n2
                mu_y = mu2
                variance_y = variance2      
                sigma = np.sqrt(variance_y)
                #sigma = sigma/2
                variance_y = sigma**2
                
            if letter2 == 'Ea':
                mu2 = Ea2/1000/4.184 
                mu_y = mu2
                variance_y = variance2
                sigma = math.sqrt(variance_y)
                sigma = sigma*ct.gas_constant/(1000*4.184)
                #sigma = sigma/2
                variance_y = sigma**2

            
            if letter2 =='Ea' or letter1 == 'Ea':
                covariance_couple = covariance_couple*ct.gas_constant/(1000*4.184)
                if letter2=='Ea' and letter1=='Ea':
                    covariance_couple = np.sqrt(covariance_couple)
                    covariance_couple = covariance_couple*ct.gas_constant/(1000*4.184)
                    covariance_couple = covariance_couple**2

            #if letter1=='A' or letter2=='A':
                #covariance_couple = np.exp(covariance_couple)
                #covariance_couple  = covariance_couple/2
                #covariance_couple = np.log(covariance_couple)
                
                
           
            x = np.linspace(mu1 - 3*np.sqrt(variance_x), mu1 + 3*np.sqrt(variance_x),1000)
            y = np.linspace(mu2 - 3*np.sqrt(variance_y), mu2 + 3*np.sqrt(variance_y),1000)
            
            #x = np.linspace(mu1 - 2*np.sqrt(variance_x), mu1 + 2*np.sqrt(variance_x),1000)
            #y = np.linspace(mu2 - 2*np.sqrt(variance_y), mu2 + 2*np.sqrt(variance_y),1000)
            #TEST

            
            
            X,Y = np.meshgrid(x,y)
            #X, Y = np.meshgrid(x,y)


            
            
            pos = np.empty(X.shape + (2,))
            pos[:, :, 0] = X; pos[:, :, 1] = Y
            rv = multivariate_normal([mu_x, mu_y], [[variance_x, covariance_couple], [covariance_couple, variance_y]])
            print(couple,[mu_x, mu_y], [[variance_x, covariance_couple], [covariance_couple, variance_y]])
            fig = plt.figure()

            ax = fig.gca(projection='3d')
            ax.plot_surface(X, Y, rv.pdf(pos),cmap='viridis',linewidth=0)
            ax.set_xlabel(couple[0])
            ax.set_ylabel(couple[1])
            ax.set_zlabel('Z axis')
            plt.show()

            additional_dictionary = {'A_5':{'reaction':'H2O2 + M = 2OH + M','our_value':np.log(4.99999e8),'hong_value':np.log(5.60e8)},
                                     'A_6':{'reaction':'OH + H2O2 = H2O + HO2','our_value':np.log(5624842396127.52),'hong_value':np.log(6.93e12)},
                                     'A_7':{'reaction': 'OH + HO2 = H2O + O2' , 'our_value':np.log(16646221572429.6),'hong_value':np.log(1.82e13)},
                                     'A_8':{'reaction':'2HO2 = H2O2 + O2','our_value':np.log(806831822530.157),'hong_value':np.log(3.17e12)},
                                     'A_11':{'reaction':'2OH = H2O + O','our_value':np.log(1730749579423.63),'hong_value':np.log(2.355e12)},
                                     'Sigma_1':{'reaction':'sigma H2O2','our_value':-.03846,'hong_value':0},
                                     'Sigma_2':{'reaction':'sigma_HO2','our_value':.0721,'hong_value':0}}
           
            additional_dictionary = {'A_5':{'reaction':'H2O2 + M = 2OH + M','our_value':np.log(4.99999e8),'hong_value':np.log(5.60e8)},
                                     'A_6':{'reaction':'OH + H2O2 = H2O + HO2','our_value':np.log(5917630773605.197),'hong_value':np.log(6.93e12)},
                                     'A_7':{'reaction': 'OH + HO2 = H2O + O2' , 'our_value':np.log(18236369573049.9),'hong_value':np.log(1.82e13)},
                                     'A_8':{'reaction':'2HO2 = H2O2 + O2','our_value':np.log(863643827140.3533),'hong_value':np.log(3.17e12)},
                                     'A_11':{'reaction':'2OH = H2O + O','our_value':np.log(1734217478483.0261),'hong_value':np.log(2.355e12)},
                                     'Sigma_1':{'reaction':'sigma H2O2','our_value':-.03846,'hong_value':0},
                                     'Sigma_2':{'reaction':'sigma_HO2','our_value':.0721,'hong_value':0}}
            
            error_dictonary =  {'A_5':{'reaction':'H2O2 + M = 2OH + M','our_value':None,'hong_value':0},
                                     'A_6':{'reaction':'OH + H2O2 = H2O + HO2','our_value':np.log(5624842396127.52),'hong_value':0},
                                     'A_7':{'reaction': 'OH + HO2 = H2O + O2' , 'our_value':np.log(16646221572429.6),'hong_value':0},
                                     'A_8':{'reaction':'2HO2 = H2O2 + O2','our_value':np.log(806831822530.157),'hong_value':0},
                                     'A_11':{'reaction':'2OH = H2O + O','our_value':np.log(1730749579423.63),'hong_value':0},
                                     'Sigma_1':{'reaction':'sigma H2O2','our_value':-.03846,'hong_value':0},
                                     'Sigma_2':{'reaction':'sigma_HO2','our_value':.0721,'hong_value':0}}
            Z = rv.pdf(pos)
            plt.figure()
            levels = [.65,.95,.99]
            #contour = plt.contour(X, Y, Z, levels, colors='k')
            #plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
           # plt.colorbar(contour_filled)
            plt.contour(X,Y,Z)
            plt.xlabel(couple[0])
            plt.ylabel(couple[1])

        
#            
#            plt.figure()
#            
#            Z_test = mlab.bivariate_normal(X, Y,np.sqrt(covariance_couple),np.sqrt(covariance_couple),mu_x,mu_y)
#            z1 = mlab.bivariate_normal(0, 1 * np.sqrt(covariance_couple), np.sqrt(covariance_couple), np.sqrt(covariance_couple),mu_x,mu_y)
#            z2 = mlab.bivariate_normal(0, 2 * np.sqrt(covariance_couple), np.sqrt(covariance_couple), np.sqrt(covariance_couple),mu_x,mu_y)
#            z3 = mlab.bivariate_normal(0, 3 * np.sqrt(covariance_couple), np.sqrt(covariance_couple), np.sqrt(covariance_couple),mu_x,mu_y)
#            
##plot Gaussian:
#            im = plt.imshow(Z_test,interpolation='bilinear', origin='lower',
#                 extent=(-50,50,-50,50),cmap=cm.gray)
##Plot contours at whatever z values we want:
#            CS = plt.contour(Z_test, [z1, z2, z3], origin='lower', extent=(-50,50,-50,50),colors='red')            
            
            
            
            
            
            if bool(additional_dictionary):
                plt.xlabel(additional_dictionary[couple[0]]['reaction'])
                plt.ylabel(additional_dictionary[couple[1]]['reaction'])
                x_error = (additional_dictionary[couple[0]]['hong_value'])*(error_dictonary[couple[0]]['hong_value'])
                print(x_error,'this is the x error')

                y_error = (additional_dictionary[couple[1]]['hong_value'])*(error_dictonary[couple[1]]['hong_value'])
                print(y_error,'this is the y error')
                plt.errorbar(additional_dictionary[couple[0]]['hong_value'],additional_dictionary[couple[1]]['hong_value'],xerr=x_error,yerr=y_error)
                
                plt.scatter(additional_dictionary[couple[0]]['hong_value'],additional_dictionary[couple[1]]['hong_value'],zorder=4,label='Hong Values From Table')
                
                plt.scatter(additional_dictionary[couple[0]]['our_value'],additional_dictionary[couple[1]]['our_value'],zorder=4,marker='x',label='MSI Values')
                plt.legend()

            if bool(joint_data_csv):
                df2 = pd.read_csv(joint_data_csv)
                #plt.figure()
                plt.scatter(df2[couple[0]], df2[couple[1]])
                
                plt.savefig(self.working_directory+'/'+couple[0]+'_'+couple[1]+'_distribution'+'_.pdf',bbox_inches='tight')


    def difference_plotter(self,
                           paramter_list,
                           optimized_cti_file='',
                           pdf_distribution_file=''):
                        
            
        all_parameters = self.shock_tube_instance.posterior_diag_df['parameter'].tolist()
        df = self.shock_tube_instance.posterior_diag_df
        gas_optimized = ct.Solution(optimized_cti_file)
        
        for parameter in paramter_list:
            indx = all_parameters.index(parameter)
            variance = df['value'][indx]
            letter,number = parameter.split('_')
            number = int(number)
            A=gas_optimized.reaction(number).rate.pre_exponential_factor
            n=gas_optimized.reaction(number).rate.temperature_exponent
            Ea=gas_optimized.reaction(number).rate.activation_energy
            
            if letter =='A':
                mu = np.log(A*1000)
                sigma = math.sqrt(variance)
                sigma = sigma
                
            if letter == 'n':
                mu = n
                sigma = math.sqrt(variance)
                #sigma = sigma/2
            if letter == 'Ea':
                mu=Ea/1000/4.184            
                sigma = math.sqrt(variance)
                sigma = sigma*ct.gas_constant/(1000*4.184)
                #sigma = sigma/2

            
            
            x = np.linspace(mu - 6*sigma, mu + 6*sigma, 100)
            #plt.figure()
            #plt.plot(x, stats.norm.pdf(x, mu, sigma))
           # plt.xlabel(parameter)
           # plt.ylabel('pdf')
           # plt.savefig(self.working_directory+'/'+parameter+'_distribution'+'_.pdf',bbox_inches='tight')

            if bool(pdf_distribution_file):
                df2 = pd.read_csv(pdf_distribution_file)
                #temp = np.log(np.exp(df2[parameter].values)/9.33e13)
                #plt.plot(temp,df2['pdf_'+parameter])
                interp_y = np.interp(df2[parameter],x,stats.norm.pdf(x, mu, sigma))
                plt.figure()
                plt.plot(df2[parameter],interp_y)
                plt.plot(df2[parameter],df2['pdf_'+parameter])
                interp_x = np.interp(df2['pdf_'+parameter],stats.norm.pdf(x,mu,sigma),x)
                y_shift = np.divide((df2['pdf_'+parameter] - interp_y),df2['pdf_'+parameter])
                x_shift = np.divide((df2[parameter] - interp_x),df2[parameter])
                plt.figure()
                plt.title('Percent Difference In Y')
                plt.plot(y_shift)
                plt.xlabel(parameter)
                plt.figure()
                plt.plot(x_shift)
                plt.title('Percent Difference In X')
                plt.xlabel(parameter)
              
     def plotting_histograms_of_MSI_simulations(self,experiments_want_to_plot_data_from=[],bins='auto',directory_to_save_images=''):
        s_shape = self.S_matrix.shape[1]
        if self.k_target_value_S_matrix.any():
            target_values_for_s = self.k_target_value_S_matrix
            s_shape = s_shape+target_values_for_s.shape[0]
        y_shape = self.y_matrix.shape[0]
        difference = y_shape-s_shape
        y_values = self.y_matrix[0:difference,0]
        Y_values = self.Y_matrix[0:difference,0]
        self.lengths_of_experimental_data()

        #plotting_Y Histagrams 
        if bool(experiments_want_to_plot_data_from):
            y_values = []
            Y_values = []
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        Y_values.append(temp)
                        temp2 = self.y_matrix[start:stop,:]
                        y_values.append(temp2)
                    
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]              
          
                    
                    
                    
                    
            Y_values = np.vstack((Y_values))
            y_values = np.vstack((y_values))
            plt.figure()            
            plt.subplot(2,2,1)

            n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
            min_value = min(Y_values)
            max_value=max(Y_values)
            plt.xlim([min_value,max_value])
            plt.xlabel('Y')
            plt.suptitle('Including Experiments_'+ str(experiments_want_to_plot_data_from), fontsize=10)

            plt.subplot(2,2,2)
            plt.hist(y_values,bins=bins,align='mid')
            plt.xlabel('y')

            plt.subplot(2,2,3)
            plt.hist(Y_values,bins=bins,density=True,align='mid')
            plt.xlabel('Y')
            plt.ylabel('normalized')

            plt.subplot(2,2,4)
            plt.hist(y_values,bins=bins,density=True,align='mid')
            plt.xlabel('y')      
            plt.ylabel('normalized')
            
            
            
            
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
            plt.savefig(directory_to_save_images+'/'+'Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_4.pdf',dpi=1000,bbox_inches='tight')
            
            
            
            #plotting two fold plots 
            plt.figure()            
            plt.subplot(2,1,1)
            plt.title('Including Experiments_'+ str(experiments_want_to_plot_data_from))

            n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
            plt.xlabel('Y')
            #plt.xlim([-1,1])

            plt.subplot(2,1,2)
            plt.hist(y_values,bins=bins,align='mid')
            plt.xlabel('y')
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
            plt.savefig(directory_to_save_images+'/'+'Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2.pdf',dpi=1000,bbox_inches='tight')

#plotting normalized values
            plt.figure()            
            plt.subplot(2,1,1)
            n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid',density=True)
            plt.xlabel('Y')
            plt.title('Including Experiments_'+ str(experiments_want_to_plot_data_from))
            plt.ylabel('normalized')

            plt.subplot(2,1,2)
            plt.hist(y_values,bins=bins,align='mid',density=True)
            plt.xlabel('y')
            plt.ylabel('normalized')
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
            plt.savefig(directory_to_save_images+'/'+'Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2_normalized.pdf',dpi=1000,bbox_inches='tight')

            



        else:
            plt.figure()            
            plt.subplot(2,2,1)
            min_value = min(Y_values)
            max_value=max(Y_values)
            plt.xlim([min_value,max_value])
            n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
            #plt.xlim([min_value,max_value])
            plt.xlabel('Y')
            plt.suptitle("Including All Experiments", fontsize=10)

            plt.subplot(2,2,2)
            plt.hist(y_values,bins=bins,align='mid')
            plt.xlabel('y')

            plt.subplot(2,2,3)
            plt.hist(Y_values,bins=bins,density=True,align='mid')
            plt.xlabel('Y')
            plt.ylabel('normalized')

            plt.subplot(2,2,4)
            plt.hist(y_values,bins=bins,density=True,align='mid')
            plt.xlabel('y')      
            plt.ylabel('normalized')
            
            
            
            
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
            plt.savefig(directory_to_save_images+'/'+'Including all Experiments'+'_Yy_hist_4.pdf',dpi=1000,bbox_inches='tight')
            
            
            
            #plotting two fold plots 
            plt.figure()            
            plt.subplot(2,1,1)
            min_value = np.min(Y_values)
            max_value = np.max(Y_values)
            plt.title('Including all Experiments')

            n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
            plt.xlabel('Y')
            #plt.xlim([-1,1])

            plt.subplot(2,1,2)
            plt.hist(y_values,bins=bins,align='mid')
            plt.xlabel('y')
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
            plt.savefig(directory_to_save_images+'/'+'Including all Experiments'+'_Yy_hist_2.pdf',dpi=1000,bbox_inches='tight')
            
#plotting normalized values
            plt.figure()            
            plt.subplot(2,1,1)
            min_value = np.min(Y_values)
            max_value = np.max(Y_values)
            n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid',density=True)
            plt.xlabel('Y')
            plt.title('Including all Experiments')
            plt.ylabel('normalized')

            plt.subplot(2,1,2)
            plt.hist(y_values,bins=bins,align='mid',density=True)
            plt.xlabel('y')
            plt.ylabel('normalized')
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
            plt.savefig(directory_to_save_images+'/'+'Including all Experiments'+'_Yy_hist_2_normalized.pdf',dpi=1000,bbox_inches='tight')

    def plotting_T_and_time_full_simulation(self,experiments_want_to_plot_data_from=[],directory_to_save_images=''):
        init_temperature_list = []
        for exp in self.exp_dict_list_original:
            init_temperature_list.append(exp['simulation'].temperature)
        tottal_times = []
        temperature_list_full_simulation = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            single_exp_dict = []
            temp_list_single_experiment = []
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_exp_dict.append(exp['experimental_data'][observable_counter]['Time']*1e3)
                    interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                    temp_list_single_experiment.append(interploated_temp)

                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    single_exp_dict.append(exp['experimental_data'][observable_counter]['Time']*1e3)
                    interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                    temp_list_single_experiment.append(interploated_temp) 
                    #print(interploated_temp.shape ,exp['experimental_data'][observable_counter]['Time'].shape )


                    observable_counter+=1
                    
            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_exp_dict.append(exp['absorbance_experimental_data'][k]['time']*1e3)
                    interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                    temp_list_single_experiment.append(interploated_temp)
                    #print(interploated_temp.shape, exp['absorbance_experimental_data'][k]['time'].shape )
  

            tottal_times.append(single_exp_dict)
            temperature_list_full_simulation.append(temp_list_single_experiment)
            
            
        if bool(experiments_want_to_plot_data_from)==False:
            experiments_want_to_plot_data_from = np.arange(0,len(self.exp_dict_list_optimized))
        else:
            experiments_want_to_plot_data_from = experiments_want_to_plot_data_from
           
        y_values = []
        Y_values = []
        temperature_values_list = []
        time_values_list = []
        full_temperature_range_list = []
        start = 0
        stop = 0             
        for x in range(len(self.simulation_lengths_of_experimental_data)):
            single_experiment_Y =[]
            single_experiment_y =[]
            single_experiment_temperature_values_list=[]
            single_experiment_time_values_list=[]
            single_experiment_full_temp_range=[]
            for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                stop = self.simulation_lengths_of_experimental_data[x][y] + start
                if x in experiments_want_to_plot_data_from:
                    temp = self.Y_matrix[start:stop,:]
                    single_experiment_Y.append(temp)
                    temp2 = self.y_matrix[start:stop,:]
                    single_experiment_y.append(temp2)
                    intial_temp = np.array(([init_temperature_list[x]]*temp.shape[0]))
                    intial_temp = intial_temp.reshape((intial_temp.shape[0],1))
                    single_experiment_temperature_values_list.append(intial_temp)
                    
                    
                    time_values = tottal_times[x][y].values
                    time_values = time_values.reshape((time_values.shape[0],1))
                    single_experiment_time_values_list.append(time_values)
                    
                    temperature_full = temperature_list_full_simulation[x][y]
                    temperature_full = temperature_full.reshape((temperature_full.shape[0],1))
                    single_experiment_full_temp_range.append(temperature_full)

                    start = start + self.simulation_lengths_of_experimental_data[x][y]
                else:
                    start = start + self.simulation_lengths_of_experimental_data[x][y] 
            Y_values.append(single_experiment_Y)
            y_values.append(single_experiment_y)
            temperature_values_list.append(single_experiment_temperature_values_list)
            time_values_list.append(single_experiment_time_values_list)
            full_temperature_range_list.append(single_experiment_full_temp_range)
        
        x = np.arange(10)
        ys = [i+x+(i*x)**2 for i in range(10)]
        colors=cm.rainbow(np.linspace(0,1,30))

        #colors = cm.rainbow(np.linspace(0, 1, len(ys)))

        plt.figure() 
        for x,simulation_list in enumerate(Y_values):
            for y,lst in enumerate(Y_values[x]):
                plt.subplot(2,1,1)
                plt.xlabel('Y')
                plt.ylabel('Time')
                plt.scatter(Y_values[x][y],time_values_list[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color=colors[x])
                plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                plt.subplot(2,1,2)
                plt.scatter(y_values[x][y],time_values_list[x][y],color=colors[x])
                plt.xlabel('y')
                plt.ylabel('Time')
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                plt.savefig(directory_to_save_images+'/'+'Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_vs_time.pdf',dpi=1000,bbox_inches='tight')

    
                
                
        plt.figure() 
    
        for x,simulation_list in enumerate(Y_values):
            for y,lst in enumerate(Y_values[x]):
                plt.subplot(2,1,1)
                plt.scatter(Y_values[x][y],temperature_values_list[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color=colors[x])
                plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                plt.xlabel('Y')
                plt.ylabel('Initial Simulation Temp')
                plt.subplot(2,1,2)
                plt.scatter(y_values[x][y],temperature_values_list[x][y],color=colors[x])    
                plt.xlabel('y')
                plt.ylabel('Initial Simulation Temp')
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                plt.savefig(directory_to_save_images+'/'+'Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_vs_init_temp.pdf',dpi=1000,bbox_inches='tight')

        plt.figure() 
        for x,simulation_list in enumerate(Y_values):
            for y,lst in enumerate(Y_values[x]):
                plt.subplot(2,1,1)
                plt.scatter(Y_values[x][y],full_temperature_range_list[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color=colors[x])
                plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))

                plt.xlabel('Y')
                plt.ylabel('Temperature')
                plt.subplot(2,1,2)
                plt.scatter(y_values[x][y],full_temperature_range_list[x][y],color=colors[x])      
                plt.xlabel('y')
                plt.ylabel('Temperature')
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                plt.savefig(directory_to_save_images+'/'+'Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_vs_temperature.pdf',dpi=1000,bbox_inches='tight')

        return 

#working here 
    def plotting_histograms_of_individual_observables(self,experiments_want_to_plot_data_from,bins='auto',directory_to_save_images='',csv=''):
        s_shape = self.S_matrix.shape[1]
        if self.k_target_value_S_matrix.any():
            target_values_for_s = self.k_target_value_S_matrix
            s_shape = s_shape+target_values_for_s.shape[0]
        y_shape = self.y_matrix.shape[0]
        difference = y_shape-s_shape
        y_values = self.y_matrix[0:difference,0]
        Y_values = self.Y_matrix[0:difference,0]

        self.lengths_of_experimental_data()

        #plotting_Y Histagrams 
        #obserervable_list = []
        
        observables_tottal = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]

        
        if bool(experiments_want_to_plot_data_from):
            print('inside here')
            y_values = []
            Y_values = []
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)

                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]              
          
                    
                    
                    
                    
            for i,observable in enumerate(empty_nested_observable_list_Y):
                if bool(observable):
                    Y_values = np.vstack((observable))
                    y_values = np.vstack((empty_nested_observable_list_y[i]))
                    
                    plt.figure()            
                    plt.subplot(2,2,1)
        
                    n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
                    min_value = min(Y_values)
                    max_value=max(Y_values)
                    plt.xlim([min_value,max_value])
                    plt.xlabel('Y')
                    plt.suptitle(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from), fontsize=10)
        
                    plt.subplot(2,2,2)
                    plt.hist(y_values,bins=bins,align='mid')
                    plt.xlabel('y')
        
                    plt.subplot(2,2,3)
                    plt.hist(Y_values,bins=bins,density=True,align='mid')
                    plt.xlabel('Y')
                    plt.ylabel('normalized')
        
                    plt.subplot(2,2,4)
                    plt.hist(y_values,bins=bins,density=True,align='mid')
                    plt.xlabel('y')      
                    plt.ylabel('normalized')
            
            
            
                    
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                    #plt.savefig(directory_to_save_images+'/'+str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_4.pdf',dpi=1000,bbox_inches='tight')
                    
                    
                    
                    #plotting two fold plots 
                    plt.figure()            
                    plt.subplot(2,1,1)
                    plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
        
                    n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
                    plt.xlabel('Y')
                    #plt.xlim([-1,1])
        
                    plt.subplot(2,1,2)
                    plt.hist(y_values,bins=bins,align='mid')
                    plt.xlabel('y')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                    #plt.savefig(directory_to_save_images+'/'+str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2.pdf',dpi=1000,bbox_inches='tight')
        
        #plotting normalized values
                    plt.figure()            
                    plt.subplot(2,1,1)
                    n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid',density=True)
                    plt.xlabel('Y')
                    plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    plt.ylabel('normalized')
        
                    plt.subplot(2,1,2)
                    plt.hist(y_values,bins=bins,align='mid',density=True)
                    plt.xlabel('y')
                    plt.ylabel('normalized')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                    

                    #plotting two fold plots 
                    plt.figure()            
                    plt.subplot(2,1,1)
                    plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
        
                    n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid')
                    plt.xlabel('Y')
                    #plt.xlim([-1,1])
        
                    plt.subplot(2,1,2)
                    plt.hist(y_values,bins=bins,align='mid')
                    plt.xlabel('y')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                    #plt.savefig(directory_to_save_images+'/'+str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2.pdf',dpi=1000,bbox_inches='tight')
        
        #plotting normalized values
                    plt.figure()            
                    plt.subplot(2,1,1)
                    n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid',density=True)
                    plt.xlabel('Y')
                    plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    plt.ylabel('normalized')
        
                    plt.subplot(2,1,2)
                    plt.hist(y_values,bins=bins,align='mid',density=True)
                    plt.xlabel('y')
                    plt.ylabel('normalized')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                   # plt.savefig(directory_to_save_images+'/'+str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2_normalized.pdf',dpi=1000,bbox_inches='tight')

    def plotting_histograms_of_individual_observables_for_paper_2(self,experiments_want_to_plot_data_from,experiments_want_to_plot_data_from_2=[],bins='auto',directory_to_save_images='',csv=''):
        s_shape = self.S_matrix.shape[1]
        if self.k_target_value_S_matrix.any():
            target_values_for_s = self.k_target_value_S_matrix
            s_shape = s_shape+target_values_for_s.shape[0]
        y_shape = self.y_matrix.shape[0]
        difference = y_shape-s_shape
        y_values = self.y_matrix[0:difference,0]
        Y_values = self.Y_matrix[0:difference,0]
        self.lengths_of_experimental_data()

        #plotting_Y Histagrams 
        #edit this part 
        #obserervable_list = []
        
        observables_tottal = []
        
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        
        empty_nested_observable_list_Z = [[] for x in range(len(observables_unique))]
        
        empty_nested_observable_list_Y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_Z_2 = [[] for x in range(len(observables_unique))]

        if bool(experiments_want_to_plot_data_from):
            print('inside here')
            y_values = []
            Y_values = []
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)

                        temp3 = self.z_matrix[start:stop,:]
                        empty_nested_observable_list_Z[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]              
          
        if bool(experiments_want_to_plot_data_from_2):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from_2:
                        print(x)
                        print(current_observable,'this is current')

                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y_2[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y_2[observables_unique.index(current_observable)].append(temp2)

                        temp3 = self.z_matrix[start:stop,:]
                        empty_nested_observable_list_Z_2[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]                     
                    
        import matplotlib.gridspec as gridspec
    
        fig = plt.figure(figsize=(6,7))

        gs = gridspec.GridSpec(3, 1,height_ratios=[3,3,3],wspace=0.1,hspace=0.1)
        gs.update(wspace=0, hspace=0.7)
        ax1=plt.subplot(gs[0])
        ax2=plt.subplot(gs[1])
        ax3=plt.subplot(gs[2])  
        for i,observable in enumerate(empty_nested_observable_list_Y):
            new_Y_test_2 =[]
            if bool(observable):
                Y_values = np.vstack((observable))
                y_values = np.vstack((empty_nested_observable_list_y[i]))
                z_values = np.vstack((empty_nested_observable_list_Z[i]))
                indecies = np.argwhere(z_values > 100)
                new_y_test = copy.deepcopy(Y_values)
                new_y_test = np.delete(new_y_test,indecies)
             #   print(indecies.shape)
            #    print(indecies)
            #    print(i)
                if bool(experiments_want_to_plot_data_from_2) and bool(empty_nested_observable_list_y_2[i]):
                    Y_values_2 = np.vstack((empty_nested_observable_list_Y_2[i]))
                    y_values_2 = np.vstack((empty_nested_observable_list_y_2[i]))
                    z_values_2 = np.vstack((empty_nested_observable_list_Z_2[i]))
                    indecies_2 = np.argwhere(z_values_2 > 100)
                    new_Y_test_2 = copy.deepcopy(Y_values_2)
                    new_Y_test_2 = np.delete(new_Y_test_2,indecies_2)
                    
            
                #plt.figure()            
                #plt.subplot(1,1,1)
                #plt.subplots(3,1,1)
                #n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid',density=True,label='Hong Experiments')
                test = [-0.06402874, -0.05325865, -0.04248857, -0.03171848, -0.02094839, -0.0101783,
                        0.00059179,  0.01136188,  0.02213197,  0.03290205,  0.04367214,  0.05444223,
                        0.06521232,  0.07598241,  0.0867525,   0.09752259,  0.10829268]
                if i ==0:
                #n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    #ax1.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n,bins_test_1,patches = ax1.hist(new_y_test,bins=bins ,align='mid',density=True,label='#1')
                    ax1.set_xlim(left=-.3, right=.3, emit=True, auto=False)
                    ax1.set_ylim(top=15,bottom=0)

                    ax1.set_xlabel('Y')
                    ax1.set_xlabel('Relative Difference')
                #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    ax1.set_title(str(observables_unique[i]))
                    ax1.set_ylabel('pdf')

                #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                   # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        #ax1.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        ax1.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='#2')

                    if bool(csv):
                        df = pd.read_csv(csv)
                        #ax1.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')
                        #ax1.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='#3')

                    ax1.legend()
                if i ==1:
                #n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n,bins_test_2,patches = ax2.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    ax2.set_xlim(left=-.08, right=.08, emit=True, auto=False)
                    ax2.set_ylim(top=28,bottom=0)
                    ax2.set_xlabel('Y')
                    ax2.set_xlabel('Relative Difference')
                #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    #ax2.set_title(str(observables_unique[i]))
                    ax2.set_title(r'H$_2$O')
                    ax2.set_ylabel('pdf')
                    
                #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                   # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        ax2.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        
                    if bool(csv):
                        df = pd.read_csv(csv)
                        #ax2.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')
                
                
                if i ==3:
                #n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n,bins_test_3,patches = ax3.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    ax3.set_xlim(left=-.15, right=.15, emit=True, auto=False)
                    ax3.set_ylim(top=12,bottom=0)

                    ax3.set_xlabel('Y')
                    ax3.set_xlabel('Relative Difference')
                    ax3.set_ylabel('pdf')
                #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    ax3.set_title(str(observables_unique[i]))
                    ax3.set_title('Absorbance '+ str(observables_unique[i])+ ' nm')

                #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                        print('inside here')
                        print(experiments_want_to_plot_data_from_2)
                   # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        ax3.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        
                    if bool(csv):
                        df = pd.read_csv(csv)
                        #ax3.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')                        
                
                
                

                    plt.savefig(directory_to_save_images+'/'+str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2_normalized.pdf',dpi=1000,bbox_inches='tight')    
    
    
    
    def plotting_histograms_of_individual_observables_for_paper(self,experiments_want_to_plot_data_from,experiments_want_to_plot_data_from_2=[],bins='auto',directory_to_save_images='',csv=''):
        s_shape = self.S_matrix.shape[1]
        if self.k_target_value_S_matrix.any():
            target_values_for_s = self.k_target_value_S_matrix
            s_shape = s_shape+target_values_for_s.shape[0]
        y_shape = self.y_matrix.shape[0]
        difference = y_shape-s_shape
        y_values = self.y_matrix[0:difference,0]
        Y_values = self.Y_matrix[0:difference,0]
        self.lengths_of_experimental_data()

        #plotting_Y Histagrams 
        #obserervable_list = []
        
        observables_tottal = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        
        empty_nested_observable_list_Z = [[] for x in range(len(observables_unique))]
        
        empty_nested_observable_list_Y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_Z_2 = [[] for x in range(len(observables_unique))]

        if bool(experiments_want_to_plot_data_from):
            print('inside here')
            y_values = []
            Y_values = []
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)

                        temp3 = self.z_matrix[start:stop,:]
                        empty_nested_observable_list_Z[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]              
          
        if bool(experiments_want_to_plot_data_from_2):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from_2:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y_2[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y_2[observables_unique.index(current_observable)].append(temp2)

                        temp3 = self.z_matrix[start:stop,:]
                        empty_nested_observable_list_Z_2[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]                     
                    
            import matplotlib.gridspec as gridspec
       
            for i,observable in enumerate(empty_nested_observable_list_Y):
                if bool(observable):
                    Y_values = np.vstack((observable))
                    y_values = np.vstack((empty_nested_observable_list_y[i]))
                    z_values = np.vstack((empty_nested_observable_list_Z[i]))
                    indecies = np.argwhere(z_values > 100)
                    new_y_test = copy.deepcopy(Y_values)
                    new_y_test = np.delete(new_y_test,indecies)
                    print(indecies.shape)
                    print(indecies)
                    if bool(experiments_want_to_plot_data_from_2) and bool(empty_nested_observable_list_y_2[i]):
                        
                        Y_values_2 = np.vstack((empty_nested_observable_list_Y_2[i]))
                        y_values_2 = np.vstack((empty_nested_observable_list_y_2[i]))
                        z_values_2 = np.vstack((empty_nested_observable_list_Z_2[i]))
                        indecies_2 = np.argwhere(z_values_2 > 100)
                        new_Y_test_2 = copy.deepcopy(Y_values_2)
                        new_Y_test_2 = np.delete(new_Y_test_2,indecies_2)
                
                    plt.figure()            
                    plt.subplot(1,1,1)
                    #plt.subplots(3,1,1)
                    #n, bins2, patches = plt.hist(Y_values,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')

                    plt.xlabel('Y')
                    #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    plt.title(str(observables_unique[i]))
                    #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                       # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        plt.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')

                    if bool(csv):
                        df = pd.read_csv(csv)
                        plt.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')
                    plt.legend()
        return 

    def plotting_T_and_time_full_simulation_individual_observables(self,experiments_want_to_plot_data_from,bins='auto',directory_to_save_images=''):
#working_here

        observables_tottal = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)

                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        #print(interploated_temp.shape, exp['absorbance_experimental_data'][k]['time'].shape )
  
        
        x = np.arange(10)
        ys = [i+x+(i*x)**2 for i in range(10)]
        colors=cm.rainbow(np.linspace(0,1,30))

        #colors = cm.rainbow(np.linspace(0, 1, len(ys)))

        
        for x,observable in enumerate(empty_nested_observable_list_Y):
            if bool(observable):
                plt.figure()
                for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        plt.subplot(2,1,1)
                        plt.xlabel('Y')
                        plt.ylabel('Time')
                        plt.scatter(empty_nested_observable_list_Y[x][y],empty_nested_observable_list_time[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color=colors[x])
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        plt.title(observables_unique[x])

                        plt.subplot(2,1,2)
                        plt.scatter(empty_nested_observable_list_y[x][y],empty_nested_observable_list_time[x][y],color=colors[x])
                        plt.xlabel('y')
                        plt.ylabel('Time')
                        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)

                plt.savefig(directory_to_save_images+'/'+str(observables_unique[x])+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_time.pdf',dpi=1000,bbox_inches='tight')

        for x,observable in enumerate(empty_nested_observable_list_Y):
            if bool(observable):
                plt.figure()
                for y,array in enumerate(empty_nested_observable_list_Y[x]):
                    plt.subplot(2,1,1)
                    plt.scatter(empty_nested_observable_list_Y[x][y],empty_nested_observable_list_temperature[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color=colors[x])
                    #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                    plt.xlabel('Y')
                    plt.ylabel('Temperature')
                    plt.title(observables_unique[x])

                    plt.subplot(2,1,2)
                    
                    plt.scatter(empty_nested_observable_list_y[x][y],empty_nested_observable_list_temperature[x][y],color=colors[x])    
                    plt.xlabel('y')
                    plt.ylabel('Temperature')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                plt.savefig(directory_to_save_images+'/'+str(observables_unique[x])+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_temperature.pdf',dpi=1000,bbox_inches='tight')
        for x,observable in enumerate(empty_nested_observable_list_Y):
            if bool(observable):
                plt.figure()
                for y,array in enumerate(empty_nested_observable_list_Y[x]):
                    plt.subplot(2,1,1)
                    plt.scatter(empty_nested_observable_list_Y[x][y],empty_nested_observable_list_initial_temperature[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color=colors[x])
                    #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                    plt.xlabel('Y')
                    plt.ylabel('Initial Temperature')
                    plt.title(observables_unique[x])

                    plt.subplot(2,1,2)
                    
                    plt.scatter(empty_nested_observable_list_y[x][y],empty_nested_observable_list_initial_temperature[x][y],color=colors[x])    
                    plt.xlabel('y')
                    plt.ylabel('Initial Temperature')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)    
                plt.savefig(directory_to_save_images+'/'+str(observables_unique[x])+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_initial_temperature.pdf',dpi=1000,bbox_inches='tight')





    def plotting_T_and_time_full_simulation_individual_observables_for_paper(self,experiments_want_to_plot_data_from,
                                                                             bins='auto',
                                                                             directory_to_save_images='',csv='',experiments_want_to_plot_data_from_2=[]):
#working_here

        observables_tottal = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)

                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
####################################################################################################################################################################################################################    
        empty_nested_observable_list_Y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature_2 = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from_2):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from_2:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y_2[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y_2[observables_unique.index(current_observable)].append(temp2)

                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from_2:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time_2[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time_2[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from_2:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time_2[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])    
###################################################################################################################################################################################################################  
        
        x = np.arange(10)
        ys = [i+x+(i*x)**2 for i in range(10)]
        colors=cm.rainbow(np.linspace(0,1,30))

        #colors = cm.rainbow(np.linspace(0, 1, len(ys)))

        
        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])
            if bool(observable):
                plt.figure()
                if bool(csv):
                    df = pd.read_csv(csv)
                    plt.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=0.5,color='k',zorder=4)      
                for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        plt.subplot(1,1,1)
                        plt.xlabel('Y')
                        plt.ylabel('Time')
                        plt.scatter(empty_nested_observable_list_Y[x][y],empty_nested_observable_list_time[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color='blue')
                        

                        if y<length_of_2nd_list:
                            plt.scatter(empty_nested_observable_list_Y_2[x][y],empty_nested_observable_list_time_2[x][y],color='red',zorder=4)

                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        plt.title(observables_unique[x])
                

#                        plt.subplot(2,1,2)
#                        plt.scatter(empty_nested_observable_list_y[x][y],empty_nested_observable_list_time[x][y],color=colors[x])
#                        plt.xlabel('y')
#                        plt.ylabel('Time')
#                        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)

                #plt.savefig(directory_to_save_images+'/'+str(observables_unique[x])+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_time.pdf',dpi=1000,bbox_inches='tight')

        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])

            if bool(observable):
                plt.figure()
                if bool(csv):
                    df = pd.read_csv(csv)
                    plt.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),alpha=0.5,color='k',zorder=4)                    
                for y,array in enumerate(empty_nested_observable_list_Y[x]):
                    plt.subplot(1,1,1)
                    plt.scatter(empty_nested_observable_list_Y[x][y],empty_nested_observable_list_temperature[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color='blue')
                    #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                    plt.xlabel('Y')
                    plt.ylabel('Temperature')
                    plt.title(observables_unique[x])
                    
                    if y<length_of_2nd_list:
                        plt.scatter(empty_nested_observable_list_Y_2[x][y],empty_nested_observable_list_temperature_2[x][y],color='red',zorder=4)


#                    plt.subplot(2,1,2)
#                    
#                    plt.scatter(empty_nested_observable_list_y[x][y],empty_nested_observable_list_temperature[x][y],color=colors[x])    
#                    plt.xlabel('y')
#                    plt.ylabel('Temperature')
#                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)
                #plt.savefig(directory_to_save_images+'/'+str(observables_unique[x])+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_temperature.pdf',dpi=1000,bbox_inches='tight')
        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])
            
            if bool(observable):
                plt.figure()
                if bool(csv):
                    df = pd.read_csv(csv)
                    plt.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_initial_Temperature'].dropna(),alpha=0.5,color='k',zorder=4)                 
                for y,array in enumerate(empty_nested_observable_list_Y[x]):
                    plt.subplot(1,1,1)
                    plt.scatter(empty_nested_observable_list_Y[x][y],empty_nested_observable_list_initial_temperature[x][y],label='Experiment_'+str(x)+'_observable_'+str(y),color='blue')
                    #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                    plt.xlabel('Y')
                    plt.ylabel('Initial Temperature')
                    plt.title(observables_unique[x])
                    if y<length_of_2nd_list:
                        plt.scatter(empty_nested_observable_list_Y_2[x][y],empty_nested_observable_list_initial_temperature_2[x][y],color='red',zorder=4)
#                    plt.subplot(2,1,2)
#                    
#                    plt.scatter(empty_nested_observable_list_y[x][y],empty_nested_observable_list_initial_temperature[x][y],color=colors[x])    
#                    plt.xlabel('y')
#                    plt.ylabel('Initial Temperature')
#                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.5, hspace=.5)    
                #plt.savefig(directory_to_save_images+'/'+str(observables_unique[x])+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_initial_temperature.pdf',dpi=1000,bbox_inches='tight')

    def plotting_T_and_time_full_simulation_individual_observables_for_paper_2(self,experiments_want_to_plot_data_from,
                                                                             bins='auto',
                                                                             directory_to_save_images='',csv='',experiments_want_to_plot_data_from_2=[]):
#working_here

        observables_tottal = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_z = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)
                        
                        temp3 = self.z_matrix[start:stop,:]
                        empty_nested_observable_list_z[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
####################################################################################################################################################################################################################    
        empty_nested_observable_list_Y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_z_2 = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from_2):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from_2:
                        temp = self.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y_2[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = self.y_matrix[start:stop,:]
                        empty_nested_observable_list_y_2[observables_unique.index(current_observable)].append(temp2)
                        
                        temp3 = self.z_matrix[start:stop,:]
                        empty_nested_observable_list_z_2[observables_unique.index(current_observable)].append(temp3)
                        
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from_2:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time_2[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time_2[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from_2:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time_2[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])    
###################################################################################################################################################################################################################  
        
        x = np.arange(10)
        ys = [i+x+(i*x)**2 for i in range(10)]
        colors=cm.rainbow(np.linspace(0,1,30))

        #colors = cm.rainbow(np.linspace(0, 1, len(ys)))
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(6,7))
        gs = gridspec.GridSpec(3, 1,height_ratios=[3,3,3],wspace=0.025,hspace=0.1)
        gs.update(wspace=0, hspace=0.7)
        ax1=plt.subplot(gs[0])
        ax2=plt.subplot(gs[1])
        ax3=plt.subplot(gs[2]) 
        
        fig2 = plt.figure(figsize=(6,7))
        gs2 = gridspec.GridSpec(3, 1,height_ratios=[3,3,3],wspace=0.025,hspace=0.1)
        gs2.update(wspace=0, hspace=0.7)
        ax4=plt.subplot(gs2[0])
        ax5=plt.subplot(gs2[1])
        ax6=plt.subplot(gs2[2]) 


        
        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])
            if bool(observable):

                print(x)
                if x ==0:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax1.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=1,color='g',zorder=4,label='_nolegend_')      
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                            
                        
                            z_values = empty_nested_observable_list_z[x][y]
                            indecies = np.argwhere(z_values > 100)
                            new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                            new_y_test = np.delete(new_y_test,indecies)
                            new_time_test = copy.deepcopy(empty_nested_observable_list_time[x][y])
                            new_time_test = new_time_test.values
                            new_time_test = new_time_test.reshape((new_time_test.shape[0],1))
                            new_time_test  = np.delete(new_time_test,indecies)
                            ax1.scatter(new_y_test,new_time_test, c='#1f77b4',alpha=1)
                            ax1.set_xlabel('Relative Difference')
                            ax1.set_ylabel('Time (ms)')
                            
    
                            if y<length_of_2nd_list:
                                z_values_2 = empty_nested_observable_list_z_2[x][y]
                                indecies_2 = np.argwhere(z_values_2 > 100)
                                new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                                new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                                new_time_test_2 = copy.deepcopy(empty_nested_observable_list_time_2[x][y])
                                new_time_test_2 = new_time_test_2.values
                                new_time_test_2 = new_time_test_2.reshape((new_time_test_2.shape[0],1))
                                new_time_test_2  = np.delete(new_time_test_2,indecies_2)
                                ax1.scatter(new_y_test_2,new_time_test_2,color='orange',zorder=3,alpha=.15)
    
                            ax1.set_title(observables_unique[x])
                            ax1.set_xlim(left=-.25, right=.25, emit=True, auto=False)
                    ax1.scatter([],[],c='#1f77b4',label='#1')
                    ax1.scatter([],[],color='orange',label='#2')                                
                    #ax1.scatter([],[],color='green',label='#3')
                    ax1.legend(frameon=False)




                if x ==1:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax2.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=1,color='g',zorder=4)      
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                            
                        
                            z_values = empty_nested_observable_list_z[x][y]
                            indecies = np.argwhere(z_values > 100)
                            new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                            new_y_test = np.delete(new_y_test,indecies)
                            new_time_test = copy.deepcopy(empty_nested_observable_list_time[x][y])
                            new_time_test = new_time_test.values
                            new_time_test = new_time_test.reshape((new_time_test.shape[0],1))
                            new_time_test  = np.delete(new_time_test,indecies)
                            ax2.scatter(new_y_test,new_time_test, c='#1f77b4',alpha=1)
                            ax2.set_xlabel('Relative Difference')
                            ax2.set_ylabel('Time (ms)')
                            ax2.set_xlim(left=-.09, right=.09, emit=True, auto=False)
                            
    
                            if y<length_of_2nd_list:
                                z_values_2 = empty_nested_observable_list_z_2[x][y]
                                indecies_2 = np.argwhere(z_values_2 > 100)
                                new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                                new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                                new_time_test_2 = copy.deepcopy(empty_nested_observable_list_time_2[x][y])
                                new_time_test_2 = new_time_test_2.values
                                new_time_test_2 = new_time_test_2.reshape((new_time_test_2.shape[0],1))
                                new_time_test_2  = np.delete(new_time_test_2,indecies_2)
                                ax2.scatter(new_y_test_2,new_time_test_2,color='orange',zorder=3,alpha=.15)
    
                            #ax2.set_title(observables_unique[x])
                            ax2.set_title(r'H$_2$O')
                            
                            


                if x ==3:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax3.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=1,color='g',zorder=4,)      
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                            
                        
                            z_values = empty_nested_observable_list_z[x][y]
                            indecies = np.argwhere(z_values > 100)
                            new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                            new_y_test = np.delete(new_y_test,indecies)
                            new_time_test = copy.deepcopy(empty_nested_observable_list_time[x][y])
                            new_time_test = new_time_test.values
                            new_time_test = new_time_test.reshape((new_time_test.shape[0],1))
                            new_time_test  = np.delete(new_time_test,indecies)
                            ax3.scatter(new_y_test,new_time_test,c='#1f77b4',alpha=1)
                            ax3.set_xlabel('Relative Difference')
                            ax3.set_ylabel('Time (ms)')
                            ax3.set_xlim(left=-.3, right=.3, emit=True, auto=False)


                            
    
                            if y<length_of_2nd_list:
                                z_values_2 = empty_nested_observable_list_z_2[x][y]
                                indecies_2 = np.argwhere(z_values_2 > 100)
                                new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                                new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                                new_time_test_2 = copy.deepcopy(empty_nested_observable_list_time_2[x][y])
                                new_time_test_2 = new_time_test_2.values
                                new_time_test_2 = new_time_test_2.reshape((new_time_test_2.shape[0],1))
                                new_time_test_2  = np.delete(new_time_test_2,indecies_2)
                                ax3.scatter(new_y_test_2,new_time_test_2,color='orange',zorder=3,alpha=.15)
    
                            #ax3.set_title(observables_unique[x])
                            ax3.set_title('Absorbance '+ str(observables_unique[x])+str(' nm'))

                fig.savefig(directory_to_save_images+'/'+'Three_pannel_plot_'+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_time.pdf',dpi=1000,bbox_inches='tight')

        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])

            if bool(observable):
                
                if x==0:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax4.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),label='_nolegend_',alpha=1,color='green',zorder=4)                    
                    
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        
                        z_values = empty_nested_observable_list_z[x][y]
                        indecies = np.argwhere(z_values > 100)
                        new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                        new_y_test = np.delete(new_y_test,indecies)
                        new_temperature_test = copy.deepcopy(empty_nested_observable_list_temperature[x][y])
                        new_temperature_test = new_temperature_test.reshape((new_temperature_test.shape[0],1))
                        new_temperature_test  = np.delete(new_temperature_test,indecies)                    
                        
                        ax4.scatter(new_y_test,new_temperature_test,c='#1f77b4',alpha=1,label='_nolegend_')
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        ax4.set_xlabel('Relative Difference')
                        ax4.set_ylabel('Temperature (K)')
                        ax4.set_title(observables_unique[x])
                        ax4.set_xlim(left=-.25, right=.25, emit=True, auto=False)
                        
                        if y<length_of_2nd_list:
                            z_values_2 = empty_nested_observable_list_z_2[x][y]
                            indecies_2 = np.argwhere(z_values_2 > 100)
                            new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                            new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                            new_temperature_test_2 = copy.deepcopy(empty_nested_observable_list_temperature_2[x][y])
                            new_temperature_test_2 = new_temperature_test_2.reshape((new_temperature_test_2.shape[0],1))
                            new_temperature_test_2  = np.delete(new_temperature_test_2,indecies_2)                              
                            ax4.scatter(new_y_test_2,new_temperature_test_2,color='orange',zorder=3,alpha=1,label='_nolegend_')
                        
                    ax4.scatter([],[],c='#1f77b4',label='#1')
                    ax4.scatter([],[],c='orange',label='#2')
                    #ax4.scatter([],[],c='green',label='#3')
                    ax4.legend(frameon=False)
   
                                    
                if x==1:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax5.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),alpha=1,color='green',zorder=4)                    
                    
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        
                        z_values = empty_nested_observable_list_z[x][y]
                        indecies = np.argwhere(z_values > 100)
                        new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                        new_y_test = np.delete(new_y_test,indecies)
                        new_temperature_test = copy.deepcopy(empty_nested_observable_list_temperature[x][y])
                        new_temperature_test = new_temperature_test.reshape((new_temperature_test.shape[0],1))
                        new_temperature_test  = np.delete(new_temperature_test,indecies)                    
                        
                        ax5.scatter(new_y_test,new_temperature_test,c='#1f77b4',alpha=1)
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        ax5.set_xlabel('Relative Difference')
                        ax5.set_ylabel('Temperature (K)')
                        #ax5.set_title(observables_unique[x])
                        ax5.set_title(r'H$_2$O')
                        ax5.set_xlim(left=-.09, right=.09, emit=True, auto=False)

                        
                        if y<length_of_2nd_list:
                            z_values_2 = empty_nested_observable_list_z_2[x][y]
                            indecies_2 = np.argwhere(z_values_2 > 100)
                            new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                            new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                            new_temperature_test_2 = copy.deepcopy(empty_nested_observable_list_temperature_2[x][y])
                            new_temperature_test_2 = new_temperature_test_2.reshape((new_temperature_test_2.shape[0],1))
                            new_temperature_test_2  = np.delete(new_temperature_test_2,indecies_2)                              
                            ax5.scatter(new_y_test_2,new_temperature_test_2,color='orange',zorder=3,alpha=1)

                if x==3:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax6.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),alpha=1,color='green',zorder=4)                    
                    
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        
                        z_values = empty_nested_observable_list_z[x][y]
                        indecies = np.argwhere(z_values > 100)
                        new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                        new_y_test = np.delete(new_y_test,indecies)
                        new_temperature_test = copy.deepcopy(empty_nested_observable_list_temperature[x][y])
                        new_temperature_test = new_temperature_test.reshape((new_temperature_test.shape[0],1))
                        new_temperature_test  = np.delete(new_temperature_test,indecies)                    
                        
                        ax6.scatter(new_y_test,new_temperature_test,label='Experiment_'+str(x)+'_observable_'+str(y),c='#1f77b4',alpha=1)
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        ax6.set_xlabel('Relative Difference')
                        ax6.set_ylabel('Temperature (K)')
                        ax6.set_title('Absorbance '+ str(observables_unique[x])+str(' nm'))
                        ax6.set_xlim(left=-.3, right=.3, emit=True, auto=False)

                        if y<length_of_2nd_list:
                            z_values_2 = empty_nested_observable_list_z_2[x][y]
                            indecies_2 = np.argwhere(z_values_2 > 100)
                            new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                            new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                            new_temperature_test_2 = copy.deepcopy(empty_nested_observable_list_temperature_2[x][y])
                            new_temperature_test_2 = new_temperature_test_2.reshape((new_temperature_test_2.shape[0],1))
                            new_temperature_test_2  = np.delete(new_temperature_test_2,indecies_2)                              
                            ax6.scatter(new_y_test_2,new_temperature_test_2,color='orange',zorder=3,alpha=.5)
  
                fig2.savefig(directory_to_save_images+'/'+'Three_pannel_plot'+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_initial_temperature.pdf',dpi=1000,bbox_inches='tight')

    def residual_sum_of_squares(self):
        overall_list = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            single_exp_dict = []
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    difference = exp['experimental_data'][observable_counter][observable].values - (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    square_of_differences = np.square(difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    difference = exp['experimental_data'][observable_counter][observable+'_ppm'].values - (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    square_of_differences = np.square(difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)

                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    difference = exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values- exp['absorbance_model_data'][wl]
                    square_of_differences = np.square(difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)
            
            overall_list.append(single_exp_dict)
        
        return overall_list
                    
                    
    def sum_of_squares_of_Y(self):
        overall_list = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            single_exp_dict = []
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                plt.figure()
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    mean_calculated_experimental = np.mean(exp['experimental_data'][observable_counter][observable].values)
                    #mean_calculated_predicted = (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    
                    difference = exp['experimental_data'][observable_counter][observable].values - mean_calculated_experimental
                    square_of_differences = np.square(difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    mean_calculated_experimental = np.mean(exp['experimental_data'][observable_counter][observable+'_ppm'].values)
                    difference = exp['experimental_data'][observable_counter][observable+'_ppm'].values - mean_calculated_experimental
                    square_of_differences = np.square(difference)
                    sum_of_squares = sum(square_of_differences)
                    #print(sum_of_squares)
                    single_exp_dict.append(sum_of_squares)

                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    mean_calculated_experimental = np.mean(exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values)
                    difference = exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values- mean_calculated_experimental
                    square_of_differences = np.square(difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)
            
            overall_list.append(single_exp_dict)
        
        return overall_list                    



    
    def calculating_R_squared(self,RSS_list,SYY_list):
        overall_list = []
        for i,lst in enumerate(RSS_list):
            single_exp_dict = []
            for j,value in enumerate(RSS_list[i]):
                single_exp_dict.append(1-(RSS_list[i][j]/SYY_list[i][j]))
            overall_list.append(single_exp_dict)

        
        return overall_list        
    
    
    def weighted_sum_of_squares(self):
        
        def uncertainty_calc(relative_uncertainty,absolute_uncertainty,data,experimental_data):

            if 'W' in list(experimental_data.columns):
                weighting_factor = experimental_data['W'].values
                if 'Relative_Uncertainty' in list(experimental_data.columns):
                    time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                    un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                    total_uncertainty = time_dependent_uncertainty/weighting_factor
                    
                else:
                    length_of_data = data.shape[0]
                    relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty) 
                    un_weighted_uncertainty = copy.deepcopy(relative_uncertainty_array)
                    total_uncertainty = un_weighted_uncertainty/weighting_factor
                
            
            
            elif 'Relative_Uncertainty' in list(experimental_data.columns):  
                
                time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                #do we need to take the natrual log of this?
                time_dependent_uncertainty = np.log(time_dependent_uncertainty+1)
                #do we need to take the natrual log of this?
                length_of_data = data.shape[0]
                un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                total_uncertainty = np.divide(time_dependent_uncertainty,(1/length_of_data**.5) )
#                

               
            else:
                length_of_data = data.shape[0]
                relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty)
                
                if absolute_uncertainty != 0:
                #check if this weighting factor is applied in the correct place 
                #also check if want these values to be the natural log values 
                    absolute_uncertainty_array = np.divide(data,absolute_uncertainty)
                    total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array) + np.square(absolute_uncertainty_array)))
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                     #weighting factor
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )
                
                else:
                    #total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array)))
                    total_uncertainty = relative_uncertainty_array
                    #weighting factor
                    
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )

            #make this return a tuple 
            return total_uncertainty,un_weighted_uncertainty
        
        
        overall_list = []
        overall_list_maximum_deviation = []
        overall_list_weighted_uncertainty = []
        overall_percent_difference = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            single_exp_dict = []
            single_exp_dict_maximum_deviation = []
            single_experiment_weighted_uncertainty = []
            single_experiment_percent_difference = []
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    difference = exp['experimental_data'][observable_counter][observable].values - (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['mole_fraction_relative_uncertainty'][observable_counter],
                        exp['uncertainty']['mole_fraction_absolute_uncertainty'][observable_counter],
                        exp['experimental_data'][observable_counter][observable].values,exp['experimental_data'][observable_counter])
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))
                    
                    divided = np.divide(difference,un_weighted_uncertainty)
                    numerator = exp['experimental_data'][observable_counter][observable].values
                    numerator = numerator.reshape((numerator.shape[0],
                                                                               1))   
                    
                    percent_difference = np.divide(divided,np.divide(numerator,un_weighted_uncertainty))
                    single_experiment_percent_difference.append(np.max(np.absolute(percent_difference)))
                    
                    divided2 = np.divide(difference, total_uncertainty)
                    
                    single_exp_dict_maximum_deviation.append(np.max(np.absolute(divided)))
                    square_of_differences = np.square(divided)
                    square_of_differences2 = np.square(divided2)

                    
                    sum_of_squares = sum(square_of_differences)
                    sum_of_squares2 = sum(square_of_differences2)

                    sqrt_of_difference = np.sqrt(sum_of_squares)                    
                    single_exp_dict.append(sqrt_of_difference)
                    single_experiment_weighted_uncertainty.append(sum_of_squares2)
                    
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    difference = exp['experimental_data'][observable_counter][observable+'_ppm'].values - (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['concentration_relative_uncertainty'][observable_counter],
                         exp['uncertainty']['concentration_absolute_uncertainty'][observable_counter],
                         exp['experimental_data'][observable_counter][observable+'_ppm'].values,exp['experimental_data'][observable_counter])    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))                    
                    divided = np.divide(difference,un_weighted_uncertainty)
                    numerator = exp['experimental_data'][observable_counter][observable+'_ppm'].values
                    numerator = numerator.reshape((numerator.shape[0],
                                                                               1))   
                    
                    percent_difference = np.divide(divided,np.divide(numerator,un_weighted_uncertainty))
                    single_experiment_percent_difference.append(np.max(np.absolute(percent_difference)))                    
                    
                    divided2 = np.divide(difference, total_uncertainty)

                    single_exp_dict_maximum_deviation.append(np.max(np.absolute(divided)))

                    square_of_differences = np.square(divided)
                    square_of_differences2 = np.square(divided2)
                    sum_of_squares2 = sum(square_of_differences2)

                    sum_of_squares = sum(square_of_differences)
                    sqrt_of_difference = np.sqrt(sum_of_squares)                    
                    single_exp_dict.append(sqrt_of_difference)
                    single_experiment_weighted_uncertainty.append(sum_of_squares2)

                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    difference = exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values- exp['absorbance_model_data'][wl]
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['absorbance_relative_uncertainty'][k],
                                                             exp['uncertainty']['absorbance_absolute_uncertainty'][k],
                                                             exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp['absorbance_experimental_data'][k])                    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))                    
                    
                    divided = np.divide(difference,un_weighted_uncertainty)
                    numerator =  exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values
                    numerator = numerator.reshape((numerator.shape[0],
                                                                               1))   
                    
                    percent_difference = np.divide(divided,np.divide(numerator,un_weighted_uncertainty))
                    single_experiment_percent_difference.append(np.max(np.absolute(percent_difference)))                    
                    single_exp_dict_maximum_deviation.append(np.max(np.absolute(divided)))
                    divided2 = np.divide(difference, total_uncertainty)
                    square_of_differences2 = np.square(divided2)
                    sum_of_squares2 = sum(square_of_differences2)
                    
                    square_of_differences = np.square(divided)
                    sum_of_squares = sum(square_of_differences)
                    sqrt_of_difference = np.sqrt(sum_of_squares)                    
                    single_exp_dict.append(sqrt_of_difference)
                    single_experiment_weighted_uncertainty.append(sum_of_squares2)

            overall_list.append(single_exp_dict)
            overall_list_maximum_deviation.append(single_exp_dict_maximum_deviation)
            overall_list_weighted_uncertainty.append(single_experiment_weighted_uncertainty)
            overall_percent_difference.append(single_experiment_percent_difference)
        return overall_list,overall_list_maximum_deviation,overall_list_weighted_uncertainty,overall_percent_difference
    
    def weighted_sum_of_squares_of_Y(self):
        
        def uncertainty_calc(relative_uncertainty,absolute_uncertainty,data,experimental_data):

            if 'W' in list(experimental_data.columns):
                weighting_factor = experimental_data['W'].values
                if 'Relative_Uncertainty' in list(experimental_data.columns):
                    time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                    un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                    total_uncertainty = time_dependent_uncertainty/weighting_factor
                    
                else:
                    length_of_data = data.shape[0]
                    relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty) 
                    un_weighted_uncertainty = copy.deepcopy(relative_uncertainty_array)
                    total_uncertainty = un_weighted_uncertainty/weighting_factor
                
            
            
            elif 'Relative_Uncertainty' in list(experimental_data.columns):  
                
                time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                #do we need to take the natrual log of this?
                time_dependent_uncertainty = np.log(time_dependent_uncertainty+1)
                #do we need to take the natrual log of this?
                length_of_data = data.shape[0]
                un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                total_uncertainty = np.divide(time_dependent_uncertainty,(1/length_of_data**.5) )
#                

               
            else:
                length_of_data = data.shape[0]
                relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty)
                
                if absolute_uncertainty != 0:
                #check if this weighting factor is applied in the correct place 
                #also check if want these values to be the natural log values 
                    absolute_uncertainty_array = np.divide(data,absolute_uncertainty)
                    total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array) + np.square(absolute_uncertainty_array)))
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                     #weighting factor
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )
                
                else:
                    #total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array)))
                    total_uncertainty = relative_uncertainty_array
                    #weighting factor
                    
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )

            #make this return a tuple 
            return total_uncertainty,un_weighted_uncertainty        
        
        overall_list = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            single_exp_dict = []
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    mean_calculated_experimental = np.mean(exp['experimental_data'][observable_counter][observable].values)
                    #mean_calculated_predicted = (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['mole_fraction_relative_uncertainty'][observable_counter],
                        exp['uncertainty']['mole_fraction_absolute_uncertainty'][observable_counter],
                        exp['experimental_data'][observable_counter][observable].values,exp['experimental_data'][observable_counter])                    
                    difference = exp['experimental_data'][observable_counter][observable].values - mean_calculated_experimental
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))                    
                    weighted_difference = np.divide(difference,total_uncertainty)
                    square_of_differences = np.square(weighted_difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    mean_calculated_experimental = np.mean(exp['experimental_data'][observable_counter][observable+'_ppm'].values)
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['concentration_relative_uncertainty'][observable_counter],
                         exp['uncertainty']['concentration_absolute_uncertainty'][observable_counter],
                         exp['experimental_data'][observable_counter][observable+'_ppm'].values,exp['experimental_data'][observable_counter])    
                    
                    difference = exp['experimental_data'][observable_counter][observable+'_ppm'].values - mean_calculated_experimental
                    
                                        
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))                    
                    weighted_difference = np.divide(difference,total_uncertainty)

                    square_of_differences = np.square(weighted_difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)

                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    mean_calculated_experimental = np.mean(exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values)
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['absorbance_relative_uncertainty'][k],
                                                             exp['uncertainty']['absorbance_absolute_uncertainty'][k],
                                                             exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp['absorbance_experimental_data'][k])                      
                    difference = exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values- mean_calculated_experimental
                    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))  
                    weighted_difference = np.divide(difference,total_uncertainty)
                    square_of_differences = np.square(weighted_difference)
                    sum_of_squares = sum(square_of_differences)
                    single_exp_dict.append(sum_of_squares)
            
            overall_list.append(single_exp_dict)
        
        return overall_list    


    def plotting_individual_histograms(self,experimental_dict_list,parsed_yaml_list,directory_to_save_images=''):
        
        def uncertainty_calc(relative_uncertainty,absolute_uncertainty,data,experimental_data):

            if 'W' in list(experimental_data.columns):
                weighting_factor = experimental_data['W'].values
                if 'Relative_Uncertainty' in list(experimental_data.columns):
                    time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                    un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                    total_uncertainty = time_dependent_uncertainty/weighting_factor
                    
                else:
                    length_of_data = data.shape[0]
                    relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty) 
                    un_weighted_uncertainty = copy.deepcopy(relative_uncertainty_array)
                    total_uncertainty = un_weighted_uncertainty/weighting_factor
                
            
            
            elif 'Relative_Uncertainty' in list(experimental_data.columns):  
                
                time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                #do we need to take the natrual log of this?
                time_dependent_uncertainty = np.log(time_dependent_uncertainty+1)
                #do we need to take the natrual log of this?
                length_of_data = data.shape[0]
                un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                total_uncertainty = np.divide(time_dependent_uncertainty,(1/length_of_data**.5) )
#                

               
            else:
                length_of_data = data.shape[0]
                relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty)
                
                if absolute_uncertainty != 0:
                #check if this weighting factor is applied in the correct place 
                #also check if want these values to be the natural log values 
                    absolute_uncertainty_array = np.divide(data,absolute_uncertainty)
                    total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array) + np.square(absolute_uncertainty_array)))
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                     #weighting factor
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )
                
                else:
                    #total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array)))
                    total_uncertainty = relative_uncertainty_array
                    #weighting factor
                    
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )

            #make this return a tuple 
            return total_uncertainty,un_weighted_uncertainty
        
        
        for i,exp in enumerate(experimental_dict_list):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if observable in exp['mole_fraction_observables']:
                    difference = exp['experimental_data'][observable_counter][observable].values - (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    log_difference = np.log(exp['experimental_data'][observable_counter][observable].values)-np.log((exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6)
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['mole_fraction_relative_uncertainty'][observable_counter],
                        exp['uncertainty']['mole_fraction_absolute_uncertainty'][observable_counter],
                        exp['experimental_data'][observable_counter][observable].values,exp['experimental_data'][observable_counter])
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))
                    
                    divided = np.divide(difference,un_weighted_uncertainty)
                    numerator = exp['experimental_data'][observable_counter][observable].values
                    numerator = numerator.reshape((numerator.shape[0],
                                                                               1))   
                    
                    log_difference = log_difference.reshape((log_difference.shape[0],
                                                                               1))                       
                    
                    weighted_log_difference = np.divide(log_difference,total_uncertainty)
                    
                    plt.subplot(2,1,1)
                    plt.hist(difference)
                    plt.xlabel('Y')
                    plt.title(observable)
                    
                    plt.subplot(2,1,2)
                    plt.xlabel('y')
                    plt.hist(weighted_log_difference)
                    plt.title(observable)
                    
                    
                    plt.subplot(2,1,1)
                    plt.figure()
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,difference)
                    plt.ylabel('Y')
                    plt.xlabel('Time (ms)')
                    plt.title(observable)
                    plt.subplot(2,1,2)
                    plt.figure(exp['experimental_data'][observable_counter]['Time']*1e3,weighted_log_difference)
                    plt.xlabel('Time (ms)')
                    plt.ylabel('y')
                    
                    
                    plt.subplot(2,1,1)
                    plt.figure(exp['simulation'].timeHistoryInterpToExperiment['temperature'].dropna().values, difference)
                    plt.ylabel('Y')
                    plt.xlabel('Temperature')                    
                    plt.subplot(2,1,2)
                    plt.figure(exp['simulation'].timeHistoryInterpToExperiment['temperature'].dropna().values, weighted_log_difference)
                    plt.ylabel('y')
                    plt.xlabel('Temperature')                         
                    
                    
                    
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    difference = exp['experimental_data'][observable_counter][observable+'_ppm'].values - (exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6
                    log_difference = np.log(exp['experimental_data'][observable_counter][observable+'_ppm'].values) - np.log((exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6)
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['concentration_relative_uncertainty'][observable_counter],
                         exp['uncertainty']['concentration_absolute_uncertainty'][observable_counter],
                         exp['experimental_data'][observable_counter][observable+'_ppm'].values,exp['experimental_data'][observable_counter])    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))                    
                    divided = np.divide(difference,un_weighted_uncertainty)
                    numerator = exp['experimental_data'][observable_counter][observable+'_ppm'].values
                    numerator = numerator.reshape((numerator.shape[0],
                                                                               1))   
                    
                    log_difference = log_difference.reshape((log_difference.shape[0],
                                                                               1))                       
                    
                    weighted_log_difference = np.divide(log_difference,total_uncertainty)
                    
                    plt.figure()
                    plt.subplot(2,1,1)
                    plt.hist(difference)
                    plt.xlabel('Y')
                    plt.title(observable + ' ' + 'Experiment Number' + ' '+str(i))
                    
                    plt.subplot(2,1,2)
                    plt.xlabel('y')
                    plt.hist(weighted_log_difference)
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.5)

                    plt.savefig(directory_to_save_images+'/'+observable+'_'+'Experiment Number_' +str(i)+'_hist.pdf',dpi=1000,bbox_inches='tight')
                    plt.savefig(directory_to_save_images+'/'+observable+'_'+'Experiment Number_' +str(i)+'_hist.png',dpi=1000,bbox_inches='tight')
                    
                    #plt.subplots_adjust( hspace=.2)

                    plt.figure()
                    plt.subplot(2,1,1)
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,difference)
                    plt.ylabel('Y')
                    plt.title(observable+ ' ' + 'Experiment Number' + ' '+str(i))
                    plt.subplot(2,1,2)
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,weighted_log_difference)
                    plt.xlabel('Time (ms)')
                    plt.ylabel('y')
                    plt.savefig(directory_to_save_images+'/'+observable+'_'+'Experiment Number_' +str(i)+'_timeVSy.pdf',dpi=1000,bbox_inches='tight')
                    plt.savefig(directory_to_save_images+'/'+observable+'_'+'Experiment Number_'+str(i)+'_timeVSy.png',dpi=1000,bbox_inches='tight')
                    
                    
                    plt.figure()
                    plt.subplot(2,1,1)
                    interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                    plt.plot(interploated_temp, difference)
                    plt.ylabel('Y')
                    plt.title(observable+ ' ' + 'Experiment Number' + ' '+str(i))
                    plt.subplot(2,1,2)
                    plt.plot(interploated_temp, weighted_log_difference)
                    plt.ylabel('y')
                    plt.xlabel('Temperature') 
                    plt.savefig(directory_to_save_images+'/'+observable+'_'+'Experiment Number_' +str(i)+'_tempVSy.pdf',dpi=1000,bbox_inches='tight')
                    plt.savefig(directory_to_save_images+'/'+observable+'_'+'Experiment Number_'+str(i)+'_tempVSy.png',dpi=1000,bbox_inches='tight')

                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = parsed_yaml_list[i]['absorbanceCsvWavelengths']
                plt.figure()
                for k,wl in enumerate(wavelengths):
                    difference = exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values- exp['absorbance_model_data'][wl]
                    log_difference = np.log(exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values) - np.log(exp['absorbance_model_data'][wl])
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['absorbance_relative_uncertainty'][k],
                                                             exp['uncertainty']['absorbance_absolute_uncertainty'][k],
                                                             exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp['absorbance_experimental_data'][k])                    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))                    
                    
                    divided = np.divide(difference,un_weighted_uncertainty)
                    numerator =  exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values
                    numerator = numerator.reshape((numerator.shape[0],
                                                                               1))   
                    log_difference = log_difference.reshape((log_difference.shape[0],
                                                                               1))                       
                    
                    weighted_log_difference = np.divide(log_difference,total_uncertainty)
                    
                    plt.figure()
                    plt.subplot(2,1,1)
                    plt.title('Absorbance'+ ' ' +str(wl)+ ' ' + 'Experiment Number' + ' '+str(i))
                    plt.xlabel('Y')

                    plt.hist(log_difference)
                    
                    plt.subplot(2,1,2)
                    plt.hist(weighted_log_difference)
                    plt.xlabel('y')
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.5)
                    plt.savefig(directory_to_save_images+'/'+'Absorbance'+ '_' +str(wl)+ '_' + 'Experiment Number' + '_'+str(i)+'_hist.pdf',dpi=1000,bbox_inches='tight')
                    plt.savefig(directory_to_save_images+'/'+'Absorbance'+ '_' +str(wl)+ '_' + 'Experiment Number' + '_'+str(i)+'_hist.png',dpi=1000,bbox_inches='tight')                     

                    plt.figure()
                    plt.subplot(2,1,1)
                    plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,difference)
                    plt.ylabel('Y')
                    plt.title('Absorbance'+' '+str(wl)+ ' ' + 'Experiment Number' + ' '+str(i))
                    plt.subplot(2,1,2)
                    plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,weighted_log_difference)
                    plt.xlabel('Time (ms)')
                    plt.ylabel('y')
                    plt.savefig(directory_to_save_images+'/'+'Absorbance'+ '_' +str(wl)+ ' ' + 'Experiment Number' + '_'+str(i)+'_timeVSy.pdf',dpi=1000,bbox_inches='tight')
                    plt.savefig(directory_to_save_images+'/'+'Absorbance'+ '_' +str(wl)+ ' ' + 'Experiment Number' + '_'+str(i)+'_timeVSy.png',dpi=1000,bbox_inches='tight')
 
                    
                    plt.figure()                    
                    plt.subplot(2,1,1)
                    interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                    plt.plot(interploated_temp,difference)
                    plt.title('Absorbance'+' '+str(wl)+ ' ' + 'Experiment Number' + ' '+str(i))
                    plt.ylabel('Y')
                    plt.subplot(2,1,2)
                    plt.plot(interploated_temp,weighted_log_difference)
                    plt.xlabel('Temperature') 
                    plt.ylabel('y')
                    plt.savefig(directory_to_save_images+'/'+'Absorbance'+ '_' +str(wl)+ '_' + 'Experiment Number' + '_'+str(i)+'_tempVSy.pdf',dpi=1000,bbox_inches='tight')
                    plt.savefig(directory_to_save_images+'/'+'Absorbance'+ '_' +str(wl)+ '_' + 'Experiment Number' + '_'+str(i)+'_tempVSy.png',dpi=1000,bbox_inches='tight')

        return 
    
    
    def objective_functions(self,exp_dict_list):
        
        def uncertainty_calc(relative_uncertainty,absolute_uncertainty,data,experimental_data):

            if 'W' in list(experimental_data.columns):
                weighting_factor = experimental_data['W'].values
                if 'Relative_Uncertainty' in list(experimental_data.columns):
                    time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                    un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                    total_uncertainty = time_dependent_uncertainty/weighting_factor
                    
                else:
                    length_of_data = data.shape[0]
                    relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty) 
                    un_weighted_uncertainty = copy.deepcopy(relative_uncertainty_array)
                    total_uncertainty = un_weighted_uncertainty/weighting_factor
                
            
            
            elif 'Relative_Uncertainty' in list(experimental_data.columns):  
                
                time_dependent_uncertainty = experimental_data['Relative_Uncertainty'].values
                #do we need to take the natrual log of this?
                time_dependent_uncertainty = np.log(time_dependent_uncertainty+1)
                #do we need to take the natrual log of this?
                length_of_data = data.shape[0]
                un_weighted_uncertainty = copy.deepcopy(time_dependent_uncertainty)
                total_uncertainty = np.divide(time_dependent_uncertainty,(1/length_of_data**.5) )
#                

               
            else:
                length_of_data = data.shape[0]
                relative_uncertainty_array = np.full((length_of_data,1),relative_uncertainty)
                
                if absolute_uncertainty != 0:
                #check if this weighting factor is applied in the correct place 
                #also check if want these values to be the natural log values 
                    absolute_uncertainty_array = np.divide(data,absolute_uncertainty)
                    total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array) + np.square(absolute_uncertainty_array)))
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                     #weighting factor
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )
                
                else:
                    #total_uncertainty = np.log(1 + np.sqrt(np.square(relative_uncertainty_array)))
                    total_uncertainty = relative_uncertainty_array
                    #weighting factor
                    
                    un_weighted_uncertainty = copy.deepcopy(total_uncertainty)
                    total_uncertainty = np.divide(total_uncertainty,(1/length_of_data**.5) )

            #make this return a tuple 
            return total_uncertainty,un_weighted_uncertainty        
        
        objective_function_w_inside = []
        objective_function_w_outside = []
        
        for i,exp in enumerate(exp_dict_list):
            single_exp_dict_w_inside = []
            single_exp_dict_w_outside = []
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['mole_fraction_relative_uncertainty'][observable_counter],
                        exp['uncertainty']['mole_fraction_absolute_uncertainty'][observable_counter],
                        exp['experimental_data'][observable_counter][observable].values,exp['experimental_data'][observable_counter])                    
                    difference = np.log(exp['experimental_data'][observable_counter][observable].values) - np.log((exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6)
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))                    
                    weighted_difference = np.divide(difference,total_uncertainty)
                    square_unweighted_differences = np.square(difference)
                    square_of_differences = np.square(weighted_difference)
                    sum_of_squares = sum(square_of_differences)
                    square_unweighted_differences = difference*total_uncertainty
                    sum_of_squares_w_outside = sum(square_unweighted_differences)
                    single_exp_dict_w_outside.append((sum_of_squares_w_outside[0],np.shape(total_uncertainty[0])))
                    single_exp_dict_w_inside.append((sum_of_squares[0],np.shape(total_uncertainty)[0]))
                    
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['concentration_relative_uncertainty'][observable_counter],
                         exp['uncertainty']['concentration_absolute_uncertainty'][observable_counter],
                         exp['experimental_data'][observable_counter][observable+'_ppm'].values,exp['experimental_data'][observable_counter])    
                    
                    difference = np.log(exp['experimental_data'][observable_counter][observable+'_ppm'].values) - np.log((exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6)
                    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))                    
                    weighted_difference = np.divide(difference,total_uncertainty)
                    square_unweighted_differences = np.square(difference)
                    square_of_differences = np.square(weighted_difference)
                    sum_of_squares = sum(square_of_differences)
                    square_unweighted_differences = difference*total_uncertainty
                    sum_of_squares_w_outside = sum(square_unweighted_differences)
                    single_exp_dict_w_outside.append((sum_of_squares_w_outside[0],np.shape(difference)[0]))
                    single_exp_dict_w_inside.append((sum_of_squares[0],np.shape(difference)[0]))
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp['uncertainty']['absorbance_relative_uncertainty'][k],
                                                             exp['uncertainty']['absorbance_absolute_uncertainty'][k],
                                                             exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp['absorbance_experimental_data'][k])                      
                    difference = np.log(exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values)- np.log(exp['absorbance_model_data'][wl])
                    
                    difference = difference.reshape((difference.shape[0],
                                                                               1))
                                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                                               1))                    
                    weighted_difference = np.divide(difference,total_uncertainty)
                    square_unweighted_differences = np.square(difference)
                    square_of_differences = np.square(weighted_difference)
                    sum_of_squares = sum(square_of_differences)
                    square_unweighted_differences = difference*total_uncertainty
                    sum_of_squares_w_outside = sum(square_unweighted_differences)
                    single_exp_dict_w_outside.append((sum_of_squares_w_outside[0],np.shape(difference)[0]))
                    single_exp_dict_w_inside.append((sum_of_squares[0],np.shape(difference)[0]))
                    
            
            objective_function_w_inside.append(single_exp_dict_w_inside)
            objective_function_w_outside.append(single_exp_dict_w_outside)
        
        return objective_function_w_inside,objective_function_w_outside
    
    def post_proessing_computing_cost_function(self,objective_function_list,experiments_to_consider=[]):
        #experiments_to_consider = [0,1,2,3,4,5,6,7,8]
        
        overall_list = []
        total_sum = []
        for i,value in enumerate(experiments_to_consider):
            for j,tupl in enumerate(objective_function_list[value]):
                overall_list.append(tupl[0])
                total_sum.append(tupl[1])
        

        objective_function_weighted = sum(overall_list)/sum(total_sum)
        objective_function_not_weighted = sum(overall_list)
        return objective_function_weighted , objective_function_not_weighted
    
    
    
    
    
    
    
    
    
    
    
    def plotting_histograms_of_individual_observables_for_hong_data(self,MSI_instance_one,MSI_instance_two,
                                                                    experiments_want_to_plot_data_from,
                                                                    experiments_want_to_plot_data_from_2=[],
                                                                    bins='auto',directory_to_save_images='',csv=''):
        s_shape = MSI_instance_one.S_matrix.shape[1]
        #if MSI_instance_one.k_target_value_S_matrix.any():
            #target_values_for_s = MSI_instance_one.k_target_value_S_matrix
            #s_shape = s_shape+target_values_for_s.shape[0]
        y_shape = MSI_instance_one.y_matrix.shape[0]
        difference = y_shape-s_shape
        y_values = MSI_instance_one.y_matrix[0:difference,0]
        Y_values = MSI_instance_one.Y_matrix[0:difference,0]
        self.lengths_of_experimental_data()

        #plotting_Y Histagrams 
        #edit this part 
        #obserervable_list = []
        
        observables_tottal = []
        
        for i,exp in enumerate(self.exp_dict_list_original):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        
        empty_nested_observable_list_Z = [[] for x in range(len(observables_unique))]
        
        empty_nested_observable_list_Y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_Z_2 = [[] for x in range(len(observables_unique))]

        if bool(experiments_want_to_plot_data_from):
            print('inside here')
            y_values = []
            Y_values = []
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = MSI_instance_two.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = MSI_instance_two.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)

                        temp3 = MSI_instance_two.z_matrix[start:stop,:]
                        empty_nested_observable_list_Z[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]              
          
        if bool(experiments_want_to_plot_data_from_2):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from_2:
                        #print(x)
                        #print(current_observable,'this is current')

                        temp = MSI_instance_one.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y_2[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = MSI_instance_one.y_matrix[start:stop,:]
                        empty_nested_observable_list_y_2[observables_unique.index(current_observable)].append(temp2)

                        temp3 = MSI_instance_one.z_matrix[start:stop,:]
                        empty_nested_observable_list_Z_2[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]                     
                    
        import matplotlib.gridspec as gridspec
    
        fig = plt.figure(figsize=(6,7))

        gs = gridspec.GridSpec(3, 1,height_ratios=[3,3,3],wspace=0.1,hspace=0.1)
        gs.update(wspace=0, hspace=0.7)
        ax1=plt.subplot(gs[0])
        ax2=plt.subplot(gs[1])
        ax3=plt.subplot(gs[2])  
        for i,observable in enumerate(empty_nested_observable_list_Y):
            new_Y_test_2 =[]
            if bool(observable):
                Y_values = np.vstack((observable))
                y_values = np.vstack((empty_nested_observable_list_y[i]))
                z_values = np.vstack((empty_nested_observable_list_Z[i]))
                indecies = np.argwhere(z_values > 100)
                new_y_test = copy.deepcopy(Y_values)
                new_y_test = np.delete(new_y_test,indecies)
             #   print(indecies.shape)
            #    print(indecies)
            #    print(i)
                if bool(experiments_want_to_plot_data_from_2) and bool(empty_nested_observable_list_y_2[i]):
                    Y_values_2 = np.vstack((empty_nested_observable_list_Y_2[i]))
                    y_values_2 = np.vstack((empty_nested_observable_list_y_2[i]))
                    z_values_2 = np.vstack((empty_nested_observable_list_Z_2[i]))
                    indecies_2 = np.argwhere(z_values_2 > 100)
                    new_Y_test_2 = copy.deepcopy(Y_values_2)
                    new_Y_test_2 = np.delete(new_Y_test_2,indecies_2)
                    

                if i ==0:
                #n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    #ax1.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n,bins_test_1,patches = ax1.hist(new_y_test,bins=bins ,align='mid',density=True,label='#1')
                    ax1.set_xlim(left=-.3, right=.3, emit=True, auto=False)
                    ax1.set_ylim(top=15,bottom=0)

                    ax1.set_xlabel('Y')
                    ax1.set_xlabel('Relative Difference')
                #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    ax1.set_title(str(observables_unique[i]))
                    ax1.set_ylabel('pdf')

                #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                   # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        #ax1.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        ax1.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='#2')

                    if bool(csv):
                        df = pd.read_csv(csv)
                        #ax1.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')
                        #ax1.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='#3')

                    ax1.legend()
                if i ==1:
                #n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n,bins_test_2,patches = ax2.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    ax2.set_xlim(left=-.08, right=.08, emit=True, auto=False)
                    ax2.set_ylim(top=28,bottom=0)
                    ax2.set_xlabel('Y')
                    ax2.set_xlabel('Relative Difference')
                #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    #ax2.set_title(str(observables_unique[i]))
                    ax2.set_title(r'H$_2$O')
                    ax2.set_ylabel('pdf')
                    
                #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                   # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        ax2.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        
                    if bool(csv):
                        df = pd.read_csv(csv)
                        #ax2.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')
                
                
                if i ==3:
                #n, bins2, patches = plt.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    n,bins_test_3,patches = ax3.hist(new_y_test,bins=bins ,align='mid',density=True,label='Hong Experiments')
                    ax3.set_xlim(left=-.15, right=.15, emit=True, auto=False)
                    ax3.set_ylim(top=12,bottom=0)

                    ax3.set_xlabel('Y')
                    ax3.set_xlabel('Relative Difference')
                    ax3.set_ylabel('pdf')
                #plt.title(str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from))
                    ax3.set_title(str(observables_unique[i]))
                    ax3.set_title('Absorbance '+ str(observables_unique[i])+ ' nm')

                #plt.ylabel('normalized')
                    if bool(experiments_want_to_plot_data_from_2):
                        print('inside here')
                        print(experiments_want_to_plot_data_from_2)
                   # plt.hist(Y_values_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        ax3.hist(new_Y_test_2,bins=bins ,align='mid',density=True,alpha=0.5,label='Extra Experiments')
                        
                    if bool(csv):
                        df = pd.read_csv(csv)
                        #ax3.hist(df[str(observables_unique[i])+'_Y'].dropna()*-1,bins=bins ,align='mid',density=True,alpha=0.5,label='Hong vs. Hong')                        
                
                
                

                    #plt.savefig(directory_to_save_images+'/'+str(observables_unique[i])+'_Including Experiments_'+ str(experiments_want_to_plot_data_from)+'_Yy_hist_2_normalized.pdf',dpi=1000,bbox_inches='tight')    
        
    
    
    
    
    def plotting_T_and_time_full_simulation_individual_observables_for_hong_data(self,MSI_instance_one,MSI_instance_two,experiments_want_to_plot_data_from,
                                                                             bins='auto',
                                                                             directory_to_save_images='',csv='',experiments_want_to_plot_data_from_2=[]):
#working_here

        observables_tottal = []
        for i,exp in enumerate(self.exp_dict_list_original):
            observable_counter=0
            single_experiment = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                
                if observable == None:
                    continue
                
                if observable in exp['mole_fraction_observables']:
                    single_experiment.append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:

                    single_experiment.append(observable)
                    
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_experiment.append(wl)
                    
            observables_tottal.append(single_experiment)
        
        observables_flatten = [item for sublist in observables_tottal for item in sublist]
        from collections import OrderedDict
        observables_unique = list(OrderedDict.fromkeys(observables_flatten))
        
        empty_nested_observable_list_Y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_z = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from:
                        temp = MSI_instance_one.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = MSI_instance_one.y_matrix[start:stop,:]
                        empty_nested_observable_list_y[observables_unique.index(current_observable)].append(temp2)
                        
                        temp3 = MSI_instance_one.z_matrix[start:stop,:]
                        empty_nested_observable_list_z[observables_unique.index(current_observable)].append(temp3)
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_original):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
####################################################################################################################################################################################################################    
        empty_nested_observable_list_Y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_y_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_time_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_temperature_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_initial_temperature_2 = [[] for x in range(len(observables_unique))]
        empty_nested_observable_list_z_2 = [[] for x in range(len(observables_unique))]


        if bool(experiments_want_to_plot_data_from_2):
            start = 0
            stop = 0 
            for x in range(len(self.simulation_lengths_of_experimental_data)):
                for y in range(len(self.simulation_lengths_of_experimental_data[x])):
                    current_observable = observables_tottal[x][y]
                    stop = self.simulation_lengths_of_experimental_data[x][y] + start
                    if x in experiments_want_to_plot_data_from_2:
                        temp = MSI_instance_two.Y_matrix[start:stop,:]
                        empty_nested_observable_list_Y_2[observables_unique.index(current_observable)].append(temp)
                        
                        temp2 = MSI_instance_two.y_matrix[start:stop,:]
                        empty_nested_observable_list_y_2[observables_unique.index(current_observable)].append(temp2)
                        
                        temp3 = MSI_instance_two.z_matrix[start:stop,:]
                        empty_nested_observable_list_z_2[observables_unique.index(current_observable)].append(temp3)
                        
                        
                        start = start + self.simulation_lengths_of_experimental_data[x][y]
                    else:
                        start = start + self.simulation_lengths_of_experimental_data[x][y]  


        for i,exp in enumerate(self.exp_dict_list_original):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if i in experiments_want_to_plot_data_from_2:
                    if observable in exp['mole_fraction_observables']:
                        empty_nested_observable_list_time_2[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
                        observable_counter+=1
                        
                    if observable in exp['concentration_observables']:
                        empty_nested_observable_list_time_2[observables_unique.index(observable)].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                        interploated_temp = np.interp(exp['experimental_data'][observable_counter]['Time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(observable)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(observable)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])
    
    
                        observable_counter+=1
            if i in experiments_want_to_plot_data_from_2:        
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    for k,wl in enumerate(wavelengths):
                        empty_nested_observable_list_time_2[observables_unique.index(wl)].append(exp['absorbance_experimental_data'][k]['time']*1e3)
    
                        interploated_temp = np.interp(exp['absorbance_experimental_data'][k]['time'],exp['simulation'].timeHistories[0]['time'],exp['simulation'].timeHistories[0]['temperature'])
                        empty_nested_observable_list_temperature_2[observables_unique.index(wl)].append(interploated_temp)
                        empty_nested_observable_list_initial_temperature_2[observables_unique.index(wl)].append([self.exp_dict_list_original[i]['simulation'].temperature]*np.shape(interploated_temp)[0])    
###################################################################################################################################################################################################################  
        
        x = np.arange(10)
        ys = [i+x+(i*x)**2 for i in range(10)]
        colors=cm.rainbow(np.linspace(0,1,30))

        #colors = cm.rainbow(np.linspace(0, 1, len(ys)))
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(6,7))
        gs = gridspec.GridSpec(3, 1,height_ratios=[3,3,3],wspace=0.025,hspace=0.1)
        gs.update(wspace=0, hspace=0.7)
        ax1=plt.subplot(gs[0])
        ax2=plt.subplot(gs[1])
        ax3=plt.subplot(gs[2]) 
        
        fig2 = plt.figure(figsize=(6,7))
        gs2 = gridspec.GridSpec(3, 1,height_ratios=[3,3,3],wspace=0.025,hspace=0.1)
        gs2.update(wspace=0, hspace=0.7)
        ax4=plt.subplot(gs2[0])
        ax5=plt.subplot(gs2[1])
        ax6=plt.subplot(gs2[2]) 


        
        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])
            if bool(observable):

                print(x)
                if x ==0:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax1.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=1,color='g',zorder=4,label='_nolegend_')      
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                            
                        
                            z_values = empty_nested_observable_list_z[x][y]
                            indecies = np.argwhere(z_values > 100)
                            new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                            new_y_test = np.delete(new_y_test,indecies)
                            new_time_test = copy.deepcopy(empty_nested_observable_list_time[x][y])
                            new_time_test = new_time_test.values
                            new_time_test = new_time_test.reshape((new_time_test.shape[0],1))
                            new_time_test  = np.delete(new_time_test,indecies)
                            ax1.scatter(new_y_test,new_time_test, c='#1f77b4',alpha=1)
                            ax1.set_xlabel('Relative Difference')
                            ax1.set_ylabel('Time (ms)')
                            
    
                            if y<length_of_2nd_list:
                                z_values_2 = empty_nested_observable_list_z_2[x][y]
                                indecies_2 = np.argwhere(z_values_2 > 100)
                                new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                                new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                                new_time_test_2 = copy.deepcopy(empty_nested_observable_list_time_2[x][y])
                                new_time_test_2 = new_time_test_2.values
                                new_time_test_2 = new_time_test_2.reshape((new_time_test_2.shape[0],1))
                                new_time_test_2  = np.delete(new_time_test_2,indecies_2)
                                ax1.scatter(new_y_test_2,new_time_test_2,color='orange',zorder=3,alpha=.25)
    
                            ax1.set_title(observables_unique[x])
                            ax1.set_xlim(left=-.25, right=.25, emit=True, auto=False)
                    ax1.scatter([],[],c='#1f77b4',label='#1')
                    ax1.scatter([],[],color='orange',label='#2')                                
                    #ax1.scatter([],[],color='green',label='#3')
                    ax1.legend(frameon=False)




                if x ==1:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax2.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=1,color='g',zorder=4)      
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                            
                        
                            z_values = empty_nested_observable_list_z[x][y]
                            indecies = np.argwhere(z_values > 100)
                            new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                            new_y_test = np.delete(new_y_test,indecies)
                            new_time_test = copy.deepcopy(empty_nested_observable_list_time[x][y])
                            new_time_test = new_time_test.values
                            new_time_test = new_time_test.reshape((new_time_test.shape[0],1))
                            new_time_test  = np.delete(new_time_test,indecies)
                            ax2.scatter(new_y_test,new_time_test, c='#1f77b4',alpha=1)
                            ax2.set_xlabel('Relative Difference')
                            ax2.set_ylabel('Time (ms)')
                            ax2.set_xlim(left=-.09, right=.09, emit=True, auto=False)
                            
    
                            if y<length_of_2nd_list:
                                z_values_2 = empty_nested_observable_list_z_2[x][y]
                                indecies_2 = np.argwhere(z_values_2 > 100)
                                new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                                new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                                new_time_test_2 = copy.deepcopy(empty_nested_observable_list_time_2[x][y])
                                new_time_test_2 = new_time_test_2.values
                                new_time_test_2 = new_time_test_2.reshape((new_time_test_2.shape[0],1))
                                new_time_test_2  = np.delete(new_time_test_2,indecies_2)
                                ax2.scatter(new_y_test_2,new_time_test_2,color='orange',zorder=3,alpha=.25)
    
                            #ax2.set_title(observables_unique[x])
                            ax2.set_title(r'H$_2$O')
                            
                            


                if x ==3:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax3.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_time'].dropna()*1e3,alpha=1,color='g',zorder=4,)      
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                            
                        
                            z_values = empty_nested_observable_list_z[x][y]
                            indecies = np.argwhere(z_values > 100)
                            new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                            new_y_test = np.delete(new_y_test,indecies)
                            new_time_test = copy.deepcopy(empty_nested_observable_list_time[x][y])
                            new_time_test = new_time_test.values
                            new_time_test = new_time_test.reshape((new_time_test.shape[0],1))
                            new_time_test  = np.delete(new_time_test,indecies)
                            ax3.scatter(new_y_test,new_time_test,c='#1f77b4',alpha=1)
                            ax3.set_xlabel('Relative Difference')
                            ax3.set_ylabel('Time (ms)')
                            ax3.set_xlim(left=-.3, right=.3, emit=True, auto=False)


                            
    
                            if y<length_of_2nd_list:
                                z_values_2 = empty_nested_observable_list_z_2[x][y]
                                indecies_2 = np.argwhere(z_values_2 > 100)
                                new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                                new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                                new_time_test_2 = copy.deepcopy(empty_nested_observable_list_time_2[x][y])
                                new_time_test_2 = new_time_test_2.values
                                new_time_test_2 = new_time_test_2.reshape((new_time_test_2.shape[0],1))
                                new_time_test_2  = np.delete(new_time_test_2,indecies_2)
                                ax3.scatter(new_y_test_2,new_time_test_2,color='orange',zorder=3,alpha=.25)
    
                            #ax3.set_title(observables_unique[x])
                            ax3.set_title('Absorbance '+ str(observables_unique[x])+str(' nm'))

                #fig.savefig(directory_to_save_images+'/'+'Three_pannel_plot_'+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_time.pdf',dpi=1000,bbox_inches='tight')

        for x,observable in enumerate(empty_nested_observable_list_Y):
            length_of_2nd_list = len(empty_nested_observable_list_Y_2[x])

            if bool(observable):
                
                if x==0:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax4.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),label='_nolegend_',alpha=1,color='green',zorder=4)                    
                    
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        
                        z_values = empty_nested_observable_list_z[x][y]
                        indecies = np.argwhere(z_values > 100)
                        new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                        new_y_test = np.delete(new_y_test,indecies)
                        new_temperature_test = copy.deepcopy(empty_nested_observable_list_temperature[x][y])
                        new_temperature_test = new_temperature_test.reshape((new_temperature_test.shape[0],1))
                        new_temperature_test  = np.delete(new_temperature_test,indecies)                    
                        
                        ax4.scatter(new_y_test,new_temperature_test,c='#1f77b4',alpha=1,label='_nolegend_')
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        ax4.set_xlabel('Relative Difference')
                        ax4.set_ylabel('Temperature (K)')
                        ax4.set_title(observables_unique[x])
                        ax4.set_xlim(left=-.25, right=.25, emit=True, auto=False)
                        
                        if y<length_of_2nd_list:
                            z_values_2 = empty_nested_observable_list_z_2[x][y]
                            indecies_2 = np.argwhere(z_values_2 > 100)
                            new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                            new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                            new_temperature_test_2 = copy.deepcopy(empty_nested_observable_list_temperature_2[x][y])
                            new_temperature_test_2 = new_temperature_test_2.reshape((new_temperature_test_2.shape[0],1))
                            new_temperature_test_2  = np.delete(new_temperature_test_2,indecies_2)                              
                            ax4.scatter(new_y_test_2,new_temperature_test_2,color='orange',zorder=3,alpha=.25,label='_nolegend_')
                        
                    ax4.scatter([],[],c='#1f77b4',label='#1')
                    ax4.scatter([],[],c='orange',label='#2')
                    #ax4.scatter([],[],c='green',label='#3')
                    ax4.legend(frameon=False)
   
                                    
                if x==1:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax5.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),alpha=1,color='green',zorder=4)                    
                    
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        
                        z_values = empty_nested_observable_list_z[x][y]
                        indecies = np.argwhere(z_values > 100)
                        new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                        new_y_test = np.delete(new_y_test,indecies)
                        new_temperature_test = copy.deepcopy(empty_nested_observable_list_temperature[x][y])
                        new_temperature_test = new_temperature_test.reshape((new_temperature_test.shape[0],1))
                        new_temperature_test  = np.delete(new_temperature_test,indecies)                    
                        
                        ax5.scatter(new_y_test,new_temperature_test,c='#1f77b4',alpha=1)
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        ax5.set_xlabel('Relative Difference')
                        ax5.set_ylabel('Temperature (K)')
                        #ax5.set_title(observables_unique[x])
                        ax5.set_title(r'H$_2$O')
                        ax5.set_xlim(left=-.09, right=.09, emit=True, auto=False)

                        
                        if y<length_of_2nd_list:
                            z_values_2 = empty_nested_observable_list_z_2[x][y]
                            indecies_2 = np.argwhere(z_values_2 > 100)
                            new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                            new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                            new_temperature_test_2 = copy.deepcopy(empty_nested_observable_list_temperature_2[x][y])
                            new_temperature_test_2 = new_temperature_test_2.reshape((new_temperature_test_2.shape[0],1))
                            new_temperature_test_2  = np.delete(new_temperature_test_2,indecies_2)                              
                            ax5.scatter(new_y_test_2,new_temperature_test_2,color='orange',zorder=3,alpha=.15)

                if x==3:
                    if bool(csv):
                        df = pd.read_csv(csv)
                        ax6.scatter(df[str(observables_unique[x])+'_Y'].dropna()*-1,df[str(observables_unique[x])+'_Temperature'].dropna(),alpha=1,color='green',zorder=4)                    
                    
                    for y,array in enumerate(empty_nested_observable_list_Y[x]):
                        
                        z_values = empty_nested_observable_list_z[x][y]
                        indecies = np.argwhere(z_values > 100)
                        new_y_test = copy.deepcopy(empty_nested_observable_list_Y[x][y])
                        new_y_test = np.delete(new_y_test,indecies)
                        new_temperature_test = copy.deepcopy(empty_nested_observable_list_temperature[x][y])
                        new_temperature_test = new_temperature_test.reshape((new_temperature_test.shape[0],1))
                        new_temperature_test  = np.delete(new_temperature_test,indecies)                    
                        
                        ax6.scatter(new_y_test,new_temperature_test,label='Experiment_'+str(x)+'_observable_'+str(y),c='#1f77b4',alpha=1)
                        #plt.legend(ncol=2,bbox_to_anchor=(1, 0.5))
                        ax6.set_xlabel('Relative Difference')
                        ax6.set_ylabel('Temperature (K)')
                        ax6.set_title('Absorbance '+ str(observables_unique[x])+str(' nm'))
                        ax6.set_xlim(left=-.3, right=.3, emit=True, auto=False)

                        if y<length_of_2nd_list:
                            z_values_2 = empty_nested_observable_list_z_2[x][y]
                            indecies_2 = np.argwhere(z_values_2 > 100)
                            new_y_test_2 = copy.deepcopy(empty_nested_observable_list_Y_2[x][y])
                            new_y_test_2 = np.delete(new_y_test_2,indecies_2)
                            new_temperature_test_2 = copy.deepcopy(empty_nested_observable_list_temperature_2[x][y])
                            new_temperature_test_2 = new_temperature_test_2.reshape((new_temperature_test_2.shape[0],1))
                            new_temperature_test_2  = np.delete(new_temperature_test_2,indecies_2)                              
                            ax6.scatter(new_y_test_2,new_temperature_test_2,color='orange',zorder=3,alpha=.25)
  
                #fig2.savefig(directory_to_save_images+'/'+'Three_pannel_plot'+'_Including Experiments_'+str(experiments_want_to_plot_data_from)+'_Yy_vs_initial_temperature.pdf',dpi=1000,bbox_inches='tight')

    def plotting_3_figre_observables(self,experiment_number_want_to_plot,sigmas_original=[],sigmas_optimized=[]):
        
        
        df = pd.read_csv('MSI/data/klip_optimization_comparison/graph_read_hong_data_1072.csv')
        df = pd.read_csv('MSI/data/klip_optimization_comparison/graph_read_hong_1283.csv')

        for i,exp in enumerate(self.exp_dict_list_optimized):
            if i==experiment_number_want_to_plot:
                observable_counter=0
                plt.figure(figsize=(10,10))
               # plt.figure()
                for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                    if observable == None:
                        continue

                        
                    if observable in exp['concentration_observables']:
                        if observable == 'H2O':
                            plt.subplot(3,1,1)
                            #plt.xlim(-.005,3.5)
                            #plt.ylim(1000,3500)
                            
                            plt.xlim(-.005,1)
                            plt.plot(df['hong_h2o_time'],df['hong_h2o'],'g:')
                    
                        if observable =='OH':
                            plt.subplot(3,1,2)
                            #plt.xlim(-.005,1)
                            #plt.ylim(0,25)
                           
                            plt.xlim(-.005,1)
                            plt.plot(df['hong_oh_time'],df['hong_oh'],'g:')
                            #plt.plot(df['hong_oh_time_2'],df['hong_oh_2'],'g:')
                        if observable == 'H2O2':
                            plt.subplot(3,1,1)
                            plt.xlim(-.005,.8)

                        plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable]*1e6,'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable]*1e6,'r',label= "$\it{a}$ $\it{priori}$ model")
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Hong et al. experiment') 
                       # plt.xlabel('Time (ms)')
                        if observable == 'H2O':
                            plt.ylabel(r'H$_2$O [ppm]', fontsize=18)
                        if observable == 'OH':
                            plt.ylabel('OH [ppm]',fontsize=18)
                        if observable == 'H2O2':
                            plt.ylabel(r'H$_2$O$_2$ [ppm]', fontsize=18)
                            plt.xlabel('Time [ms]',fontsize=18)

                        else:                        
                            plt.ylabel(str(observable) +' '+ '[ppm]')
                        #plt.title('Experiment_'+str(i+1))
                        
                        if bool(sigmas_optimized)==True:
                            high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                            high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            low_error_optimized = np.exp(np.array(sigmas_optimized[i][observable_counter])*-1)
                            low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,'b--')
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,'b--')                    
                            
        
        
                            high_error_original = np.exp(sigmas_original[i][observable_counter])
                            high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            low_error_original = np.exp(np.array(sigmas_original[i][observable_counter])*-1)
                            low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            

                        key_list = []
                        for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():
                            
                            #plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                            key_list.append(key)
                       
                        #plt.legend(handlelength=3)
                        #plt.legend(ncol=1)
                        sp = '_'.join(key_list)
                        #print(sp)
                        #plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K'+'_'+str(self.exp_dict_list_original[i]['simulation'].pressure)+'_'+sp+'_'+'.pdf', bbox_inches='tight')
                        
                        #stub
                        plt.tick_params(direction='in')
                        
                        
                        #plt.savefig('Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                        #plt.savefig('Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)
                        #plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight',dpi=1000)
    
    

                        observable_counter+=1
                        
    
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    plt.subplot(3,1,3)
                    #plt.xlim(-.005,3)
                    
                    plt.xlim(-.005,.2)
                    #plt.xlim(-.005,.6)
                    #plt.xlim(-.005,1.5)

                    plt.plot(df['hong_abs_time'],df['hong_abs'],'g:',zorder=10,label='Hong et al. model')
                    


                    for k,wl in enumerate(wavelengths):
                        if wl == 227:
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Hong et al. experiment')
                        if wl == 215:
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Kappel et al. experiment')
    
                        plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['absorbance_calculated_from_model'][wl],'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['absorbance_calculated_from_model'][wl],'r',label= "$\it{a}$ $\it{priori}$ model")
                        #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')
                        plt.xlabel('Time [ms]',fontsize=18)
                        #plt.ylabel('Absorbance'+' '+str(wl)+' nm')
                        plt.ylabel('Absorbance',fontsize=18)
                        #plt.title('Experiment_'+str(i+1))
                        
                        if bool(sigmas_optimized)==True:
                            high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])
                            high_error_optimized = np.multiply(high_error_optimized,exp['absorbance_model_data'][wl])
                            low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                            low_error_optimized = np.multiply(low_error_optimized,exp['absorbance_model_data'][wl])
                            
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,high_error_optimized,'b--')
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,low_error_optimized,'b--')
                            
                            high_error_original = np.exp(sigmas_original[i][observable_counter])
                            high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['absorbance_model_data'][wl])
                            low_error_original =  np.exp(sigmas_original[i][observable_counter]*-1)
                            low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['absorbance_model_data'][wl])
                            
                            
                            #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,high_error_original,'r--')
                            #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,low_error_original,'r--')
    
    #                    if bool(sigmas_optimized)==True and  i+1 == 11:    
    #                        plt.ylim(top=.35)
                        
                        #start here
                        
                        #plt.plot([],'w' ,label= 'T:'+ str(self.exp_dict_list_original[i]['simulation'].temperature))
                        #plt.plot([],'w', label= 'P:'+ str(self.exp_dict_list_original[i]['simulation'].pressure))
                        #for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():                        
                            #plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                            
    
                        #plt.legend(handlelength=3)
                        plt.legend(ncol=2, prop={'size': 12})
                        plt.savefig(self.working_directory+'/'+'Exp_'+str(experiment_number_want_to_plot)+'_three_figure_obervables.pdf', bbox_inches='tight')
                        plt.tick_params(direction='in')
   

                        #plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                        #plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)