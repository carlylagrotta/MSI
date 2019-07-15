#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 11:56:22 2018

@author: carly
"""
import numpy as np
import pandas as pd

class Master_Equation(object):
    def __init__(self):
        self.matrix = None

    def multiply_by_sensitivites(self,array,array_of_sensitivities,pressure_and_temp_array):
        sensitivity_multiplied_array = np.zeros((array_of_sensitivities.shape[0],array_of_sensitivities.shape[1],pressure_and_temp_array.shape[0]))
        for sensitivity in np.nditer(array_of_sensitivities,order='F'):
            i,j= np.where(array_of_sensitivities == sensitivity)
            temp_coef = []
            pres_coef = []
            for p,value in enumerate(pressure_and_temp_array['temperature']):
                
                t_coef = self.chebyshev_specific_poly(i[0],self.calc_reduced_T(value))
                temp_coef.append(t_coef)
                p_coef = self.chebyshev_specific_poly(j[0],self.calc_reduced_P(pressure_and_temp_array['pressure'][p]))
                pres_coef.append(p_coef)
                    
                    
            temp_coef = np.array(temp_coef)
            temp_coef = temp_coef.reshape((temp_coef.shape[0],1))
            pres_coef = np.array(pres_coef)
            pres_coef = pres_coef.reshape((pres_coef.shape[0],1))
                
            mapped_array = np.multiply(pres_coef,temp_coef)
            mapped_array = np.multiply(mapped_array,array)   
            sensitivity_multiplied_array[i[0],j[0],:] = mapped_array.flatten()
            
        return sensitivity_multiplied_array
    
    def array_reshape(self,three_d_array):
        temp_array = []
        for row in range(three_d_array.shape[0]):
            for column in range(three_d_array.shape[1]):
                alpha = three_d_array[row,column,:]
                alpha = alpha.reshape((alpha.shape[0],1))
                temp_array.append(alpha)
        temp_array = np.hstack((temp_array))
        return temp_array

    
    def chebyshev_specific_poly(self,order_needed,value):
        #this function starts indexing at 1 
        order_needed  = order_needed + 2 
        coef = [1 for value in range(order_needed)]
        
        x = np.polynomial.chebyshev.chebval(value,
                                       coef,
                                       tensor = False)
        
        y = np.polynomial.chebyshev.chebval(value,
                                       coef[0:-1],
                                       tensor = False)  
        return x-y
        
    def calc_reduced_T(self,T,T_min=200,T_max=2500):
        numerator = (2*(T**-1))-((T_min)**-1) - ((T_max)**-1)
        denominator = ((T_min)**-1) - ((T_max)**-1)
        T_reduced = np.divide(numerator,denominator)
           
        return T_reduced
        
    def calc_reduced_P(self,P,P_min=1013.25,P_max=1.013e+6):
        numerator = 2*np.log(P) - np.log(P_min) - np.log(P_max)
        denominator = np.log(P_max) - np.log(P_min)
        P_reduced = np.divide(numerator,denominator)
        return P_reduced


    
    
    def map_to_alpha(self,sensitivty_dict:dict,
                 exp_dict_list:list,
                 parsed_yaml_file_list,
                 master_equation_reactions:list):
    
        nested_list = []
        def slicing_out_reactions(reaction_string,array):
            reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
            index_of_reaction_in_cti = reactions_in_cti_file.index(reaction_string)
            column_of_array = array[:,index_of_reaction_in_cti]
            column_of_array = column_of_array.reshape((column_of_array.shape[0],
                                                                  1))                  
            return column_of_array
        mapped_to_alpha_full_simulation = []
        for i, exp in enumerate(exp_dict_list):
            simulation = []
            single_experiment = []
            
            if parsed_yaml_file_list[i]['moleFractionObservables'][0] != None or parsed_yaml_file_list[i]['concentrationObservables'][0] != None:
                As = exp['ksens']['A']
                for xx,observable in enumerate(As):
                    temp = []
                    observable_list = []
                    for reaction in master_equation_reactions:
                        column = slicing_out_reactions(reaction,observable)
                        single_reaction_array = self.array_reshape(self.multiply_by_sensitivites(column,sensitivty_dict[reaction][0],exp['simulation'].pressureAndTemperatureToExperiment[xx]))
                        temp.append(single_reaction_array)
                        observable_list.append(single_reaction_array)
                        
                    simulation.append(observable_list)
                   
                    single_experiment.append(np.hstack((temp)))
                    
                 
            if 'absorbance_observables' in list(exp.keys()):
                wavelengths = parsed_yaml_file_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    temp = []
                    observable_list = []
                    for reaction in master_equation_reactions:
                        column = slicing_out_reactions(reaction,exp['absorbance_ksens'][wl][0])
                        single_reaction_array = self.array_reshape(self.multiply_by_sensitivites(column,sensitivty_dict[reaction][0],exp['time_history_interpolated_against_abs'][wl]))
                        temp.append(single_reaction_array)
                        observable_list.append(single_reaction_array)
                        
                    single_experiment.append(np.hstack((temp)))
                    simulation.append(observable_list)
                   
                       
            

            nested_list.append(simulation)        
            mapped_to_alpha_full_simulation.append(np.vstack((single_experiment)))
        
       
        
            
        return mapped_to_alpha_full_simulation,nested_list
    
    
    
    
    def map_to_S(self,mapped_to_alpha_full_simulation,
                 sensitivity_dict:dict,
                 master_equation_reactions:list):
        MP_nested_list_full_experiment = []
        for simulation in range(len(mapped_to_alpha_full_simulation)):
            simulations = []
            for observable in range(len(mapped_to_alpha_full_simulation[simulation])):
                observables = []
                for reaction_number,reaction in enumerate(mapped_to_alpha_full_simulation[simulation][observable]):
                    reaction = []
                    for molecular_parameter,array_of_sensitivities in enumerate(sensitivity_dict[master_equation_reactions[reaction_number]]):
                        zero_array = np.zeros((np.shape(mapped_to_alpha_full_simulation[simulation][observable][reaction_number])))
                        column_counter = 0
                        for sensitivity in np.nditer(sensitivity_dict[master_equation_reactions[reaction_number]][molecular_parameter],order='F'):
                            
                            i,j= np.where(array_of_sensitivities == sensitivity)
                            mapped_alpha_array = mapped_to_alpha_full_simulation[simulation][observable][reaction_number]
                            column_to_multiply_by_sens = mapped_alpha_array[:,column_counter]
                            column_to_multiply_by_sens = np.multiply(column_to_multiply_by_sens,sensitivity)
                            zero_array[:,column_counter] = column_to_multiply_by_sens
                            column_counter+=1
                        reaction.append(zero_array)
                    observables.append(reaction)
                simulations.append(observables)
            MP_nested_list_full_experiment.append(simulations)
            
        MP = []   
        for simulation in range(len(MP_nested_list_full_experiment)):
            simulations = []
            for observable in range(len(MP_nested_list_full_experiment[simulation])):
                observables = []
                for reaction in range(len(MP_nested_list_full_experiment[simulation][observable])):
                    reactions = []
                    for molecular_parameter,MP_array in enumerate(MP_nested_list_full_experiment[simulation][observable][reaction]):
                        #how do we sum these things up 
                        print(np.shape(MP_array))
                        Column_vector  = np.sum(MP_array,axis = 1)
                        Column_vector = Column_vector.reshape((Column_vector.shape[0],1))
                        reactions.append(Column_vector)
                    observables.append(np.hstack(reactions))
                simulations.append(np.vstack(observables))
            MP.append(np.vstack(simulations))
        MP = np.vstack((MP))
        print(np.shape(MP))
        return MP    

                


    def surrogate_model_molecular_parameters_chevy( master_equation_sensitivites:dict,
                                              master_equation_reactions:list,
                                              delta_x_molecular_params_by_reaction_dict,
                                              exp_dict_list):
        
        reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
        number_of_reactions = len(reactions_in_cti_file)
        
        updates = []
        for reaction in range(len(master_equation_reactions)):
            temp = []
            for i,sensitivty in enumerate(master_equation_sensitivites[reaction]):
                temp.append(np.multiply(sensitivty*delta_x_molecular_params_by_reaction_dict['R_'+str(reaction)][i]))
            
            updates.append(sum(temp))
            
        Keys = []                    
        for x in np.arange(number_of_reactions-len(master_equation_reactions),number_of_reactions):
            Keys.append('r'+str(x))

        MP = dict(zip(Keys,updates))
        #units 
        return MP
    
    
    
    
    