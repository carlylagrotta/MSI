import numpy as np
import pandas as pd
import MSI.master_equation.master_equation as meq 
import copy

class OptMatrix(object):
    def __init__(self):
        self.S_matrix = None
        self.s_matrix = None
        self.Y_matrix = None
        self.y_matrix = None
        self.z_matrix = None
        self.delta_X  = None
        self.X = None
        self.sigma = None
 
#    #loads one experiment into self.matrix. Decides padding based on previous matrix or handle based on total exp num?
    def load_S(self, exp_dict_list:list,parsed_yaml_list:list,
               dk=.01,
               master_equation_reactions = [],
               mapped_master_equation_sensitivites=np.array(()),
               master_equation_uncertainty_df = None,
               master_equation_flag = False):
        
        

        
        #preprocessing for padding
        num_exp = len(exp_dict_list)
        pert_coef = {} #build a dict matching pert_coef to their experiment and wavelength.
                       #length of the dict gives padding information
        list_to_keep_order_of_coef = []
        for exp in exp_dict_list:
            if 'perturbed_coef' not in exp.keys():
                continue
            perturbed_for_exp = exp['perturbed_coef']
            for x in perturbed_for_exp:
                if x[0][2] not in pert_coef.keys():
                    pert_coef[x[0][2]] = [x[1]]
                else:

                    pert_coef[x[0][2]].append(x[1])

                    
                if x[0][2] not in list_to_keep_order_of_coef:
                    list_to_keep_order_of_coef.append(x[0][2])
            
        num_ind_pert_coef = len(pert_coef)
        #print(pert_coef.keys())
        
        #print(num_ind_pert_coef," sigmas")
        #establish # of independent pert before hand, to proper pad the observables, put in list, make a dict of cc,
        # values will be a list of tabs data?
        # use the list to get the padding size
        k_sens_for_whole_simulation = []
        p_sens_for_whole_simulation = []
        abs_coef_sens_for_whole_simulation = []
        
        temps = []
        for i,exp in enumerate(exp_dict_list):
            ttl_kinetic_observables_for_exp = []
            obs_counter =0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                
                single_obs_matrix = np.hstack((exp['ksens']['A'][obs_counter],
                                        exp['ksens']['N'][obs_counter],
                                        exp['ksens']['Ea'][obs_counter]))
                
               
                ttl_kinetic_observables_for_exp.append(single_obs_matrix)
                obs_counter +=1
                
            if 'perturbed_coef' in exp.keys():
                wavelengths = parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    single_obs_matrix = np.hstack((exp['absorbance_ksens'][wl][0],
                                                   exp['absorbance_ksens'][wl][1],
                                                   exp['absorbance_ksens'][wl][2]))
                    
                    ttl_kinetic_observables_for_exp.append(single_obs_matrix)
                                       
            ttl_kinetic_observables_for_exp = np.vstack((ttl_kinetic_observables_for_exp))              
            k_sens_for_whole_simulation.append(ttl_kinetic_observables_for_exp)
            ####vstack  ttl_kinetic_observables_for_exp   and append somwehre else
            if exp['simulation'].physicalSens ==1:
                ttl_phsycal_obs_for_exp = []
                for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                    obs_counter = 0
                    if observable == None:
                        continue
                    temperature_sensitivity = exp['temperature'][observable].dropna().values
                    temperature_sensitivity = temperature_sensitivity.reshape((temperature_sensitivity.shape[0],
                                                          1))
                    
                    pressure_sensitivity = exp['pressure'][observable].dropna().values
                    pressure_sensitivity = pressure_sensitivity.reshape((pressure_sensitivity.shape[0],
                                                          1))
                    species_sensitivty = []
                    for df in exp['species']:
                        single_species_sensitivty = df[observable].dropna().values
                        single_species_sensitivty = single_species_sensitivty.reshape((single_species_sensitivty.shape[0]
                                                           ,1))
                        species_sensitivty.append(single_species_sensitivty)
                        
                    species_sensitivty = np.hstack((species_sensitivty))
                        
                    
                    single_obs_physical = np.hstack((temperature_sensitivity,pressure_sensitivity,species_sensitivty))
                    ttl_phsycal_obs_for_exp.append(single_obs_physical)
                    obs_counter +=1
                if 'perturbed_coef' in exp.keys():
                    wavelengths = parsed_yaml_list[i]['absorbanceCsvWavelengths']                    
                    for k,wl in enumerate(wavelengths):
                        physical_sens = []
                        for p_sens in exp['absorbance_psens']:                            
                            array = p_sens[wl]
                            array = array.reshape((array.shape[0],1))
                            physical_sens.append(array)
                        physical_sens = np.hstack((physical_sens))
                        ttl_phsycal_obs_for_exp.append(physical_sens)
                    
                ttl_phsycal_obs_for_exp = np.vstack((ttl_phsycal_obs_for_exp))
                p_sens_for_whole_simulation.append(ttl_phsycal_obs_for_exp)
#######################################################################################################################################################               

            
            
            if 'perturbed_coef' in exp.keys():
                ttl_absorbance_obs_for_exp = []
                wavelengths = parsed_yaml_list[i]['absorbanceCsvWavelengths']
                for k,wl in enumerate(wavelengths):
                    perturbed_coefficeints = []
                    index_list = []
                    for xx in range(len(parsed_yaml_list[i]['coupledCoefficients'])):
                        for yy in range(len(parsed_yaml_list[i]['coupledCoefficients'][xx])):
                            ff = parsed_yaml_list[i]['functionalForm'][xx][yy]
                            #temp = list(parsed_yaml_list[i]['coupledCoefficients'][xx][yy])
                            for zz in range(len(parsed_yaml_list[i]['coupledCoefficients'][xx][yy])):
                                temp = list(parsed_yaml_list[i]['coupledCoefficients'][xx][yy])
                                coefficent = parsed_yaml_list[i]['coupledCoefficients'][xx][yy][zz]    
                                if coefficent!=0:
                                    perturbed_coefficent=coefficent+coefficent*dk
                                    if zz==1 and ff =='F':
                                        #change back tab
                                        perturbed_coefficent = coefficent + .01*coefficent
                                    temp[zz] = perturbed_coefficent

                                    key = tuple(temp)

                                    indx = list_to_keep_order_of_coef.index(key)
                                    
                                    index_list.append(indx)
                                    
                                    exp_index_sigma = temps.count(key)
                                    temps.append(key)

                                    array = pert_coef[key][exp_index_sigma][wl]
                                    
                                    array = array.reshape((array.shape[0],1))
                                    perturbed_coefficeints.append(array)

                                
                    
                    
                    missing_sigmas = []  
                    for indp_sigma in range(len(list_to_keep_order_of_coef)):
                        if indp_sigma not in index_list:
                            missing_sigmas.append(indp_sigma)
                            
                    perturbed_coefficents_padded_with_zeros = []
                    count_sigma=0
                    for indp_sigma in range(len(list_to_keep_order_of_coef)):
                        if indp_sigma in missing_sigmas:
                            zero_array = np.zeros((perturbed_coefficeints[0].shape[0],1))
                            perturbed_coefficents_padded_with_zeros.append(zero_array)                            
                        else:
                            perturbed_coefficents_padded_with_zeros.append(perturbed_coefficeints[count_sigma])
                            count_sigma +=1
                    
                    
                        
                        
                        
                        
                    
                    perturbed_coefficents_padded_with_zeros = np.hstack((perturbed_coefficents_padded_with_zeros))
                    ttl_absorbance_obs_for_exp.append(perturbed_coefficents_padded_with_zeros) 
                    
                ttl_absorbance_obs_for_exp = np.vstack((ttl_absorbance_obs_for_exp))
                abs_coef_sens_for_whole_simulation.append(ttl_absorbance_obs_for_exp)
            
                #vstack ttl_absorbance_obs_for_exp and append somehwere else 
            else:
                 abs_coef_sens_for_whole_simulation.append(0)
                
               
######################################################################################################################################################

                
        #assembling the S matrix from the individual experiments 
        #master_equation = False
        if master_equation_flag == True:
            S_ksens = np.vstack((k_sens_for_whole_simulation))           
            A_k = np.hsplit(S_ksens,3)[0]
            N_k = np.hsplit(S_ksens,3)[1]
            Ea_k  = np.hsplit(S_ksens,3)[2]
            
            number_of_master_equation_reactions = len(master_equation_reactions)
            
            A_k = A_k[:,:-number_of_master_equation_reactions]
            N_k = N_k[:,:-number_of_master_equation_reactions]
            Ea_k = Ea_k[:,:-number_of_master_equation_reactions]
            
           
            S_ksens = np.hstack((A_k,N_k,Ea_k))
            #print(np.shape(S_ksens),'this is the shape of the S matrix before MP')
            S_ksens = np.hstack((S_ksens,mapped_master_equation_sensitivites))
            

        else:
            S_ksens = np.vstack((k_sens_for_whole_simulation))
        
        
        def sum_of_zeros(idx,array,column_list):
            rows_behind = array.shape[0]
            rows_infront = array.shape[0]
            columns_behind = sum(column_list[:idx])
            columns_infront = sum(column_list[idx+1:])  
            behind_tuple = (rows_behind,columns_behind)
            infront_tuple = (rows_infront,columns_infront)
            
            return (behind_tuple,infront_tuple)
        
        
        if exp_dict_list[0]['simulation'].physicalSens ==1:
            number_of_columns_in_psens_arrays = []
            number_of_rows_in_psens_arrays=[]
            for i,array in enumerate(p_sens_for_whole_simulation):                
                number_of_rows_in_psens_arrays.append(array.shape[0])
                number_of_columns_in_psens_arrays.append(array.shape[1])
                
            p_sens_whole_simulation_with_padding = []
            for i,array in enumerate(p_sens_for_whole_simulation):
                zero_array_behind = np.zeros(sum_of_zeros(i,array,number_of_columns_in_psens_arrays)[0])
                if zero_array_behind.shape[1] != 0:
                    array = np.hstack((zero_array_behind,array))
                zero_array_infront = np.zeros(sum_of_zeros(i,array,number_of_columns_in_psens_arrays)[1])
                if zero_array_infront.shape[1] != 0:
                    array = np.hstack((array,zero_array_infront))
                
                p_sens_whole_simulation_with_padding.append(array)
            
            S_psens = np.vstack((p_sens_whole_simulation_with_padding))
            
            
            
##############################################################################################         
        absorb_coef_whole_simulation_with_padding = []
        
        for i,exp in enumerate(exp_dict_list):
            single_experiment_absorption = []
            if exp['mole_fraction_observables'][0] != None or exp['concentration_observables'][0] != None:                
                if 'perturbed_coef' not in exp.keys():
                    zero_array_for_observables_padding = np.zeros((number_of_rows_in_psens_arrays[i],
                                           num_ind_pert_coef))    
                    
                    single_experiment_absorption.append(zero_array_for_observables_padding)
                    
                    
            if 'perturbed_coef' in exp.keys():  
                zero_padded_aborption_coef_array = abs_coef_sens_for_whole_simulation[i]   
                combined = abs_coef_sens_for_whole_simulation[i] 

                if exp['mole_fraction_observables'][0] != None or exp['concentration_observables'][0] != None:
                    zero_array_for_observables_padding = np.zeros((number_of_rows_in_psens_arrays[i]-zero_padded_aborption_coef_array.shape[0],
                                       num_ind_pert_coef))  
                    combined = np.vstack((zero_array_for_observables_padding,zero_padded_aborption_coef_array))
                
                
                
                single_experiment_absorption.append(combined)     
            
            single_experiment_absorption = np.vstack((single_experiment_absorption))
            absorb_coef_whole_simulation_with_padding.append(single_experiment_absorption) 
            
        
        
        absorb_coef_whole_simulation_with_padding = np.vstack((absorb_coef_whole_simulation_with_padding))  
        S_abs_coef  = absorb_coef_whole_simulation_with_padding

        

        S_matrix = np.hstack((S_ksens,S_psens,S_abs_coef))
        shape = np.shape(S_matrix)[1]
        #append identy matrix
        identity_matrix = np.identity(shape)
        S_matrix = np.vstack((S_matrix,identity_matrix))
        self.S_matrix = S_matrix
        S_matrix_wo_k_targets = copy.deepcopy(self.S_matrix)
        self.S_matrix_wo_k_targets = S_matrix_wo_k_targets
        S_matrix_df = pd.DataFrame(S_matrix)
        return S_matrix 

                                


    def load_Y(self, exp_dict_list:list,parsed_yaml_file_list:list,
               loop_counter:int = 0,
               X:dict={},
               master_equation_reactions = [],
               master_equation_uncertainty_df = None,
               master_equation_flag = False):
        
        
        def natural_log_difference(experiment,model):
            natural_log_diff = np.log(experiment) - np.log(model)
            return natural_log_diff
        
        Y = []
        Y_data_Frame = []
        for i,exp_dic in enumerate(exp_dict_list):
            counter = 0
            for j,observable in enumerate((exp_dic['mole_fraction_observables']+
                                           exp_dic['concentration_observables'])):
                if observable == None:
                    pass
                else:
                    #if you need to add something with concentration add it here 
                    if 'ppm' in exp_dic['experimental_data'][counter].columns.tolist()[1]:
                        natural_log_diff = natural_log_difference(exp_dic['experimental_data'][counter][observable+'_ppm'].values,
                                                                  (exp_dic['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)*1e6)
                        
                        
                        natural_log_diff =  natural_log_diff.reshape((natural_log_diff.shape[0],
                                                          1))


                        
                    else:
                        natural_log_diff = natural_log_difference(exp_dic['experimental_data'][counter][observable].values,
                                                                  exp_dic['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        natural_log_diff =  natural_log_diff.reshape((natural_log_diff.shape[0],
                                                          1))
                
                    tempList = [observable+'_'+'experiment'+str(i)]*np.shape(natural_log_diff)[0]
                    Y_data_Frame.extend(tempList)
                
                    Y.append(natural_log_diff)
                    counter+=1
            if 'absorbance_observables' in list(exp_dic.keys()):
                wavelengths = parsed_yaml_file_list[i]['absorbanceCsvWavelengths']
                
                for k,wl in enumerate(wavelengths):
                    natural_log_diff = natural_log_difference(exp_dic['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp_dic['absorbance_model_data'][wl])
                    natural_log_diff =  natural_log_diff.reshape((natural_log_diff.shape[0],
                                                  1))
                    
                    tempList = [str(wl)+'_'+'experiment'+'_'+str(i)]*np.shape(natural_log_diff)[0]
                    Y_data_Frame.extend(tempList)
                    
                    
                    Y.append(natural_log_diff)
        
        Y = np.vstack((Y))
        
       
        #YdataFrame = pd.DataFrame({'value': YdataFrame,'ln_difference': Y})
        
        reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
       
            #assembling the target values portion of the Y matrix 
            #getting the size of the cti file from the first simulation because 
            #they all use the same cti file and it shouldn't matter 
            
            
            # add in a conditional statment for if there is master equation data
            #which is getting included in the simulation

        if master_equation_flag ==True:
            A_n_Ea_length = int((len(reactions_in_cti_file) - len(master_equation_reactions))*3)            
            number_of_molecular_parameters_list = []
            for col in master_equation_uncertainty_df:
                number_of_molecular_parameters_list.append(len(master_equation_uncertainty_df[col].dropna().values))
                
            number_of_molecular_parameters = sum(number_of_molecular_parameters_list) 
            #print('we do not have master equation installed yet')
            #subtract out the necessary target values and add the other ones in 
        else:
            A_n_Ea_length = len(reactions_in_cti_file)*3
            
            #addint the zeros to the Y array 
            
            #adding the strings to the dictonary 
            ## making a,n and Ea zero list 
        A_n_Ea_zeros = np.zeros((A_n_Ea_length,1))  
        
        if master_equation_flag ==True:
            molecular_paramter_zeros = np.zeros((number_of_molecular_parameters,1))
        
        
        for variable in range(A_n_Ea_length//3):
            Y_data_Frame.append('A'+'_'+str(variable))
        for variable in range(A_n_Ea_length//3):
            Y_data_Frame.append('n'+'_'+str(variable))
        for variable in range(A_n_Ea_length//3):
            Y_data_Frame.append('Ea'+'_'+str(variable))
            
        if master_equation_flag == True:
            for i,value in enumerate(number_of_molecular_parameters_list):
                for j,parameter in enumerate(range(value)):
                    Y_data_Frame.append('R'+'_'+str(i)+'P'+'_'+str(j))
               
            
        
        
        if loop_counter == 0:
            Y = np.vstack((Y,A_n_Ea_zeros))
            if master_equation_flag ==True:
                Y = np.vstack((Y,molecular_paramter_zeros))
        
             
        else:
            #print('we do not have loop counter installed yet')
            #need to check what we would need to do here 
            #should be tottal X ?
            
            #clean this part of the code up here
            temp_array = np.array(X['As_ns_Eas'])*-1
            temp_array = temp_array.reshape((temp_array.shape[0],
                                                      1))
            
            Y = np.vstack((Y, temp_array))
            #clean this part of the code up here
            #tab
            if master_equation_flag == True:
                temp_array = np.array(X['molecular_parameters'])*-1
                temp_array = temp_array.reshape((temp_array.shape[0],
                                                      1))
                Y = np.vstack((Y,temp_array))
                
            
           
    #Assembling the phsycial portion of the Y matrix 
        if exp_dict_list[0]['simulation'].physicalSens ==1:
            if loop_counter == 0:
                for i,exp_dic in enumerate(exp_dict_list):
                    dic_of_conditions = exp_dic['simulation'].conditions
                    #subtract out the dilluant 
                    species_in_simulation = len(set(dic_of_conditions.keys()).difference(['Ar','AR','ar','HE','He','he','Kr','KR','kr','Xe','XE','xe','NE','Ne','ne']))
                    
                    #add two for Temperature and Pressure
                    len_of_phsycial_observables_in_simulation = species_in_simulation + 2 
                    temp_zeros = np.zeros((len_of_phsycial_observables_in_simulation,1))
                    #stacking the zeros onto the Y array 
                    Y = np.vstack((Y,temp_zeros))
                    Y_data_Frame.append('T'+'_'+'experiment'+'_'+str(i))
                    Y_data_Frame.append('P'+'_'+'experiment'+'_'+str(i))
                    for variable in range(species_in_simulation):
                        Y_data_Frame.append('X'+'_'+str(variable)+'_'+'experiment'+'_'+str(i))
            else:
                for i,exp_dic in enumerate(exp_dict_list):
                    dic_of_conditions = exp_dic['simulation'].conditions
                    #subtract out the dilluant 
                    species_in_simulation = len(set(dic_of_conditions.keys()).difference(['Ar','AR','ar','HE','He','he','Kr','KR','kr','Xe','XE','xe','NE','Ne','ne']))
                    
                    Y_data_Frame.append('T'+'_'+'experiment'+'_'+str(i))
                    Y_data_Frame.append('P'+'_'+'experiment'+'_'+str(i))
                    for variable in range(species_in_simulation):
                        Y_data_Frame.append('X'+'_'+str(variable)+'_'+'experiment'+'_'+str(i))
                
                
                
                temp_array = np.array(X['physical_observables'])*-1
                temp_array = temp_array.reshape((temp_array.shape[0],
                                                      1))
                
                Y = np.vstack((Y,temp_array))
            
        #Assembling the portion of the Y matrix for the absorbance coefficient sensitiviteis 
        pert_coef = {} #build a dict matching pert_coef to their experiment and wavelength.             
               #length of the dict gives padding information
        for exp in exp_dict_list:
            if 'perturbed_coef' not in exp.keys():
                continue
            perturbed_for_exp = exp['perturbed_coef']
            for x in perturbed_for_exp:
                if x[0][2] not in pert_coef.keys():
                    pert_coef[x[0][2]] = [x[1]]
                else:
                    pert_coef[x[0][2]].append(x[1])
                    
        num_ind_pert_coef = len(pert_coef)         
        temp_zeros = np.zeros((num_ind_pert_coef,1))
        if loop_counter == 0:
            Y = np.vstack((Y,temp_zeros))
        else:  
            if 'absorbance_coefficent_observables' in X.keys():
                
                #temp_array = np.array(X['absorbance_coefficent_observables'])
                temp_array = X['absorbance_coefficent_observables'] 
                temp_array = [a for a in temp_array if a != 'null']
               
                #temp_array = temp_array[temp_array!=0]
                #temp_array = temp_array[temp_array!=0]
                temp_array = np.array(temp_array)
                temp_array = np.array(temp_array)*-1
                
                temp_array = temp_array.reshape((temp_array.shape[0],
                                                      1))
                
                Y = np.vstack((Y,temp_array))
                
        for x in range(num_ind_pert_coef):
            Y_data_Frame.append('Sigma'+'_'+str(x))

 
        Y_data_Frame = pd.DataFrame({'value': Y_data_Frame,'ln_difference': Y.reshape((Y.shape[0],))})  
        self.Y_matrix = Y
        return Y, Y_data_Frame
        
      

    def build_Z(self, exp_dict_list:list,
                parsed_yaml_file_list:list,
                loop_counter:int = 0,
                reaction_uncertainty=None,
                master_equation_uncertainty_df=None,
                master_equation_reaction_list=[],
                master_equation_flag = False):
        Z = []
        Z_data_Frame = [] 
        sigma = []
        #need to append to sigma
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
        #tab, start working here tomorrow with how we want to read in csv file     
        for i,exp_dic in enumerate(exp_dict_list):
            counter = 0
            for j,observable in enumerate((exp_dic['mole_fraction_observables']+
                                           exp_dic['concentration_observables'])):
                if observable == None:
                    pass
                elif observable in exp_dic['mole_fraction_observables']:
                    ## add ppm statment here ? check if it exists? and add concentration statment below just for parcing 
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp_dic['uncertainty']['mole_fraction_relative_uncertainty'][counter],
                        exp_dic['uncertainty']['mole_fraction_absolute_uncertainty'][counter],
                        exp_dic['experimental_data'][counter][observable].values,exp_dic['experimental_data'][counter])
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                      1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))
                    
                else:                
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp_dic['uncertainty']['concentration_relative_uncertainty'][counter],
                         exp_dic['uncertainty']['concentration_absolute_uncertainty'][counter],
                         exp_dic['experimental_data'][counter][observable+'_ppm'].values,exp_dic['experimental_data'][counter])
                   
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                      1))
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))
                    
                    
                    
                    
                    
                    
                    Z.append(total_uncertainty)
                    sigma.append(un_weighted_uncertainty)
                    tempList = [observable+'_'+'experiment'+str(i)]*np.shape(total_uncertainty)[0]
                    Z_data_Frame.extend(tempList)
                    counter+=1
            if 'absorbance_observables' in list(exp_dic.keys()):
                wavelengths = parsed_yaml_file_list[i]['absorbanceCsvWavelengths']
                    
                for k,wl in enumerate(wavelengths):
                    total_uncertainty,un_weighted_uncertainty = uncertainty_calc(exp_dic['uncertainty']['absorbance_relative_uncertainty'][k],
                                                             exp_dic['uncertainty']['absorbance_absolute_uncertainty'][k],
                                                             exp_dic['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp_dic['absorbance_experimental_data'][k])
                        
                    total_uncertainty = total_uncertainty.reshape((total_uncertainty.shape[0],
                                                  1))
                    
                    un_weighted_uncertainty = un_weighted_uncertainty.reshape((un_weighted_uncertainty.shape[0],
                                                                               1))
                    
                    tempList = [str(wl)+'_'+'experiment'+'_'+str(i)]*np.shape(total_uncertainty)[0]
                    Z_data_Frame.extend(tempList)                   
                    Z.append(total_uncertainty)
                    sigma.append(un_weighted_uncertainty)
                        
        Z = np.vstack((Z))
        
        sigma = np.vstack((sigma))
        
        #Here we are adding A,n,and Ea uncertainty
        #we go do not through an additional step to make sure that the A,N and Ea 
        #values are paired with the correct reactions as in the old code,
        #because we wrote a function to make the excel sheet which will arrange things in the correct order 
        #We also need to decide if we want to put this in as ln values or not in the spreadsheet 
        active_parameters = []
        reaction_uncertainty = pd.read_csv(reaction_uncertainty)
        
        if master_equation_flag ==True:
            for reaction in master_equation_reaction_list:
                index = reaction_uncertainty.loc[reaction_uncertainty['Reaction'] == reaction].index[0]
                reaction_uncertainty = reaction_uncertainty.drop([index])
                
                
                
                

        #tab fix this correctly
        uncertainty_As = reaction_uncertainty['Uncertainty A (unit)'].values 
        uncertainty_As = uncertainty_As.reshape((uncertainty_As.shape[0],
                                                  1))
        
        
        
        #uncertainty_As = np.log(uncertainty_As)
        Z = np.vstack((Z,uncertainty_As))
        sigma = np.vstack((sigma,uncertainty_As))
        for variable in range(uncertainty_As.shape[0]):
            Z_data_Frame.append('A'+'_'+str(variable))
            active_parameters.append('A'+'_'+str(variable))
        
        
        
        uncertainty_ns = reaction_uncertainty['Uncertainty N (unit)'].values
        uncertainty_ns = uncertainty_ns.reshape((uncertainty_ns.shape[0],
                                                  1))
        #print(uncertainty_ns)
        Z = np.vstack((Z,uncertainty_ns))
        sigma = np.vstack((sigma,uncertainty_ns))
        for variable in range(uncertainty_ns.shape[0]):
            Z_data_Frame.append('n'+'_'+str(variable))
            active_parameters.append('n'+'_'+str(variable))
        
        
        uncertainty_Eas = reaction_uncertainty['Uncertainty Ea (unit)'].values
        uncertainty_Eas = uncertainty_Eas.reshape((uncertainty_Eas.shape[0],
                                                  1))
     
        
        Z = np.vstack((Z,uncertainty_Eas))
        
        sigma = np.vstack((sigma,uncertainty_Eas))
        for variable in range(uncertainty_Eas.shape[0]):
            Z_data_Frame.append('Ea'+'_'+str(variable))
            active_parameters.append('Ea'+'_'+str(variable))
        

        
        
        if master_equation_flag ==True:
            master_equation_uncertainty = []
            for col in master_equation_uncertainty_df:
                master_equation_uncertainty.append(list(master_equation_uncertainty_df[col].dropna().values))
                
                
            for i,reaction in enumerate(master_equation_uncertainty):
                for j,uncer in enumerate(reaction):
                    Z_data_Frame.append('R'+'_'+str(i)+'_'+'P'+str(j))
                    active_parameters.append(master_equation_reaction_list[i]+'_P_'+str(j))
                    #active_parameters.append('R'+'_'+str(i)+'_'+'P'+str(j))
            ##check this 
            
            master_equation_uncertainty = [item for sublist in master_equation_uncertainty for item in sublist]
            master_equation_uncertainty = np.array(master_equation_uncertainty)
            master_equation_uncertainty = master_equation_uncertainty.reshape((master_equation_uncertainty.shape[0],
                                                  1))
            Z = np.vstack((Z,master_equation_uncertainty))
            sigma = np.vstack((sigma,master_equation_uncertainty))
            
            
        if exp_dict_list[0]['simulation'].physicalSens ==1:
            for i,exp_dic in enumerate(exp_dict_list):
                experiment_physical_uncertainty = []
                #Temperature Uncertainty 
                experiment_physical_uncertainty.append(exp_dic['uncertainty']['temperature_relative_uncertainty'])
                Z_data_Frame.append('T'+'_'+'experiment'+'_'+str(i))
                active_parameters.append('T'+'_'+'experiment'+'_'+str(i))
                #Pressure Uncertainty
                experiment_physical_uncertainty.append(exp_dic['uncertainty']['pressure_relative_uncertainty'])
                Z_data_Frame.append('P'+'_'+'experiment'+'_'+str(i))
                active_parameters.append('P'+'_'+'experiment'+'_'+str(i))
                #Species Uncertainty
                species_uncertainties = exp_dic['uncertainty']['species_relative_uncertainty']['dictonary_of_values']
                species_to_loop =  exp_dic['uncertainty']['species_relative_uncertainty']['species']
                dilluant = ['Ar','AR','ar','HE','He','he','Kr','KR','kr','Xe','XE','xe','NE','Ne','ne']
                for specie in species_to_loop:
                    if specie in dilluant:
                        continue
                    experiment_physical_uncertainty.append(species_uncertainties[specie])
                    Z_data_Frame.append('X'+'_'+str(specie)+'_'+'experiment'+'_'+str(i))
                    active_parameters.append('X'+'_'+str(specie)+'_'+'experiment'+'_'+str(i))

                experiment_physical_uncertainty = np.array(experiment_physical_uncertainty)
                experiment_physical_uncertainty =  experiment_physical_uncertainty.reshape((experiment_physical_uncertainty.shape[0],
                                                  1))
                Z = np.vstack((Z,experiment_physical_uncertainty))
                sigma = np.vstack((sigma,experiment_physical_uncertainty))
        #building dictonary to keep track of independtend coupled coefficients 
        count = 0
        coef_dict = {} 
        
        uncertainties_of_coefficents = []
        for i,exp_dic in enumerate(exp_dict_list):
            if 'perturbed_coef' not in exp_dic.keys():
                continue
            dictonary_of_coef_and_uncertainty = exp_dic['uncertainty']['coupled_coef_and_uncertainty']
            
            for x in dictonary_of_coef_and_uncertainty:
                if x not in coef_dict.keys():
                    coef_dict[x] = dictonary_of_coef_and_uncertainty[x]
                    
        for x in coef_dict:       
            for y in coef_dict[x]:
                if y[0]!=0:                            #this might cause a problem in the future
                        count+=1
                        uncertainties_of_coefficents.append(y)
                        Z_data_Frame.append('Sigma'+'_'+str(count))
                        active_parameters.append('Sigma'+'_'+str(count))
                        
        uncertainties_of_coefficents = np.array(uncertainties_of_coefficents)
        
       
        if uncertainties_of_coefficents.any() == True:
            uncertainties_of_coefficents =  uncertainties_of_coefficents.reshape((uncertainties_of_coefficents.shape[0],
                                                  1))
            Z = np.vstack((Z,uncertainties_of_coefficents))            
            sigma = np.vstack((sigma,uncertainties_of_coefficents))
        
        Z_data_Frame = pd.DataFrame({'value': Z_data_Frame,'Uncertainty': Z.reshape((Z.shape[0],))})       
        self.z_matrix = Z
        self.sigma = sigma
        return Z,Z_data_Frame,sigma,active_parameters
    
    
    
    
    
    def breakup_X(self, X, 
                        exp_dict_list:list, 
                        exp_uncertainty_dict_list_original:list,
                        loop_counter:int = 0,
                        master_equation_uncertainty_df=None,
                        master_equation_reactions = [],
                        master_equation_flag = False):
        
        
        X_to_subtract_from_Y = {}
        reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
        number_of_reactions = len(reactions_in_cti_file)
        
        ####Grab off updates directly for the CTI file 
        ####need to add master equation reactions 
        

        
        ##################################################################
        if loop_counter !=0:
            X_new = X 
        
            
        else:
            X_new = X
        ##################################################################
        #print('USING BURKE X VALUES')
        #X = pd.read_csv('MSI/data/test_data/burke_X_values.csv')
        #X= X['Burke_Value'].values
        #X = X.reshape(X.shape[0],1)
        
        ################################################################
        
        
        X_new = list(X_new.flatten())            
        if exp_dict_list[0]['simulation'].kineticSens ==1:
            
            value1 = 3*(number_of_reactions - len(master_equation_reactions))
            
               
            AsNsEas = X_new[:value1]
            X_to_subtract_from_Y['As_ns_Eas'] = AsNsEas
            
            #### pickup here
            dividedBy = int(len(AsNsEas) / 3)
            
            def list_slice(S,step):
                return [S[i::step] for i in range(step)]
            
            resortedList = list_slice(AsNsEas,dividedBy)
            
            innerDict = ['A','n','Ea']
            l = [dict(zip(innerDict,resortedList[x])) for x in range(len(resortedList))]    
            Keys= []
            for xx in range(int(value1/3)):
                Keys.append('r'+str(xx))
                
                
            deltaXAsNsEas = dict(zip(Keys,l))
            
            
            innerDictNew = ['A_update','n_update','Ea_update']
            ll = [dict(zip(innerDictNew,resortedList[x])) for x in range(len(resortedList))]
            kinetic_paramter_dict = dict(zip(reactions_in_cti_file,ll))
        #molecularParams = np.array([.1,.2,.3,.4,.2,.3,.4]).flatten().tolist()
        # might need to fix this based on how lei is passing me information, check in notebook
        
        if master_equation_flag == True:
            number_of_molecular_parameters_list = []
            for col in master_equation_uncertainty_df:
                number_of_molecular_parameters_list.append(len(master_equation_uncertainty_df[col].dropna().values))
                
            sum_of_moleular_paramters = sum(number_of_molecular_parameters_list)
            value2 = sum_of_moleular_paramters     
            deltaXmolecularParams = X_new[value1:(value1+value2)]
            
            X_to_subtract_from_Y['molecular_parameters'] = deltaXmolecularParams
            molecular_paramters_by_reaction = []
            reaction_numbers = []
            start_mp = 0
            for r,number in enumerate(number_of_molecular_parameters_list):
                stop_mp = start_mp + number
                molecular_paramters_by_reaction.append(deltaXmolecularParams[start_mp:stop_mp])
                start_mp = stop_mp
                reaction_numbers.append('R_'+str(r))

            delta_x_molecular_params_by_reaction_dict = dict(zip(master_equation_reactions,molecular_paramters_by_reaction))
            list_of_mp = []
            for i,reaction in enumerate(molecular_paramters_by_reaction):
                temp=[]
                for j,value in enumerate(reaction):
                    temp.append('Paramter_'+str(j)+'_Update')
                list_of_mp.append(temp)
            inner_dict_temp = [dict(zip(list_of_mp[x],molecular_paramters_by_reaction[x])) for x in range(len(molecular_paramters_by_reaction))]
            inner_dict_temp_2 = dict(zip(master_equation_reactions,inner_dict_temp))
                    
            kinetic_paramter_dict.update(inner_dict_temp_2)

            
           
       
        else:
            value2 = 0
        
        
        physical_observables = []
        previous_value = 0
        physical_observables_for_Y = []
        if exp_dict_list[0]['simulation'].physicalSens ==1:
            for i,exp_dic in enumerate(exp_dict_list):
                dic_of_conditions = exp_dic['simulation'].conditions
                    #subtract out the dilluant 
                species_in_simulation = len(set(dic_of_conditions.keys()).difference(['Ar','AR','ar','HE','He','he','Kr','KR','kr','Xe','XE','xe','NE','Ne','ne']))
                    #add two for Temperature and Pressure
                len_of_phsycial_observables_in_simulation = species_in_simulation + 2 
                new_value = previous_value + len_of_phsycial_observables_in_simulation
                single_experiment_physical_observables = X_new[(value1+value2+previous_value):(value1+value2+new_value)]
                physical_observables_for_Y.append(single_experiment_physical_observables)
                temp_keys = []
                    #stacking the zeros onto the Y array 
                temp_keys.append('T'+'_'+'experiment'+'_'+str(i))
                temp_keys.append('P'+'_'+'experiment'+'_'+str(i))
                for variable in range(species_in_simulation):
                    temp_keys.append('X'+'_'+str(variable)+'_'+'experiment'+'_'+str(i))
                temp_dict = dict(zip(temp_keys,single_experiment_physical_observables))
                physical_observables.append(temp_dict)
                ##come back to this and do a test on paper
                previous_value = new_value
                
        physical_observables_for_Y = [item for sublist in physical_observables_for_Y for item in sublist]   
        X_to_subtract_from_Y['physical_observables'] = physical_observables_for_Y
        
        
        test_abs = []
        absorbance_coefficients_for_Y = []
        coef_dict = {}  
        coef_dict_list = []
        absorbance_coef_update_dict = {}
        for i,exp_dic in enumerate(exp_uncertainty_dict_list_original):
            if 'coupled_coef_and_uncertainty' not in exp_dic.keys():
                continue
            dictonary_of_coef_and_uncertainty = exp_dic['coupled_coef_and_uncertainty']
            #tab start working here tomorrow, need to pass in the original version of this dict 
            #dictonary_of_coef_and_uncertainty = {(140000, 0.0): ([0.7], [0.0]), (1270000, 0.0): ([0.7], [0.0])}

            for x in dictonary_of_coef_and_uncertainty:
                if x not in coef_dict.keys():
                    coef_dict[x] = dictonary_of_coef_and_uncertainty[x]
                if x not in coef_dict_list:
                    coef_dict_list.append(x)
                    
                    
                    
        start_abs = 0
        stop_abs = 1         
        for i,cof in enumerate(coef_dict_list):
            temp=[]
            temp2=[]
         #   counter=1
            for value in cof:
                if value==0:
                    temp.append([0])
                    temp2.append(['null'])
                else:
                    temp.append(X_new[(value1+value2+new_value+start_abs):(value1+value2+new_value+stop_abs)])
                    temp2.append(X_new[(value1+value2+new_value+start_abs):(value1+value2+new_value+stop_abs)])
                    start_abs = stop_abs
                    stop_abs +=1
                    
            temp = [item for sublist in temp for item in sublist]     
            temp2 = [item for sublist in temp2 for item in sublist]     
            absorbance_coef_update_dict[cof] = temp
            absorbance_coefficients_for_Y.append(temp2)
            test_abs.append(temp2)
            
            

        # return everything in a dictonary??   
        absorbance_coefficients_for_Y = [item for sublist in absorbance_coefficients_for_Y for item in sublist] 
        X_to_subtract_from_Y['absorbance_coefficent_observables'] = absorbance_coefficients_for_Y
#       
        if master_equation_flag == False:
            return deltaXAsNsEas,physical_observables,absorbance_coef_update_dict,X_to_subtract_from_Y,kinetic_paramter_dict
        else:
            return deltaXAsNsEas,physical_observables,absorbance_coef_update_dict,X_to_subtract_from_Y,delta_x_molecular_params_by_reaction_dict,kinetic_paramter_dict
    
    
    def matrix_manipulation(self,runCounter,S_matrix,Y_matrix,z_matrix,XLastItteration = np.array(()),active_parameters=[]):

    

        one_over_z = np.true_divide(1,z_matrix)
        y_matrix = Y_matrix * one_over_z
        s_matrix = S_matrix * (one_over_z.flatten()[:,np.newaxis])
        self.y_matrix = y_matrix
        
        sTimesZ = S_matrix * (z_matrix.flatten())[:,np.newaxis]
        #calculate covariance matrix 
        shape = np.shape(self.S_matrix_wo_k_targets)
        #print(shape,'shape of s matrix wo k targets')
        #print(shape[0]-len(active_parameters))
        
        #print(S_matrix.shape,'shape of s matrix')
        
        s_wo_k_targets = s_matrix[:shape[0],:shape[1]]
        identity_matrix = s_wo_k_targets[shape[0]-len(active_parameters):,:]
        #print(identity_matrix.shape)
        #print(len(active_parameters),'number of active parameters')
        
        
        try:
            if runCounter==0:
              
                c = np.dot(np.transpose(identity_matrix),identity_matrix)
                c = np.linalg.inv(c)
                prior_diag = np.diag(c)
                prior_sigmas = np.sqrt(prior_diag)
                covariance_prior_df = pd.DataFrame(c)
                covariance_prior_df.columns = active_parameters
                covariance_prior_df.reindex(labels = active_parameters)
                prior_diag_df = pd.DataFrame({'parameter': active_parameters,'value': prior_diag.reshape((prior_diag.shape[0],))})
                sorted_prior_diag = prior_diag_df.sort_values(by=['value'])
                prior_sigmas_df = pd.DataFrame({'parameter': active_parameters,'value': prior_sigmas.reshape((prior_sigmas.shape[0],))})
            else:
                c = np.dot(np.transpose(s_matrix),s_matrix)
                c = np.linalg.inv(c)
                
                covariance_posterior_df = pd.DataFrame(c)
                covariance_posterior_df.columns = active_parameters
                covariance_posterior_df.reindex(labels = active_parameters)
                posterior_diag = np.diag(c)
                posterior_sigmas = np.sqrt(posterior_diag)
                posterior_sigmas_df = pd.DataFrame({'parameter': active_parameters,'value': posterior_sigmas.reshape((posterior_sigmas.shape[0],))})
                posterior_diag_df =  pd.DataFrame({'parameter': active_parameters,'value': posterior_diag.reshape((posterior_diag.shape[0],))})
                sorted_posterior_diag  = posterior_diag_df.sort_values(by=['value'])
        except:
            #stub
            print('WE ARE IN THE EXCEPT STATMENT')
            if runCounter==0:
              
                c = -1
                c = -1
                prior_diag = -1
                prior_sigmas = -1
                covariance_prior_df = -1
                prior_diag_df = -1
                sorted_prior_diag = -1
                prior_sigmas_df = -1
            else:
                c = -1
                c =-1
                
                covariance_posterior_df = -1
                posterior_diag = -1
                posterior_sigmas = -1
                posterior_sigmas_df = -1
                posterior_diag_df =  -1
                sorted_posterior_diag  = -1
        
        
        self.covariance = c


        
    
        self.s_matrix = s_matrix
    
        psudoInverse = np.linalg.pinv(s_matrix)
        delta_X = np.dot(psudoInverse,y_matrix)
        self.delta_X = delta_X
        
    

        if runCounter == 0:
            XlastItteration = np.zeros(np.shape(delta_X))

        else:     
            XlastItteration = XLastItteration

        X = XlastItteration + delta_X

        self.X = X
       
        #STUB THIS
        try:
            X_data_frame = pd.DataFrame({'value': active_parameters,'Parameter': X.reshape((X.shape[0],))})
        except:
            X_data_frame = -1
        if runCounter==0:
            return X,c,s_matrix,y_matrix,delta_X,z_matrix,X_data_frame,prior_diag,prior_diag_df,sorted_prior_diag,covariance_prior_df,prior_sigmas_df
        else:
            return X,c,s_matrix,y_matrix,delta_X,z_matrix,X_data_frame,posterior_diag,posterior_diag_df,sorted_posterior_diag,covariance_posterior_df,posterior_sigmas_df
            
    
    
    
        
class Adding_Target_Values(meq.Master_Equation):
    def __init__(self,S_matrix,Y_matrix,z_matrix,sigma,Y_data_Frame,z_data_Frame):
        self.S_matrix = S_matrix
        self.Y_matrix = Y_matrix
        self.z_matrix = z_matrix
        self.sigma = sigma
        self.Y_data_Frame = Y_data_Frame
        self.z_data_Frame = z_data_Frame
        meq.Master_Equation.__init__(self)
        
         
        
    def target_values_Y(self,target_value_csv,exp_dict_list:list):
        import cantera as ct
        Y_df_list = []
        Y_values = []
        #make sure we put the reactions into the file in the units cantera uses
        target_value_csv = pd.read_csv(target_value_csv)
        target_reactions = target_value_csv['Reaction']
        target_temp = target_value_csv['temperature']
        target_press = target_value_csv['pressure']
        target_k = target_value_csv['k']
        bath_gas = target_value_csv['M']
        reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
        gas = ct.Solution(exp_dict_list[0]['simulation'].processor.cti_path)
        diff_in_ks_for_Y = []
        for i,reaction in enumerate(target_reactions): 
                #ask about the mixture composition
            if target_press[i] == 0:
                pressure = 1e-9
            else:
                pressure = target_press[i]
            if bath_gas[i] !=0:
                gas.TPX = target_temp[i],pressure*101325,{'H2O':.013,'O2':.0099,'H':.0000007,'Ar':.9770993}

            else:
                gas.TPX = target_temp[i],pressure*101325,{'Ar':.99}
            reaction_number_in_cti = reactions_in_cti_file.index(reaction)
            k = gas.forward_rate_constants[reaction_number_in_cti]
            
                #check and make sure we are subtracting in the correct order 
            difference = np.log(target_k[i]) - np.log(k) 

            diff_in_ks_for_Y.append(difference)
            Y_df_list.append(reaction)
            Y_values.append(difference)
            
        k_targets_for_y = np.array(diff_in_ks_for_Y)
        k_targets_for_y = k_targets_for_y.reshape((k_targets_for_y.shape[0],1))
        Y_values = np.array(Y_values)
        
        Y_df_temp = pd.DataFrame({'value': Y_df_list,'ln_difference': Y_values.reshape((Y_values.shape[0],))}) 
        self.Y_data_Frame = self.Y_data_Frame.append(Y_df_temp, ignore_index=True)
        
        return k_targets_for_y,self.Y_data_Frame
    
    def target_values_for_Z(self,target_value_csv):
        z_over_w = []
        sigma = []
        target_value_csv = pd.read_csv(target_value_csv)
        target_ln_uncertainty = target_value_csv['ln_unc_k']
        target_W = target_value_csv['W']
        target_reactions = target_value_csv['Reaction']
        z_df_list=[]
        z_values = []
        for i,value in enumerate(target_ln_uncertainty):
            temp = np.divide(value,target_W[i])
            sigma.append(value)
            z_over_w.append(temp)
            z_values.append(temp)
            z_df_list.append(target_reactions[i])
            
        k_targets_for_z = np.array(z_over_w)
        sigma = np.array(sigma)
        sigma = sigma.reshape((sigma.shape[0],1))
        z_values = np.array(z_values)
        k_targets_for_z = k_targets_for_z.reshape((k_targets_for_z.shape[0],1))
        Z_data_Frame_temp = pd.DataFrame({'value': z_df_list,'Uncertainty': z_values.reshape((z_values.shape[0],))})
        self.z_data_Frame = self.z_data_Frame.append(Z_data_Frame_temp, ignore_index=True)    
        return k_targets_for_z,sigma,self.z_data_Frame
    




    
    def target_values_for_S(self,target_value_csv,
                            exp_dict_list,
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
            
        Number_of_MP = []
        nested_reaction_list = [[] for x in range(len(master_equation_reaction_list))]
            
        for reaction in master_equation_reaction_list:
            for MP in master_equation_reaction_list[reaction]:
                nested_reaction_list[reaction].append(0)
                Number_of_MP.append(MP)
        copy.deepcopy(nested_reaction_list)            
        Number_of_MP = len(Number_of_MP)
              
        MP_stack = []
        target_values_to_stack =  []
        for i,reaction in enumerate(target_reactions):
            #temp_array = np.zeros((1,Number_of_MP))
            if reaction in master_equation_reaction_list:
                indx = master_equation_reaction_list.index(reaction)
                MP_sens_array_list = master_equation_sensitivites[reaction]
                copy.deepcopy(MP_sens_array_list)
                MP_sens_array_list_copy = MP_sens_array_list
                for j, MP_array in MP_sens_array_list:
                    #alpha_array = np.zeros(MP_array.shape)
                    for sensitivity in np.nditer(MP_array,order='F'):
                        k,l= np.where(MP_array == sensitivity)
                            #need to add reduced p and t, and check these units were using to map
                            
                            #these might not work
                        t_alpha= meq.Master_Equation.chebyshev_specific_poly(k[0],meq.Master_Equation.calc_reduced_T(target_temp[i]))
                        p_alpha = meq.Master_Equation.chebyshev_specific_poly(l[0],meq.Master_Equation.calc_reduced_P(target_press[i]))
                        
                        #these might nowt work 
                        alpha = t_alpha*p_alpha
                        MP_sens_array_list_copy[j][k,l] = alpha
                        
                        #needt to figure out how to put this in the correct spot un the array of zeros 
                    multiplied = np.multiply(MP_sens_array_list_copy[j],MP_sens_array_list[j])
                    array_sum = sum(multiplied)
                    temp = nested_reaction_list
                    temp[indx][j] = array_sum
                    MP_stack.append(temp)
                    flat_list = [item for sublist in temp for item in sublist]
                    flat_list = np.array(flat_list)
                    flat_list = flat_list.reshape((1,flat_list.shape[1]))                                                         
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
                
                
            
        S_matrix = self.S_matrix
        shape_s = S_matrix.shape
        S_target_values = []
        for i,row in enumerate(target_values_to_stack):
            if target_reactions[i] in master_equation_reaction_list:
                zero_to_append_infront = np.zeros((1,(number_of_reactions_in_cti-len(master_equation_reaction_list)*3)))
                zero_to_append_behind = np.zeros((1, shape_s[1] - (number_of_reactions_in_cti-len(master_equation_reaction_list)*3) - np.shape(row)[1] ))                
                temp_array = np.hstack((zero_to_append_infront,row,zero_to_append_behind))
                S_target_values.append(temp_array)
            else:
                zero_to_append_behind = np.zeros((1,shape_s[1]-np.shape(row)[1]))
                temp_array = np.hstack((row,zero_to_append_behind))
                S_target_values.append(temp_array)
                
                
                
            
            
        
        S_target_values = np.vstack((S_target_values))    
        return  S_target_values
    
    
    def appending_target_values(self,target_values_for_z,
                                target_values_for_y,
                                target_values_for_s,
                                sigma_target_values):
        z_matrix = np.vstack((self.z_matrix ,target_values_for_z))
        Y_matrix = np.vstack((self.Y_matrix,target_values_for_y))
        
        S_matrix = np.vstack((self.S_matrix,target_values_for_s))
        sigma = np.vstack((self.sigma,sigma_target_values))
        
        self.S_matrix = S_matrix
        self.Y_matrix = Y_matrix
        self.z_matrix = z_matrix
        self.sigma = sigma
        
        return S_matrix,Y_matrix,z_matrix,sigma
    

        

        
        
        
        
                
        
            

            
                
                
            
                
                
       
            
            
                        
                        
                        
               
            
                
                
            
            
        
        
        
              
