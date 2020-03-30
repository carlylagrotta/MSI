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
        
#        identity_matrix[1,0]=.1
#        identity_matrix[0,1]=.1
#        identity_matrix[0,20]=.1
#        identity_matrix[20,0]=.1
#        identity_matrix[39,0]=.1
#        identity_matrix[0,39]=.1
        
        ####making edits to this just for masten test 
        
        
        
        S_matrix = np.vstack((S_matrix,identity_matrix))
        self.S_matrix = S_matrix
        S_matrix_wo_k_targets = copy.deepcopy(self.S_matrix)
        self.S_matrix_wo_k_targets = S_matrix_wo_k_targets
        S_matrix_df = pd.DataFrame(S_matrix)
        return S_matrix 