import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cantera as ct
import copy
from textwrap import wrap

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
                 simulation_run=None):
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
                    length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Time'].shape[0])
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Time'].shape[0])
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
                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable],'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original['simulation'].timeHistories[0][observable],'r',label= "$\it{A priori}$ model")
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable],'o',color='black',label='Experimental Data')
                    plt.xlabel('Time (ms)')
                    plt.ylabel('Mole Fraction'+''+str(observable))
                    plt.title('Experiment_'+str(i+1))
                    
                    
                    
                    

                    
                    if bool(sigmas_optimized) == True:
                        
                        high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                        high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                        low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        plt.figure()
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,'b--')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,'b--')
                        
                        
                        
                        high_error_original = np.exp(sigmas_original[i][observable_counter])
                        high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        low_error_original = np.exp(sigmas_original[i][observable_counter]*-1)
                        low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        plt.figure()
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'r--')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'r--')
                    
                    plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight',dpi=1000)
                    

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

                    plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                    plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)

                    
                    

# make function to plot rate constants 
                    
    def plotting_rate_constants(self,optimized_cti_file='',
                                original_cti_file='',
                                initial_temperature=250,
                                final_temperature=2500):
        
        gas_optimized = ct.Solution(optimized_cti_file)
        gas_original = ct.Solution(original_cti_file)

        def unique_list(seq):
            checked = []
            for e in seq:
                if e not in checked:
                    checked.append(e)
            return checked
        
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
                k.append(gas.forward_rate_constants[reaction_number])
            return Temp,k

        def calculate_sigmas_for_rate_constants(k_target_value_S_matrix,k_target_values_parsed_csv,unique_reactions,gas,covarience):

            
            reaction_list_from_mechanism = gas.reaction_equations()
            sigma_list_for_target_ks = [[] for reaction in range(len(unique_reactions))]
            shape = k_target_value_S_matrix.shape
            for row in range(shape[0]):
                SC = np.dot(k_target_value_S_matrix[row,:],covarience)
                sigma_k = np.dot(SC,np.transpose(k_target_value_S_matrix[row,:]))
                sigma_k = np.sqrt(sigma_k)
                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
                
            return sigma_list_for_target_ks
        
        def calculating_target_value_ks_from_cantera_for_sigmas(k_target_values_parsed_csv,gas,unique_reactions):
            target_value_ks = [[] for reaction in range(len(unique_reactions))]
            
            
            
            
            target_reactions = k_target_values_parsed_csv['Reaction']
            target_temp = k_target_values_parsed_csv['temperature']
            target_press = k_target_values_parsed_csv['pressure']
            reactions_in_cti_file = gas.reaction_equations()
            for i,reaction in enumerate(target_reactions): 
                #ask about the mixture composition
                if target_press[i] == 0:
                    pressure = 1e-9
                else:
                    pressure = target_press[i]
                    
                gas.TPX = target_temp[i],pressure*101325,{'Ar':.99}
                reaction_number_in_cti = reactions_in_cti_file.index(reaction)
                k = gas.forward_rate_constants[reaction_number_in_cti]
                indx = reaction_list_from_mechanism.index(reaction)
                target_value_ks[unique_reactions.index(indx)].append(k)
                
                
                
                #check and make sure we are subtracting in the correct order 

            return target_value_ks
    
    
        if bool(self.target_value_rate_constant_csv) and self.k_target_values=='On':
            
            
            #make two unique
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
            
            
            

            sigma_list_for_target_ks_optimized = calculate_sigmas_for_rate_constants(self.k_target_value_S_matrix,k_target_value_csv,unique_reactions_optimized,gas_optimized,self.covarience)
            self.sigma_list_for_target_ks_optimized = sigma_list_for_target_ks_optimized
            sigma_list_for_target_ks_original = calculate_sigmas_for_rate_constants(self.k_target_value_S_matrix,k_target_value_csv,unique_reactions_original,gas_original,self.original_covariance)
            self.sigma_list_for_target_ks_original = sigma_list_for_target_ks_original
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
                #calculate sigmas 
                high_error_optimized = np.exp(np.array(sigma_list_for_target_ks_optimized[i]))
                high_error_optimized = np.multiply(high_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])
                
                
                low_error_optimized = np.exp(np.array(sigma_list_for_target_ks_optimized[i])*-1)
                low_error_optimized = np.multiply(low_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])    
                
                #plt.semilogy(target_value_temps_optimized[i],high_error_optimized,'b--')   
                
                a, b = zip(*sorted(zip(target_value_temps_optimized[i],high_error_optimized)))

                plt.scatter(a,b,color='blue')
               # print(a,b)
                #plt.semilogy(a,b,'b--')
                
                a, b = zip(*sorted(zip(target_value_temps_optimized[i],low_error_optimized)))  
                #plt.semilogy(a,b,'b--')
                plt.scatter(a,b,color='blue')
               # print(a,b)
                Temp_original,k_original = rate_constant_over_temperature_range_from_cantera(unique_reactions_original[unique_reactions_original.index(reaction)],
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_original,k_original,'r')
                
                high_error_original = np.exp(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction)])
                high_error_original = np.multiply(high_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction)])
                
               
                low_error_original = np.exp(np.array(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction)])*-1)
                low_error_original = np.multiply(low_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction)])  
                
                a, b = zip(*sorted(zip(target_value_temps_original[unique_reactions_original.index(reaction)],high_error_original)))  
                #plt.semilogy(a,b,'r--')
                plt.scatter(a,b,color='red')
                
                
                a, b = zip(*sorted(zip(target_value_temps_original[unique_reactions_original.index(reaction)],low_error_original)))  
                #plt.semilogy(a,b,'r--')
                plt.scatter(a,b,color='red')
                
                plt.semilogy(target_value_temps_optimized[i],target_value_ks_optimized[i],'o',color='black')
                
                plt.xlabel('Temperature (K)')
                plt.ylabel('Kmol/m^3-s')
                plt.title(reaction_list_from_mechanism[reaction])
                
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
                ccn=pd.DataFrame(cc).as_matrix()
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
                    time_profiles[i].append(exp['experimental_data'][observable_counter]['Time']*1e3)
                    observables[i].append(observable)
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    time_profiles[i].append(exp['experimental_data'][observable_counter]['Time']*1e3)        
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
                
                plt.xlabel('Temperature (K)')
                plt.ylabel('Kmol/m^3-s')
                plt.title(reaction_list_from_mechanism[reaction])
                plt.savefig(self.working_directory+'/'+reaction_list_from_mechanism[reaction]+'.pdf', bbox_inches='tight')
                plt.savefig(self.working_directory+'/'+reaction_list_from_mechanism[reaction]+'.svg', bbox_inches='tight')

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

                
                
            
        
                    
                    
                    
                    
                    
                    
                    
                    
