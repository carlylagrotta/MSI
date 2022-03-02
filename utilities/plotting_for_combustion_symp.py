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
import os 
import re








class Plotting_for_2020_symposium(object):
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
                  mapped_to_alpha_full_simulation=None,
                  optimized_cti_file='',
                  original_cti_file='',
                  sigma_ones=False,
                  T_min=200,
                  T_max=3000,
                  P_min=1013.25,
                  P_max=1.013e+6,
                  T_P_min_max_dict={}):
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
        self.mapped_to_alpha_full_simulation = mapped_to_alpha_full_simulation,
        self.new_cti=optimized_cti_file
        self.nominal_cti=original_cti_file
        self.sigma_ones = sigma_ones
        self.T_min = T_min
        self.T_max = T_max
        self.P_min = P_min
        self.P_max = P_max
        self.T_P_min_max_dict = T_P_min_max_dict
        
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
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original['simulation'].timeHistories[0][observable],'r',label= "$\it{a}$ $\it{priori}$ model")
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
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable]*1e6,'r',label= "$\it{a}$ $\it{priori}$ model")
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Hong et al. Experimental') 
                    plt.xlabel('Time (ms)')
                    if observable =='H2O':
                        plt.ylabel(r'H$_2$O (ppm)')
                    if observable == 'OH':
                        plt.ylabel('OH (ppm)')
                    else:                        
                        plt.ylabel(r'H$_2$O (ppm)')
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
                        
                        #plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'r--')
                        #plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'r--')
                    
                    #plt.plot([],'w' ,label= 'T:'+ str(self.exp_dict_list_original[i]['simulation'].temperature))
                    #plt.plot([],'w', label= 'P:'+ str(self.exp_dict_list_original[i]['simulation'].pressure))
                    key_list = []
                    for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():
                        
                        #plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                        key_list.append(key)
                   
                    #plt.legend(handlelength=3)
                    plt.legend(ncol=1)
                    sp = '_'.join(key_list)
                    #print(sp)
                    #plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K'+'_'+str(self.exp_dict_list_original[i]['simulation'].pressure)+'_'+sp+'_'+'.pdf', bbox_inches='tight')
                    
                    #stub
                    plt.tick_params(direction='in')
                    
                    
                    plt.savefig('Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                    plt.savefig('Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)
                    #plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight',dpi=1000)



                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                plt.figure()
                for k,wl in enumerate(wavelengths):
                    if wl == 227:
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Kircher et al. Experimental')
                    if wl == 215:
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Kappel et al. Experimental')

                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['absorbance_calculated_from_model'][wl],'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['absorbance_calculated_from_model'][wl],'r',label= "$\it{a}$ $\it{priori}$ model")
                    #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')
                    plt.xlabel('Time (ms)')
                    plt.ylabel('Absorbance'+' ('+str(wl)+' nm)')
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
                    plt.legend(ncol=1)
                    #plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                    plt.tick_params(direction='in')
                    
                    
                    plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                    plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)

    def plotting_3_figre_observables(self,experiment_number_want_to_plot,sigmas_original=[],sigmas_optimized=[]):
        
        OH_x_lim_low = -0.005
        OH_x_lim_high = 1.2
        OH_y_lim_low = 0
        OH_y_lim_high= 270
        
        H2O_x_lim_low = -0.005
        H2O_x_lim_high = 2.5
        H2O_y_lim_low = 1200
        H2O_y_lim_high= 4700
        
        Abs_x_lim_low = -.0002
        Abs_x_lim_high= .125
        Abs_y_lim_low= 0
        Abs_y_lim_high=.25
        
        axis_number_size = 14
        axis_label_size=16
        
        
        
        for i,exp in enumerate(self.exp_dict_list_optimized):
            if i==experiment_number_want_to_plot:
                observable_counter=0
                plt.figure()
                
                #axes=plt.gca()
                #axes.set_aspect(10)

                fig = plt.figure(figsize=(6, 7))
                
                
                for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                    if observable == None:
                        continue
                    

                        
                    if observable in exp['concentration_observables']:
                        if observable == 'H2O':
                            plt.subplot(3,1,1)
                            plt.xlim(H2O_x_lim_low,H2O_x_lim_high)
                            plt.ylim(H2O_y_lim_low,H2O_y_lim_high)
                            plt.xticks(fontsize= axis_number_size)
                            plt.yticks(fontsize= axis_number_size)
                        if observable =='OH':
                            plt.subplot(3,1,2)
                            plt.xlim(OH_x_lim_low,OH_x_lim_high)
                            plt.ylim(OH_y_lim_low,OH_y_lim_high)
                            plt.xticks(fontsize= axis_number_size)
                            plt.yticks(fontsize= axis_number_size)

                        plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable]*1e6,'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable]*1e6,'r',label= "$\it{a}$ $\it{priori}$ model")
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Hong et al. experiment') 
                       # plt.xlabel('Time (ms)')
                        
                        if observable == 'H2O':
                            plt.ylabel(r'H$_2$O [ppm]', fontsize=axis_label_size)
                        if observable == 'OH':
                            plt.ylabel('OH [ppm]',  fontsize=axis_label_size)
                            
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
                    plt.xlim(Abs_x_lim_low, Abs_x_lim_high)
                    plt.ylim(Abs_y_lim_low,Abs_y_lim_high)
                    plt.xticks(fontsize= axis_number_size)
                    plt.yticks(fontsize= axis_number_size)

                    for k,wl in enumerate(wavelengths):
                        if wl == 227:
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Hong et al. experiment')
                        if wl == 215:
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Kappel et al. experiment')
    
                        plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['absorbance_calculated_from_model'][wl],'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['absorbance_calculated_from_model'][wl],'r',label= "$\it{a}$ $\it{priori}$ model")
                        #plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')
                        plt.xlabel('Time [ms]',fontsize=axis_label_size)
                       # plt.ylabel('Absorbance'+' '+str(wl)+' nm',  fontsize=16)
                        plt.ylabel('Absorbance',fontsize=axis_label_size)
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
                        plt.legend(ncol=1,frameon=False,prop={'size': 10.5},loc="best")

                        plt.tick_params(direction='in')

                        plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                        plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight')

                        
                        #plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                        #plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)
                        
                        
                        
                        
                        
    def plotting_1_figre_observables(self,experiment_number_want_to_plot,sigmas_original=[],sigmas_optimized=[]):
        
       #kappel low temp 
        Abs_x_lim_low = -.0002
        Abs_x_lim_high= 2
        Abs_y_lim_low= .2
        Abs_y_lim_high=.5
        
        #kappel medium temp
        Abs_x_lim_low = -.0002
        Abs_x_lim_high= 2
        Abs_y_lim_low= 0
        Abs_y_lim_high=.35


        #kappel high temp
        Abs_x_lim_low = -.0002
        Abs_x_lim_high= .4
        Abs_y_lim_low= 0
        Abs_y_lim_high=.25
        
        #kircher low temp low pressure
        Abs_x_lim_low = -.0002
        Abs_x_lim_high= 5
        Abs_y_lim_low= 0
        Abs_y_lim_high=.15 
        #.1       
        
        # # #kircher high temp low pressure
        # Abs_x_lim_low = -.0002
        # Abs_x_lim_high= 22
        # Abs_y_lim_low= 0
        # Abs_y_lim_high=.125          
        
        # # #kircher low temp high pressure
        # Abs_x_lim_low = -.0002
        # Abs_x_lim_high= 22
        # Abs_y_lim_low= 0
        # Abs_y_lim_high=.1 

        # # #kircher high temp high pressure
        # Abs_x_lim_low = -.0002
        # Abs_x_lim_high= 22
        # Abs_y_lim_low= 0
        # Abs_y_lim_high=.125 
        
        # #kircher high temp high pressure
        # Abs_x_lim_low = -.0002
        # Abs_x_lim_high= 22
        # Abs_y_lim_low= 0
        # Abs_y_lim_high=.1         
        
        
        axis_number_size = 14
        axis_label_size=16
        
        df_kircher = pd.read_csv('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/lei_results_testing/Kircher_Troubleshoot/test2.csv')
        #df_kircher = pd.read_csv('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/lei_results_testing/Kircher_Troubleshoot/test3.csv')
        #df_kircher = pd.read_csv('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/lei_results_testing/Kircher_Troubleshoot/test4.csv')
        #df_kircher = pd.read_csv('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/lei_results_testing/Kircher_Troubleshoot/test5.csv')

        
        
        
        
        
        for i,exp in enumerate(self.exp_dict_list_optimized):
            
            if i==experiment_number_want_to_plot:
                print(i)
                observable_counter=0
                plt.figure()
                
                #axes=plt.gca()
                #axes.set_aspect(10)
                
                fig = plt.figure(figsize=(6, 7))
                                        
    
                if 'perturbed_coef' in exp.keys():
                    wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                    #plt.subplot(3,1,3)
                    
                    plt.subplot(3,1,1)
                    
                    
                    plt.xlim(Abs_x_lim_low, Abs_x_lim_high)
                    plt.ylim(Abs_y_lim_low,Abs_y_lim_high)
                    plt.xticks(fontsize= axis_number_size)
                    plt.yticks(fontsize= axis_number_size)
                    
                    for k,wl in enumerate(wavelengths):
                        if wl == 227:
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Kircher et al. experiment')
                            
                        if wl == 215:
                            plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Kappel et al. experiment')
    
                        plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['absorbance_calculated_from_model'][wl],'b',label='MSI')
                        plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['absorbance_calculated_from_model'][wl],'r',label= "$\it{a}$ $\it{priori}$ model")
                        
                        plt.plot(df_kircher['time'],df_kircher['abs'],'g:',label=r'Kircher et al. (k$_1$)')
                        
                        plt.xlabel('Time [ms]',fontsize=axis_label_size)
                       # plt.ylabel('Absorbance'+' '+str(wl)+' nm',  fontsize=16)
                        plt.ylabel('Absorbance',fontsize=axis_label_size)
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
                            
                            

                        plt.legend(ncol=1,frameon=False,prop={'size': 10.5},loc="best")

                        plt.tick_params(direction='in')

                        plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str('abs')+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_.pdf', bbox_inches='tight')
                        plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str('abs')+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_.svg', bbox_inches='tight')

                        
                        #plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                        #plt.savefig('Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)                        
            
                else:        
    
                    #observable_counter=0
                    for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                        if observable == None:
                            continue
                        
                        fig = plt.figure(figsize=(6, 7))
                        H2O2_x_lim_high=1
                        H2O2_x_lim_low=0
                        H2O2_y_lim_low=0
                        H2O2_y_lim_high=500
                        if observable in exp['concentration_observables']:
                            if observable == 'H2O2':
                                plt.subplot(3,1,1)
                                plt.xlim(H2O2_x_lim_low,H2O2_x_lim_high)
                                plt.ylim(H2O2_y_lim_low,H2O2_y_lim_high)
                                plt.xticks(fontsize= axis_number_size)
                                plt.yticks(fontsize= axis_number_size)
    
    
                            plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable]*1e6,'b',label='MSI')
                            plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable]*1e6,'r',label= "$\it{a}$ $\it{priori}$ model")
                            plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Alquaity et al. experiment') 
                           # plt.xlabel('Time (ms)')
                            print(observable)
                            if observable == 'H2O2':
                                plt.ylabel(r'H$_2$O$_2$ [ppm]', fontsize=axis_label_size)
                               
                            elif observable == 'OH':
                                plt.ylabel('OH [ppm]',  fontsize=axis_label_size)
                                
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
                            
                            print('here')
                            plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf', bbox_inches='tight')
                            
                            
                            print(self.working_directory+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.pdf')
                            #plt.savefig('Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'K_'+sp+'.svg', bbox_inches='tight',transparent=True)
                            #plt.savefig(self.working_directory+'/'+'Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight',dpi=1000)
        
        
        
                            observable_counter+=1


















                        
                        

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
            
    def sort_top_uncertainty_weighted_sens(self,top_sensitivity=5):
        #S_matrix_copy = copy.deepcopy(self.S_matrix)
        S_matrix_copy = copy.deepcopy(self.S_matrix_original)
        self.shorten_sigma()
        sigma_csv = self.sigma_uncertainty_weighted_sensitivity_csv
        if bool(sigma_csv):
            df = pd.read_csv(sigma_csv)
            Sig = np.array(df['Sigma'])
            Sig = Sig.reshape((Sig.shape[0],1))
        elif self.sigma_ones==True:
            shape = len(self.short_sigma)
            Sig = np.ones((shape,1))

                       
            
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
            for j,observable in enumerate(exp['mole_fraction_observables'] + 
                                          exp['concentration_observables'] +
                                          exp['ignition_delay_observables']):
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
                    elif re.match('[Ff]low[ -][Rr][eactor]',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Temperature'])
                        observables[i].append(observable)
                        observable_counter+=1                        
                elif observable in exp['concentration_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Time']*1e3)        
                        observables[i].append(observable)                                
                        observable_counter+=1
                    elif re.match('[Jj][Ss][Rr]',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Temperature'])        
                        observables[i].append(observable)                                
                        observable_counter+=1
                    elif re.match('[Ff]low[ -][Rr][eactor]',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['Temperature'])
                        observables[i].append(observable)
                        observable_counter+=1    
                elif observable in exp['ignition_delay_observables']:
                    if re.match('[Ss]hock [Tt]ube',exp['simulation_type']):
                        time_profiles[i].append(exp['experimental_data'][observable_counter]['temperature'])        
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
        
        gas = ct.Solution(self.new_cti)
        reaction_equations = gas.reaction_equations()
        
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
            
            #transform observable list
            observable_list_transformed = []
            for obs in observables_list:
                lst = obs.split('_')
                if lst[0] =='A':
                    reaction_indx = int(lst[1])
                    reaction = reaction_equations[reaction_indx]
                    observable_list_transformed.append('A_'+reaction)
                elif lst[0] =='n':
                    reaction_indx = int(lst[1])
                    reaction = reaction_equations[reaction_indx]
                    observable_list_transformed.append('n_'+reaction)                    
                elif lst[0] =='Ea':
                    reaction_indx = int(lst[1])
                    reaction = reaction_equations[reaction_indx]
                    observable_list_transformed.append('Ea_'+reaction) 
                else:
                    observable_list_transformed.append(obs)
                    
                    
                
                
        
           
        
            return observable_list_transformed
    
    
    def plotting_uncertainty_weighted_sens(self,simulation_to_plot=0):
        sensitivities,top_sensitivities = self.sort_top_uncertainty_weighted_sens()
        observables_list = self.get_observables_list()
        
        if bool(self.sigma_uncertainty_weighted_sensitivity_csv):
            
            sigma_list = self.sigma_list
        else:
            sigma_list = list(self.short_sigma)
        #start here
        time_profiles = self.getting_time_profiles_for_experiments(self.exp_dict_list_optimized)
        list_of_experiment_observables = self.observable_list
        
        
        def subplot_function(number_of_observables_in_simulation,time_profiles,
                             sensitivities,top_sensitivity_single_exp,observables_list,
                             list_of_experiment_observables,experiment_number):
            #plt.figure(figsize=(2,6))
            #stub
            colors=['k','r','b','g','m']
            lines=['-.','-','--',(0,(5,10)),':']
            #plt.figure()
            fig = plt.figure(figsize=(6, 7))
            #print('THIS IS TESTING')
            
            
            #Hong low temp
            x_right = [3.8,1,.4]
            y_bottom = [-.1,-.08,-.3]
            y_top = [.5,1.5,.8] 
            

            # # #Hong high temp
            x_left=[0,0,0]
            x_right = [2,1,.08]
            y_bottom = [-.2,-.35,-1.25]
            y_top = [.5,1.25,2]
            
            
            #Troe low temp
            x_left=[0]            
            x_right = [1.5]
            y_bottom = [-.25]
            y_top = [.3]

            #Troe middle temp
            x_left=[0]            
            x_right = [1.5]
            y_bottom = [-3]
            y_top = [2]
            
            #Troe High Temp
            x_left=[0]            
            x_right = [.65]
            y_bottom = [-6]
            y_top = [1.5]            
            
            #Kircher LTLP
            x_left=[1.5]            
            x_right = [4.75]
            y_bottom = [-.25]
            y_top = [1.5]                


            #Kircher HTLP
            x_left=[3.5]            
            x_right = [20]
            y_bottom = [-.25]
            y_top = [1.75]    


            #Kircher LTHP
            x_left=[4]            
            x_right = [20]
            y_bottom = [-1]
            y_top = [1.5]    
            
            #Kircher HTHP
            x_left=[5]            
            x_right = [20]
            y_bottom = [-1]
            y_top = [1.5]  


            
            for plot_number in range(number_of_observables_in_simulation):
                 
                for c,top_columns in enumerate(top_sensitivity_single_exp[plot_number]):
                    #comment this out
                    plt.subplot(3,1,1)
                    
                    #comment this back in
                    #plt.subplot(number_of_observables_in_simulation,1,plot_number+1)
                    
                    plt.plot(time_profiles[plot_number],sensitivities[plot_number][:,c],linestyle=lines[c],color =colors[c],
                             label = observables_list[top_columns] +'_'+str(sigma_list[top_columns])) 
                    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.5)
                    if list_of_experiment_observables[plot_number] =='H2O':
                        plt.ylabel(r'$\frac{\partial ln(H_2O)}{\partial(x_j)} \sigma_j$',fontsize=16)
                    elif list_of_experiment_observables[plot_number] =='OH':
                        plt.ylabel(r'$\frac{\partial ln(OH)}{\partial(x_j)} \sigma_j$',fontsize=16)
                    else:
                        plt.ylabel(r'$\frac{\partial ln(abs)}{\partial(x_j)} \sigma_j$',fontsize=16)                        
                    
                    top,bottom = plt.ylim()
                    left,right = plt.xlim()
                    plt.xlim(left=x_left[plot_number],right=x_right[plot_number])
                    plt.ylim(top=y_top[plot_number],bottom=y_bottom[plot_number])
                    #plt.legend(fontsize=2)
                    
                    
                    plt.legend(ncol=2,frameon=False,prop={'size': 5},loc="best")

                    plt.tick_params(direction='in')
                    plt.xticks(fontsize= 14)
                    plt.yticks(fontsize= 14)
                    #plt.legend(ncol=5, loc='upper left',bbox_to_anchor=(-.5,-.3))
                    #plt.legend(ncol=3, loc='upper left',bbox_to_anchor=(1.2,2),fontsize=2)
            
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.2)
            #plt.subplots_adjust(left=.125, bottom=.1, right=.9, top=.9, wspace=.1, hspace=.25)

            
            if self.simulation_run==None:

                plt.savefig(self.working_directory+'/'+'UWSA_Experiment_'+str(experiment_number+1)+'.pdf', bbox_inches='tight')
                plt.savefig(self.working_directory+'/'+'UWSA_Experiment_'+str(experiment_number+1)+'.svg', bbox_inches='tight')

            else:
                
                plt.title('Experiment_'+str(self.simulation_run))
                plt.savefig(self.working_directory+'/'+'UWSA_Experiment_'+str(self.simulation_run)+'.pdf', bbox_inches='tight')

               
        for x in range(len(sensitivities)):   
            if x == simulation_to_plot:
                number_of_observables_in_simulation = len(sensitivities[x])
                subplot_function(number_of_observables_in_simulation,time_profiles[x],sensitivities[x],top_sensitivities[x],observables_list,list_of_experiment_observables[x],x)
            
        return                        
                        
                        
                        