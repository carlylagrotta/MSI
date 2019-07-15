
from textwrap import wrap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors 
class overlapping_plots(object):
    def __init__(self,
                 exp_dict_list_optimized,
                 exp_dict_list_original,
                 parsed_yaml_list,
                 list_of_exp_dict_list_optimized=[],
                 list_of_exp_dict_list_original=[],
                 working_directory='/home/carly/Dropbox/Columbia'):



        self.exp_dict_list_optimized = exp_dict_list_optimized
        self.exp_dict_list_original = exp_dict_list_original
        self.parsed_yaml_list = parsed_yaml_list
        self.list_of_exp_dict_list_optimized = list_of_exp_dict_list_optimized
        self.list_of_exp_dict_list_original = list_of_exp_dict_list_original
        self.working_directory = working_directory
    
    def plotting_observables(self,sigmas_original_list=[],sigmas_optimized_list=[]):
        number = 10
        cmap = plt.get_cmap('tab10')
        cmap2 = plt.get_cmap('autumn')
        blue_colors = [cmap(i) for i in np.linspace(0, 1, number)]
        red_colors = [cmap2(i) for i in np.linspace(0, 1, number)]
        marker_list = ['+','o','D','*']
        line_type = ['-', '--', '-.', ':']
        line_width = [7,2,1,1]
     
        
        overall_observables = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            temp = []
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                temp.append(observable)
            overall_observables.append(temp)
                
        
        
            
#            # shouldn't be looping over all of this 
        for i,simulation in enumerate(overall_observables):
            if 'perturbed_coef' in self.list_of_exp_dict_list_optimized[0][i].keys():
                absorb = True
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
            else:
                absorb = False
            observable_counter=0

            for j,observable in enumerate(overall_observables[i]):
                if observable==None:
                    continue
                #print(observable,i)
                
                if observable in self.list_of_exp_dict_list_optimized[0][i]['concentration_observables']:
                    plt.figure()
                    for k,exp in enumerate(self.list_of_exp_dict_list_optimized):
                        #print(observable,k)
                       # plt.plot(self.list_of_exp_dict_list_original[k][i]['simulation'].timeHistories[0]['time']*1e3,self.list_of_exp_dict_list_original[k][i]['simulation'].timeHistories[0][observable]*1e6,linewidth=5,color = 'r',label= "$\it{A priori}$ model_"+str(k))

                        plt.plot(exp[i]['simulation'].timeHistories[0]['time']*1e3,exp[i]['simulation'].timeHistories[0][observable]*1e6,color = blue_colors[k],linestyle=line_type[k],marker=marker_list[k],markersize=1.5,linewidth=line_width[k],label='MSI_'+str(k))
                        #k+1+k*-5
                        if k ==0:
                           # plt.plot(exp[i]['experimental_data'][observable_counter]['Time']*1e3,exp[i]['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Experimental Data') 
                            plt.ylabel(observable+' ' + 'ppm')
                            plt.xlabel('Time (ms)')
                            plt.plot([],'w' ,label= 'T:'+ str(self.exp_dict_list_original[i]['simulation'].temperature))
                            plt.plot([],'w', label= 'P:'+ str(self.exp_dict_list_original[i]['simulation'].pressure))
                            for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():                        
                                plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                        plt.legend(ncol=5, loc='upper left',bbox_to_anchor=(-.2,-.2))

                        if bool(sigmas_optimized_list) == True:
                        
                            high_error_optimized = np.exp(sigmas_optimized_list[k][i][observable_counter])                   
                            high_error_optimized = np.multiply(high_error_optimized,exp[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            low_error_optimized = np.exp(sigmas_optimized_list[k][i][observable_counter]*-1)
                            low_error_optimized = np.multiply(low_error_optimized,exp[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            plt.plot(exp[i]['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,color = blue_colors[k],linestyle='--')
                            plt.plot(exp[i]['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,color = blue_colors[k],linestyle='--')
                            
                            
                            
                            #high_error_original = np.exp(sigmas_original_list[k][i][observable_counter])
                            #high_error_original = np.multiply(high_error_original,self.list_of_exp_dict_list_original[k][i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            #low_error_original = np.exp(sigmas_original_list[k][i][observable_counter]*-1)
                            #low_error_original = np.multiply(low_error_original,self.list_of_exp_dict_list_original[k][i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                            #plt.plot(exp[i]['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,color = red_colors[k+1+k*-5],linestyle='--')
                           # plt.plot(exp[i]['experimental_data'][observable_counter]['Time']*1e3,low_error_original,color = red_colors[k+1+k*-5],linestyle='--')  
                plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+'_'+str(observable)+'_'+str(self.exp_dict_list_original[i]['simulation'].temperature)+'.pdf', bbox_inches='tight')
                observable_counter+=1
                
                
            if absorb == True:
                plt.figure()
                for k,exp in enumerate(self.list_of_exp_dict_list_optimized):
                    for p,wl in enumerate(wavelengths):
                        #plt.plot(self.list_of_exp_dict_list_original[k][i]['simulation'].timeHistories[0]['time']*1e3,self.list_of_exp_dict_list_original[k][i]['absorbance_calculated_from_model'][wl],color = 'r',linewidth=5,label= "$\it{A priori}$ model_"+str(k))
                        plt.plot(exp[i]['simulation'].timeHistories[0]['time']*1e3,exp[i]['absorbance_calculated_from_model'][wl],color = blue_colors[k],linestyle=line_type[k],marker=marker_list[k],linewidth=line_width[k],markersize=1.5,label='MSI_'+str(k))                               
                        if k ==0:
                            #plt.plot(exp[i]['absorbance_experimental_data'][p]['time']*1e3,exp[i]['absorbance_experimental_data'][p]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')
                            plt.ylabel('Absorbance'+''+str(wl))
                            plt.xlabel('Time (ms)')
                            plt.plot([],'w' ,label= 'T:'+ str(self.exp_dict_list_original[i]['simulation'].temperature))
                            plt.plot([],'w', label= 'P:'+ str(self.exp_dict_list_original[i]['simulation'].pressure))
                            for key in self.exp_dict_list_original[i]['simulation'].conditions.keys():                        
                                plt.plot([],'w',label= key+': '+str(self.exp_dict_list_original[i]['simulation'].conditions[key]))
                        plt.legend(ncol=5, loc='upper left',bbox_to_anchor=(-.2,-.2))
                        
                        
                        if bool(sigmas_optimized_list)==True:
                            high_error_optimized = np.exp(sigmas_optimized_list[k][i][observable_counter])
                            high_error_optimized = np.multiply(high_error_optimized,exp[i]['absorbance_model_data'][wl])
                            low_error_optimized = np.exp(sigmas_optimized_list[k][i][observable_counter]*-1)
                            low_error_optimized = np.multiply(low_error_optimized,exp[i]['absorbance_model_data'][wl])
                            
                            plt.plot(exp[i]['absorbance_experimental_data'][p]['time']*1e3,high_error_optimized,color = blue_colors[k],linestyle='--')
                            plt.plot(exp[i]['absorbance_experimental_data'][p]['time']*1e3,low_error_optimized,color = blue_colors[k],linestyle='--')
                        
                        plt.savefig(self.working_directory+'/'+'Exp_'+str(i+1)+' '+'Absorb at'+'_'+str(wl)+str(self.exp_dict_list_original[i]['simulation'].temperature)+'.pdf', bbox_inches='tight')

                            #high_error_original = np.exp(sigmas_original_list[k][i][observable_counter])
                            #high_error_original = np.multiply(high_error_original,self.list_of_exp_dict_list_original[k][i]['absorbance_model_data'][wl])
                            #low_error_original =  np.exp(sigmas_original_list[k][i][observable_counter]*-1)
                            #low_error_original = np.multiply(low_error_original,self.list_of_exp_dict_list_original[k][i]['absorbance_model_data'][wl])
                            
                            
                            #plt.plot(exp[i]['absorbance_experimental_data'][p]['time']*1e3,high_error_original,color = red_colors[k+1+k*-5],linestyle='--')
                            #plt.plot(exp[i]['absorbance_experimental_data'][p]['time']*1e3,low_error_original,color = red_colors[k+1+k*-5],linestyle='--')

                

        
        
        
