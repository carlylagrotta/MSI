import numpy as np
import matplotlib.pyplot as plt
import copy 
import scipy.optimize as spcf
import yaml
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

def calculate_list_of_k(T_initial,T_final,reaction,dictonary):
    temperature = np.arange(T_initial,T_final)
    k_list = []
    for temp in temperature:
        k_list.append(calculate_six_parameter_fit(reaction,dictonary,temp))
    return temperature,k_list
        
    
def first_cheby_poly(x, n):
    '''Generate n-th order Chebyshev ploynominals of first kind.'''
    if n == 0: return 1
    elif n == 1: return x
    result = 2. * x * first_cheby_poly(x, 1) - first_cheby_poly(x, 0)
    m = 0
    while n - m > 2:
        result = 2. * x * result - first_cheby_poly(x, m+1)
        m += 1
    return result

def reduced_T( T, T_min, T_max):
    '''Calculate the reduced temperature.'''
    T = np.array(T)
    T_tilde = 2. * T ** (-1.0) - T_min ** (-1.0) - T_max ** (-1.0)
    T_tilde /= (T_max ** (-1.0) - T_min ** (-1.0))
    return T_tilde
    
def calc_polynomial(T,alpha):
    #calculate rate constants helper function
    T_reduced_list = reduced_T(T,200,2400)
    
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)

    return values



def fit_cheby_poly_1d(n_T, k, T_ls):
    #T needs to be a lsit 
    '''Fit the Chebyshev polynominals to rate constants.
       Input rate constants vector k should be arranged based on pressure.'''
    cheb_mat = np.zeros((len(k), n_T))
    for m, T in enumerate(T_ls):
        T_min = T_ls[0]
        T_max = T_ls[-1]
        for i in range(n_T):
            T_tilde = reduced_T(T, T_min, T_max)
            T_cheb = first_cheby_poly(T_tilde, i)

            cheb_mat[m,i] =  T_cheb
            #log_k = np.log10(np.array(k))

            
    coef,b,c,d = np.linalg.lstsq(cheb_mat,k,rcond=-1)
    
    return coef

def fit_cheby_k_from_six_parameter_fit(T_min,
                                       T_max,
                                       reaction,
                                       dictonary,
                                       number_coefficients):
    k_six_parameter_fit = []
    Temperature = []
    for temp in np.arange(T_min,T_max):
        Temperature.append(temp)
        k_six_parameter_fit.append(calculate_six_parameter_fit(reaction,dictonary,temp))
    
    coef = fit_cheby_poly_1d(number_coefficients,np.log10(np.array(k_six_parameter_fit)),Temperature)
    k_chevy = calc_polynomial(Temperature,coef)
    
    #plt.semilogy(Temperature,k_six_parameter_fit)
   # plt.semilogy(Temperature,k_chevy)
    #check to make sure that these things match, and that the
    #coefficients we chose does not overfit
    return coef,k_chevy
 
def fit_cheby_k_from_six_parameter_fit_ln(T_min,
                                       T_max,
                                       reaction,
                                       dictonary,
                                       number_coefficients):
    k_six_parameter_fit = []
    Temperature = []
    for temp in np.arange(T_min,T_max):
        Temperature.append(temp)
        k_six_parameter_fit.append(calculate_six_parameter_fit(reaction,dictonary,temp))
    
    coef = fit_cheby_poly_1d(number_coefficients,np.log(np.array(k_six_parameter_fit)),Temperature)
    k_chevy = calc_polynomial(Temperature,coef)
    
    plt.semilogy(Temperature,k_six_parameter_fit)
    plt.semilogy(Temperature,np.exp(k_chevy))
    #check to make sure that these things match, and that the
    #coefficients we chose does not overfit
    return coef,k_chevy
 
####ONLY WROTE THIS FUNCTION TO CHECK MY VALUES AGAINST LEI###### 
#USING THE ONE ABOVE IT TO DO THE FITTING#####
def fit_cheby_k_from_six_parameter_fit_ln10(T_min,
                                       T_max,
                                       reaction,
                                       dictonary,
                                       number_coefficients):
    k_six_parameter_fit = []
    Temperature = []
    for temp in np.arange(T_min,T_max):
        Temperature.append(temp)
        k_six_parameter_fit.append(calculate_six_parameter_fit(reaction,dictonary,temp))
    
    coef = fit_cheby_poly_1d(number_coefficients,np.log10(np.array(k_six_parameter_fit)),Temperature)
    k_chevy = calc_polynomial(Temperature,coef)
    
    plt.semilogy(Temperature,k_six_parameter_fit)
    plt.semilogy(Temperature,np.exp(k_chevy))
    #check to make sure that these things match, and that the
    #coefficients we chose does not overfit
    return coef,k_chevy           

        
def find_new_alpha(xdata, ydata, alphas, alpha_index_to_recompute):
    def f(xdata, new_alpha):
        ydata = []
        for x in xdata:
            value = 0.0
            for alpha_index, alpha in enumerate(alphas):
                if alpha_index == alpha_index_to_recompute:
                    alpha = new_alpha
                term = first_cheby_poly(x, alpha_index)
                value += term * alpha
            ydata.append(value)
        return ydata
    
    popt, pcov = spcf.curve_fit(f, xdata, ydata)
    return popt, pcov
    
    
def find_alpha_coefficients_sens(original_alphas,new_alphas,dk):
    change_in_coefficients = np.array(new_alphas) - np.array(original_alphas)
    sensitivity = np.true_divide(change_in_coefficients,dk)
    return sensitivity

#convert sensiivty  from log 10 to ln
def multiply_lei_values(sens_dict):
    temp_dict = {}
    for keys in sens_dict:
        temp_dict[keys] = []
    for keys in sens_dict:
        for arr in sens_dict[keys]:
            temp_arr = arr*np.log(10)
            temp_dict[keys].append(temp_arr)
    return temp_dict
            



def change_lei_values_for_testing(sens_dict,values_to_skip=None):
    temp_dict = {}
    for keys in sens_dict:
        print(keys)
        temp_dict[keys] = []
    for skp,keys in enumerate(sens_dict):
        for i,arr in enumerate(sens_dict[keys]):
            print(values_to_skip[skp], keys)
            if i in values_to_skip[skp]:
                temp_arr = arr*np.log(10)
                temp_dict[keys].append(temp_arr)
            else:
                temp_arr = arr*np.log(.01)
                temp_arr = temp_arr/.01
                temp_arr = temp_arr*np.log(10)
                temp_dict[keys].append(temp_arr)
    return temp_dict
def convert_energy_unit(energy_sens_array):
    test_factor = (1/.012)*(1/(1000*1000))
    converted = energy_sens_array * test_factor
    return converted

def write_out_dict(config2,path =''):
    with open(path,'w') as f:
        yaml.safe_dump(config2, f,default_flow_style=False)
#def loop_over_six_parameter_fit(Temp_list,reaction,dictonary):
#    k= []
#    for temp in Temp_list:
#        k.append(calculate_six_parameter_fit(reaction,dictonary,temperature))
#    return k 
#    
#    
#def perturb_parameter_find_new_alphas(dk,
#                                      reaction,
#                                      dictonary,
#                                      T_min,
#                                      T_max,
#                                      number_coefficients):
#    change_alpha_over_change_parameter = {}
#    change_alpha_over_change_parameter_by_reaction = {}
#    original_dictonary = copy.deepcopy(dictonary)
#    original_set_of_coefficients = fit_cheby_k_from_six_parameter_fit(T_min,
#                                                                         T_max, 
#                                                                         reaction,
#                                                                         original_dictonary,
#                                                                          number_coefficients)
#    for parameter in dictonary[reaction].keys():
#        new_alphas = []
#        dictonary[reaction][parameter] = 1.01*dictonary[reaction][parameter]
#        T_reduced = reduced_T(np.arange(T_min,T_max),T_min,T_max)
#        k = loop_over_six_parameter_fit(np.arange(T_min,T_max),reaction,dictonary)
#        
#        for i,alpha in enumerate(original_set_of_coefficients):
#
#            new_alpha = find_new_alpha(T_reduced,k,original_set_of_coefficients,i)
#            new_alphas.append(new_alpha)
#            # calculate new set of alphas for each A,n,Ea
#            #make sure you are passing in reduced temperature 
#            #append those to alpha
#            
#       
#            dictonary[reaction][parameter] = original_dictonary[reaction][parameter]
#
#
#
#
#
#        print(original_set_of_coefficients)
#       change_in_coefficients = np.array(new_set_of_coefficients) - np.array(original_set_of_coefficients)
#       if parameter == 'A':
#           change_in_parameter = np.log(original_dictonary[reaction][parameter]*1.01)-np.log(original_dictonary[reaction][parameter])
#           #print(change_in_parameter, parameter)
#           sensitivity_coef = np.true_divide(change_in_coefficients,change_in_parameter)
#
#       else:
#           change_in_parameter = original_dictonary[reaction][parameter]*1.01 - original_dictonary[reaction][parameter]
#          # print(change_in_parameter,parameter)
#           sensitivity_coef = np.true_divide(change_in_coefficients,change_in_parameter)
#    
#       change_alpha_over_change_parameter[parameter] = sensitivity_coef
#       change_alpha_over_change_parameter_by_reaction[reaction] = change_alpha_over_change_parameter
#
#    return change_alpha_over_change_parameter_by_reaction
#
#
#def multiply_by_six_parameter_fit_coeff(T_min,
#                                        T_max,
#                                        reaction,
#                                        dk,
#                                        number_coefficients,
#                                        dictonary,
#                                        six_parameter_fit_sensitivity_dict):
#    
#    mapped_sensitivity_dict = {}
#    change_alpha_over_change_parameter_by_reaction = perturb_parameter_find_new_alphas(dk,    
#                                                                           reaction,
#                                                                           dictonary,
#                                                                           T_min,
#                                                                           T_max,
#                                                                           number_coefficients)
#    
#    for i in range(number_coefficients):
#        alpha_key = 'alpha_'+str(i)
#        alpha_list = []
#        temp_dict = {}
#        for parameter in six_parameter_fit_sensitivity_dict[reaction].keys():
#            parameter_key = parameter
#            paramter_list = []
#            for k, molecular_parameter in enumerate(six_parameter_fit_sensitivity_dict[reaction][parameter]): 
#                change_alpha_over_change_molecular_parameter = molecular_parameter * change_alpha_over_change_parameter_by_reaction[reaction][parameter][i]
#                paramter_list.append(change_alpha_over_change_molecular_parameter)
#            temp_dict[parameter]=paramter_list
#            mapped_sensitivity_dict[alpha_key] = temp_dict
#    
#    return mapped_sensitivity_dict

    
    
    

