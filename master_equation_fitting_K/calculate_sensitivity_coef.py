import numpy as np
import matplotlib.pyplot as plt
import copy 
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
    
    coef = fit_cheby_poly_1d(number_coefficients,k_six_parameter_fit,Temperature)
    k_chevy = calc_polynomial(Temperature,coef)
    
    plt.semilogy(Temperature,k_six_parameter_fit)
    plt.semilogy(Temperature,k_chevy)
    #check to make sure that these things match, and that the
    #coefficients we chose does not overfit
    return coef
    
def perturb_parameter_find_new_alphas(dk,
                                      reaction,
                                      dictonary,
                                      T_min,
                                      T_max,
                                      number_coefficients):
    change_alpha_over_change_parameter = {}
    change_alpha_over_change_parameter_by_reaction = {}
    original_dictonary = copy.deepcopy(dictonary)
    for parameter in dictonary[reaction].keys():
       dictonary[reaction][parameter] = 1.01*dictonary[reaction][parameter]
       
       #reset dictonary
       new_set_of_coefficients = fit_cheby_k_from_six_parameter_fit(T_min,
                                                                    T_max,
                                                                    reaction,
                                                                    dictonary,
                                                                    number_coefficients)
       
       dictonary[reaction][parameter] = original_dictonary[reaction][parameter]
       original_set_of_coefficients = fit_cheby_k_from_six_parameter_fit(T_min,
                                                                         T_max, 
                                                                         reaction,
                                                                         original_dictonary,
                                                                          number_coefficients)
       change_in_coefficients = np.array(new_set_of_coefficients) - np.array(original_set_of_coefficients)
       if parameter == 'A':
           change_in_parameter = np.log(original_dictonary[reaction][parameter]*1.01)-np.log(original_dictonary[reaction][parameter])
           #print(change_in_parameter, parameter)
           sensitivity_coef = np.true_divide(change_in_coefficients,change_in_parameter)

       else:
           change_in_parameter = original_dictonary[reaction][parameter]*1.01 - original_dictonary[reaction][parameter]
          # print(change_in_parameter,parameter)
           sensitivity_coef = np.true_divide(change_in_coefficients,change_in_parameter)
    
       change_alpha_over_change_parameter[parameter] = sensitivity_coef
       change_alpha_over_change_parameter_by_reaction[reaction] = change_alpha_over_change_parameter

    return change_alpha_over_change_parameter_by_reaction


def multiply_by_six_parameter_fit_coeff(T_min,
                                        T_max,
                                        reaction,
                                        dk,
                                        number_coefficients,
                                        dictonary,
                                        six_parameter_fit_sensitivity_dict):
    
    mapped_sensitivity_dict = {}
    change_alpha_over_change_parameter_by_reaction = perturb_parameter_find_new_alphas(dk,    
                                                                           reaction,
                                                                           dictonary,
                                                                           T_min,
                                                                           T_max,
                                                                           number_coefficients)
    
    for i in range(number_coefficients):
        alpha_key = 'alpha_'+str(i)
        alpha_list = []
        temp_dict = {}
        for parameter in six_parameter_fit_sensitivity_dict[reaction].keys():
            parameter_key = parameter
            paramter_list = []
            for k, molecular_parameter in enumerate(six_parameter_fit_sensitivity_dict[reaction][parameter]): 
                change_alpha_over_change_molecular_parameter = molecular_parameter * change_alpha_over_change_parameter_by_reaction[reaction][parameter][i]
                print('alpha',i,parameter,'    ','molecular parameter:',k)
                paramter_list.append(change_alpha_over_change_molecular_parameter)
            temp_dict[parameter]=paramter_list
            mapped_sensitivity_dict[alpha_key] = temp_dict
    
    return mapped_sensitivity_dict

    
    
    

