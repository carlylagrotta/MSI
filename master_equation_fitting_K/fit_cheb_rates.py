import os 
import io
import numpy as np
import sys
import numpy as np
import matplotlib.pyplot as plt






#code written using built in numpy functions 


def one_d_cheby_fit(T,k,deg=10):
    #helper function one which calculates the coef.
    T = calculate_reduced_T(T)
    alpha= np.polynomial.chebyshev.Chebyshev.fit(T,k,deg)
    return list(alpha)

def calculate_reduced_T(T_list):
    #helper function that calculates reduced temperature 
    T_max = max(T_list)
    T_min = min(T_list)
    T_reduced = []
    for T in T_list:
        T_tilde = 2. * T ** (-1) - T_min ** (-1) - T_max ** (-1)
        T_tilde /= (T_max ** (-1) - T_min ** (-1))
        T_reduced.append(T_tilde)
    return T_reduced
        
    
def calc_polynomial(T,alpha):
    #calculate rate constants helper function
    T_reduced_list = calculate_reduced_T(T)
    
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)

    return values

def run_fitter(T,k,deg=10):
    #call this function to calculate coefficients and plot the fit for the 
    #rate constant as well as plotting it compared to the original rate constant
    alpha = one_d_cheby_fit(T,k,deg=deg)
    y  = calc_polynomial(T,alpha)
    plt.figure()
    plt.semilogy(T,k,T,y)
    

































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
    T_tilde = 2. * T ** (-1) - T_min ** (-1) - T_max ** (-1)
    T_tilde /= (T_max ** (-1) - T_min ** (-1))
    return T_tilde
    
def reduced_P(P, P_min, P_max):
    '''Calculate the reduced pressure.'''
    P_tilde = 2. * np.log10(P) - np.log10(P_min) - np.log10(P_max)
    P_tilde /= (np.log10(P_max) - np.log10(P_min))
    return P_tilde
    
def fit_cheby_poly(n_T, n_P, k, T_ls, P_ls):
    '''Fit the Chebyshev polynominals to rate constants.
       Input rate constants vector k should be arranged based on pressure.'''
    cheb_mat = np.zeros((len(k), n_T * n_P))
    for n, P in enumerate(P_ls):       # !! assume that at each presssure, we have the same temperateure range
        P_min = P_ls[0]
        P_max = P_ls[-1]
        for m, T in enumerate(T_ls):
            T_min = T_ls[0]
            T_max = T_ls[-1]
            for i in range(n_T):
                T_tilde = reduced_T(T, T_min, T_max)
                T_cheb = first_cheby_poly(T_tilde, i+1)
                for j in range(n_P):
                    P_tilde = reduced_P(P, P_min, P_max)
                    #P_tilde = 0.
                    P_cheb = first_cheby_poly(P_tilde, j+1)
                    cheb_mat[n*len(T_ls)+m, i*n_P+j] = P_cheb * T_cheb
                    
    coef = np.linalg.pinv(cheb_mat)
    coef = np.dot(coef, np.log10(np.array(k)))
    return coef


def fit_cheby_poly_1d(n_T, k, T_ls):
    '''Fit the Chebyshev polynominals to rate constants.
       Input rate constants vector k should be arranged based on pressure.'''
    cheb_mat = np.zeros((len(k), n_T))
    for m, T in enumerate(T_ls):
        T_min = T_ls[0]
        T_max = T_ls[-1]
        for i in range(n_T):
            T_tilde = reduced_T(T, T_min, T_max)
            T_cheb = first_cheby_poly(T_tilde, i)
            #print(T_cheb)
            cheb_mat[m,i] =  T_cheb
            log_k = np.log10(np.array(k))
            #cheb_mat[ i*1+1] = T_cheb
            #cheb_mat[1*len(T_ls)+m, i*1] = T_cheb
            #cheb_mat[i] = T_cheb        
            
            
    #coef = np.linalg.lstsq(cheb_mat,log_k)
    coef,b,c,d = np.linalg.lstsq(cheb_mat,np.log10(k),rcond=-1)
    #coef,b,c,d = np.linalg.lstsq(cheb_mat,log_k,rcond=-1)
    
    #coef = np.dot(coef, np.log10(np.array(k)))
    return coef,log_k,cheb_mat,b



def fit_cheby_poly_zero(n_T, n_P, k, T_ls, P_ls):
    '''Fit the Chebyshev polynominals to rate constants.
       Input rate constants vector k should be arranged based on pressure.'''
    cheb_mat = np.zeros((len(k), n_T * n_P))
    for n, P in enumerate(P_ls):       # !! assume that at each presssure, we have the same temperateure range
        P_min = P_ls[0]
        P_max = P_ls[-1]
        for m, T in enumerate(T_ls):
            T_min = T_ls[0]
            T_max = T_ls[-1]
            for i in range(n_T):
                T_tilde = reduced_T(T, T_min, T_max)
                T_cheb = first_cheby_poly(T_tilde, i)
                for j in range(n_P):
                    #P_tilde = reduced_P(P, P_min, P_max)
                    P_tilde = 0.
                    P_cheb = first_cheby_poly(P_tilde, j)
                    cheb_mat[n*len(T_ls)+m, i*n_P+j] = P_cheb * T_cheb
    
   # k = np.log10(np.array(k))
    coef,b,c,d = np.linalg.lstsq(cheb_mat,k,rcond=-1)
    
    return coef,b



def one_d_cheby_fit(T,k,deg=15):
    #log_k = np.log10(k)
    T = calculate_reduced_T(T)
    alpha= np.polynomial.chebyshev.Chebyshev.fit(T,k,deg)
    return list(alpha)

def calculate_reduced_T(T_list):
    T_max = max(T_list)

    T_min = min(T_list)

    T_reduced = []
    for T in T_list:
        T_tilde = 2. * T ** (-1) - T_min ** (-1) - T_max ** (-1)
        T_tilde /= (T_max ** (-1) - T_min ** (-1))
        T_reduced.append(T_tilde)
    return T_reduced
        
    
def calc_polynomial(T,alpha):
    T_reduced_list = calculate_reduced_T(T)
    
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)

    return values

def run_fitter(T,k,deg=10):
    alpha = one_d_cheby_fit(T,k,deg=deg)
    y  = calc_polynomial(T,alpha)
    plt.figure()
    plt.semilogy(T,k,T,y)
    
#def run_lei_fitter()
    
