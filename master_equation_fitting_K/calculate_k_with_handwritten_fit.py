import numpy as np 
import cantera as ct
import matplotlib.pyplot as plt 

T = [200,
 300,
 400,
 500,
 600,
 700,
 800,
 900,
 1000,
 1100,
 1200,
 1300,
 1400,
 1500,
 1600,
 1700,
 1800,
 1900,
 2000,
 2100,
 2200,
 2300,
 2400]

k = [3598820205677.96,
 1384365848680.28,
 918526746446.366,
 846020502047.191,
 913357292578.687,
 1064137777998.82,
 1284494707804.24,
 1572988369463.2,
 1932557956912.69,
 2367768739012.29,
 2883668329436.64,
 3485232275971.83,
 4177063735680.52,
 4963220865761.82,
 5847118271780.92,
 6831476911417.45,
 7918308618829.78,
 9108926726341.08,
 10403976846527.2,
 11803483258920.4,
 13306907175947.7,
 14913213735305,
 16620945014500.6]


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
    
def calc_polynomial(T,alpha):
    #calculate rate constants helper function
    T_reduced_list = reduced_T(T)
    
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)

    return values



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

            cheb_mat[m,i] =  T_cheb
            #log_k = np.log10(np.array(k))

            
    coef,b,c,d = np.linalg.lstsq(cheb_mat,k,rcond=-1)
    
    return coef

def run_cantera_calculate_rate_constant(T,cti_file):
    gas = ct.Solution(cti_file)
    cantera_k = []
    for Temp in T:
        gas.TPX = Temp,101325,{'Ar':1}
        cantera_k.append(gas.forward_rate_constants[0])
    return cantera_k*1000

def run_fitter(n_T,k,T_ls,cti_file):
    alpha = fit_cheby_poly_1d(15,k,T)
    k_calculated = calc_polynomial(T,alpha)
    cantera_k = run_cantera_calculate_rate_constant(T,cti_file)
    plt.figure()
    plt.semilogy(T,k,T,k_calculated)
    plt.title('k calculated with python')
    plt.figure()
    plt.semilogy(T,np.array(cantera_k))
    plt.title('k calculated with cantera')
    
    