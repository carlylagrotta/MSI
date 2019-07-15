import numpy as np 
import matplotlib.pyplot as plt 
import cantera as ct 

#temperature and k included here for ease of running code 
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

def one_d_cheby_fit(T,k,deg=10):
    #helper function one which calculates the coef for the fit 
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
    #calculate rate constants helper function and does summation over polynomials 
    T_reduced_list = calculate_reduced_T(T)
    
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)

    return values

def calc_poly_new_gridpoints(T,alpha):
    T_reduced_list = calculate_reduced_T(T)
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)
    
    
    return values
def run_cantera_calculate_rate_constant(T,cti_file):
    # calculated rate constants using cantera 
    gas = ct.Solution(cti_file)
    cantera_k = []
    cantera_T = []
    #for Temp in T:
    for Temp in np.linspace(200,2400,50):
        gas.TPX = Temp,101325,{'Ar':1}
        cantera_k.append(gas.forward_rate_constants[0]*1000)
        cantera_T.append(Temp)
        #multiply by 1000 to convert from cantera units to units rate constats are in 
    return cantera_k,cantera_T
        
    
def run_fitter(T,k,deg=15,cti_file='',og_k = []):
    #call this function to calculate coefficients and plot the fit for the 
    #rate constant as well as plotting it compared to the original rate constant
    # call this function to 
    
    alpha = one_d_cheby_fit(T,k,deg=deg)
    alpha2 = np.array(alpha)
    alpha2 = alpha2.reshape((alpha2.shape[0],1))
    k_calculated  = calc_polynomial(T,alpha)
    #k_new_calculated = calc_poly_new_gridpoints(np.linspace(200,2400,50),alpha)
    cantera_k,cantera_T = run_cantera_calculate_rate_constant(T,cti_file)
    
    plt.figure()
    plt.semilogy(T,k,np.linspace(200,2400,50),k_new_calculated)
    
    plt.figure()
    plt.semilogy(T,k,T,k_calculated)
    plt.title('k calculated with python')
    plt.figure()

    plt.semilogy(np.array(cantera_T),np.array(cantera_k),T,og_k)
    #plt.title('k calculated with cantera')
    return alpha2,k,cantera_k,k_new_calculated
    
# change the path of the cti file to run the code    
alpha,k,cantera_k,k_new_calculated = run_fitter(T,np.log10(k),deg=3,cti_file ='MSI/master_equation_fitting_K/cheb_test.cti',og_k=k )
    