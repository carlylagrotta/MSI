# -*- coding: utf-8 -*-
"""

Execute the nominal and perturbed PAPR-MESS input files,
extractor the rate constants and fit them into PLOG format.

@author: Lei Lei
"""

import Mess_executor as ME
import os, io, sys
import numpy as np
from scipy.optimize import curve_fit

class PAPR_MESS:
    
    def __init__(self, input_name, nom_dict, pert_dict, channels):
        self.nominal_file = input_name
        self.pert_dict = pert_dict
        self.nom_dict = nom_dict
        self.channels = channels
        
    def Run(self):
        """Generate the perturbed PAPR-MESS inputs, execute them, obtain the temperature-dependent rate constants
           for the specified channels."""
        model = ME.Mess_Executor(self.nominal_file, self.pert_dict, self.nom_dict)
        # !!! need to write a function to generate the folder for all possible system
        self.twd = model.new_trial_directory()
        model.new_calculation(self.channels, self.twd)
        
        if model.pert_P != model.nom_P:
            sys.exit('Error: perturbed and nominal systems have different pressure range.')
        self.pertb_ls = model.pertb
        self.nwd = model.nwd
        self.pwd = model.pwd
        self.mwd = model.mwd
        self.input_name = model.input_name
        
    def fit_Cheb_rates(self, n_P, n_T):
        """Fit the rate constants into Chebyshev polynomials."""
        # Execute the nominal and perturbed PAPR-MESS files
        self.Run()
        self.n_P = n_P      # number of perssure degree
        self.n_T = n_T      # number of temperature degree
        self.Cheb_coef = {}
        self.pert_ls = {}
        # read in the temperature-dependent rate constants for each channel of perturbed files
        for wd in [self.nwd, self.pwd]:
            os.chdir(wd)

            if wd == self.nwd:
                wpert = self.nom_dict
            else:
                wpert = self.pert_dict

            fhand = io.open('T_rate.csv', 'rb')
            lines = fhand.readlines()
            fhand.close()
            
            file_flag = False
            pres_flag = False
            channel_flag = False
            final_line = False
            
            rate = {}
            pert = {}
            coef_dict = {}
            # generate a rate constants dictionary for each perturbation
            for key in wpert.keys():
                rate[key] = {}
                pert[key] = {}
                coef_dict[key] = {}

            T_ls = []
            chan_rate = []
            self.P_ls = []
            
            for num, line in enumerate(lines):
                if num == len(lines) - 1:
                    final_line = True
                    T_ls.append(float(line.split(',')[0]))
                    chan_rate.append(float(line.split(',')[1]))
                line = line.strip()
                if line.startswith('=') or final_line:
                    file_flag = True
                    # fit Chebyshev formula and write them in files
                    if len(T_ls) > 2:
                        T_ls = np.array(T_ls)
                        chan_rate = np.array(chan_rate)
                        if self.channel in rate[self.key].keys():
                            rate[self.key][self.channel] = np.concatenate((rate[self.key][self.channel], chan_rate))
                        else:
                            rate[self.key][self.channel] = np.array(chan_rate)
                        if any(chan_rate < 0):
                            sys.exit("Error: Negative rate constants detected...")
                        self.T_ls = T_ls
                        T_ls = []
                        chan_rate = []
                    continue
                if file_flag:
                    file_flag = False
                    pres_flag = True
                    line = line[:-4]
                    self.system = line.split('_')[0]
                    self.key = '_'.join(line.split('_')[1:-1])
                    self.pertb = float(line.split('_')[-1])
                    pert[self.key] = self.pertb
                    continue
                if pres_flag:
                    pres_flag = False
                    channel_flag = True
                    self.pressure = float(line.split(' ')[0])
                    if self.pressure in self.P_ls:
                        continue
                    else:
                        self.P_ls.append(self.pressure)
                    continue
                if channel_flag:
                    channel_flag = False
                    self.channel = line.split(',')[1]
                    continue
                if not file_flag:
                    T_ls.append(float(line.split(',')[0]))
                    chan_rate.append(float(line.split(',')[1]))
                    
            if os.path.exists('Chebyshev_fit.txt'):
                fhand = io.open('Chebyshev_fit.txt', 'ab')
            else:
                fhand = io.open('Chebyshev_fit.txt', 'wb')
            for spc in rate.keys():
            	fhand.write('=' * 30 + '\n')
            	fhand.write(spc + '\n')
            	for key in rate[spc].keys():
                     k = rate[spc][key]
                     k = np.log(k)
                     fhand.write(key)
                     coef = self.fit_cheby_poly(n_P, n_T, k, self.T_ls, self.P_ls)
                     coef_dict[spc][key] = coef
                     fhand.write(str(coef) + '\n')
            fhand.close()
            # store the fitted coefficients for Chebyshev polynomials
            if wd == self.nwd:
                self.Cheb_coef['nominal'] = coef_dict
                self.pert_ls['nominal'] = pert
            else:
                self.Cheb_coef['perturbed'] = coef_dict
                self.pert_ls['perturbed'] = pert
                
            
        print("Fitting channel-specific rate conctants into Chebyshev form for system %s ..." %self.input_name.split(".")[0])
    
    def first_cheby_poly(self, x, n):
        '''Generate n-th order Chebyshev ploynominals of first kind.'''
        if n == 0: return 1
        elif n == 1: return x
        result = 2. * x * self.first_cheby_poly(x, 1) - self.first_cheby_poly(x, 0)
        m = 0
        while n - m > 2:
            result = 2. * x * result - self.first_cheby_poly(x, m+1)
            m += 1
        return result
    
    def reduced_T(self, T, T_min, T_max):
        '''Calculate the reduced temperature.'''
        T_tilde = 2. * T ** (-1) - T_min ** (-1) - T_max ** (-1)
        T_tilde /= (T_max ** (-1) - T_min ** (-1))
        return T_tilde
        
    def reduced_P(self, P, P_min, P_max):
        '''Calculate the reduced pressure.'''
        P_tilde = 2. * np.log(P) - np.log(P_min) - np.log(P_max)
        P_tilde /= (np.log(P_max) - np.log(P_min))
        return P_tilde
        
    def fit_cheby_poly(self, n_T, n_P, k, T_ls, P_ls):
        '''Fit the Chebyshev polynominals to rate constants.
           Input rate constants vector k should be arranged based on pressure.'''
        cheb_mat = np.zeros((len(k), n_T * n_P))
        for n, P in enumerate(P_ls):       # !! assume that at each presssure, we have the same temperateure range
            P_min = P_ls[0]
            P_max = P_ls[-1]
            for m, T in enumerate(T_ls):
                T_min = T_ls[0]
                T_max = T_ls[-1]
                for i in range(n_P):
                    P_tilde = self.reduced_P(P, P_min, P_max)
                    P_cheb = self.first_cheby_poly(P_tilde, i+1)
                    for j in range(n_T):
                        T_tilde = self.reduced_T(T, T_min, T_max)
                        T_cheb = self.first_cheby_poly(T_tilde, j+1)
                        cheb_mat[n*len(T_ls)+m, i*n_T+j] = P_cheb * T_cheb
                        
        coef = np.linalg.pinv(cheb_mat)
        coef = np.dot(coef, k)
        return coef
        
    def fit_Arr_perturbed_rates(self):
        """Fit rate constants into Arrhenius formula."""
        # Execute the nominal and perturbed PAPR-MESS files
        self.Run()
        
        # read in the temperature-dependent rate constants for each channel of perturbed files
        for wd in [self.nwd, self.pwd]:
            os.chdir(wd)
            
            fhand = io.open('T_rate.csv', 'rb')
            lines = fhand.readlines()
            fhand.close()
            
            file_flag = False
            pres_flag = False
            channel_flag = False
            final_line = False
            
            temp = []
            rate = []
            
            for line in lines:
                if line == lines[-1]:
                    final_line = True
                    temp.append(float(line.split(',')[0]))
                    rate.append(float(line.split(',')[1]))
                line = line.strip()
                if line.startswith('=') or final_line:
                    file_flag = True
                    # fit the three-parameter Arrhenius formula and write them in files
                    if len(temp) > 2:
                        temp = np.array(temp)
                        rate = np.array(rate)
                        if any(rate < 0):
                            sys.exit("Error: Negative rate constants detected...")
                        fit = self.log_three_para_Arr_fit(temp, rate)
                        temp = []
                        rate = []
                        if os.path.exists('Arrhenius_fit.txt'):
                            fhand = io.open('Arrhenius_fit.txt', 'ab')
                        else:
                            fhand = io.open('Arrhenius_fit.txt', 'wb')
                        fhand.write('=' * 30 + '\n')
                        fhand.write(self.system + '\n')
                        fhand.write(self.pertb + '\n')
                        fhand.write(self.pressure + '\n')
                        fhand.write(self.channel + str(fit) + '\n')
                        fhand.close()
                    continue
                if file_flag:
                    file_flag = False
                    pres_flag = True
                    line = line[:-4]
                    self.system = line.split('_')[0]
                    self.pertb = '_'.join(line.split('_')[1:])
                    continue
                if pres_flag:
                    pres_flag = False
                    channel_flag = True
                    self.pressure = line
                    continue
                if channel_flag:
                    channel_flag = False
                    self.channel = line.split(',')[1]
                    continue
                if not file_flag:
                    temp.append(float(line.split(',')[0]))
                    rate.append(float(line.split(',')[1]))
                    
        print("Fitting channel-specific rate conctants into Arrhenius form for system %s ..." %self.input_name.split(".")[0])
                
            
    def log_three_para_Arr_fit(self, Temp, rate, ini_guess = (1,1,1), max_initeration = 1000000):
        '''Fit three-parameter Arrhenius rate coefficient expression.'''
        rate = np.log(rate)
        func = lambda T, A, n, Ea: np.log(A) + n * np.log(T) - Ea / T
        fit = curve_fit(func, Temp, rate, p0 = ini_guess, maxfev = max_initeration, ftol = 1E-11, xtol = 1E-11, gtol = 1E-11)
        return fit[0]

    def Cheb_sens_coeff(self, n_P, n_T):
        """Calculate the sensitivity coefficients for the Chebyshev rate constants."""
        self.fit_Cheb_rates(n_P, n_T)   # fit the rate constants
        self.Cheb_sens = {}
        
        for key in self.pert_dict.keys():
            # write the sensitivity coefficients into file
            os.chdir(self.twd)
            if os.path.exists('Chebyshev_sens.txt'):
                fhand = io.open('Chebyshev_sens.txt', 'ab')
            else:
                fhand = io.open('Chebyshev_sens.txt', 'wb')
            fhand.write('='*30 + '\n')
            fhand.write(key + '\n')
            
            Cheb_sens = {}
            nom_key = key[:-1] + '1'
            if len(self.Cheb_coef['nominal'].keys()) == 1:
                nom_key = self.Cheb_coef['nominal'].keys()[0]
            # calculate the channel-specific sensitivity coefficients
            for chan in self.Cheb_coef['perturbed'][key].keys():
                rate_diff = self.Cheb_coef['perturbed'][key][chan] - self.Cheb_coef['nominal'][nom_key][chan]
                sens = rate_diff / (self.pert_ls['perturbed'][key] - self.pert_ls['nominal'][nom_key])
                Cheb_sens[chan] = sens
                fhand.write(str(chan))
                fhand.write(str(sens) + '\n')
            self.Cheb_sens[key] = Cheb_sens
        fhand.close()
        os.chdir(self.mwd)
            
        print("Calculating channel-specific sensitivity coefficients for Chebyshev polynomials for system %s ..." %self.input_name.split('.')[0])
        
        
