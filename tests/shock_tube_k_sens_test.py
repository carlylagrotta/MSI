#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 12:53:15 2018

@author: carly
"""

import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct

test_p = pr.Processor('MSI/data/test_data/FFCM1_custom.cti')
#test_p = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb.cti')
#test_p = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb_test.cti')

test_tube = st.shockTube(pressure=1.909,
                         temperature=1398,
                         observables=['OH','H2O'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'O2':0.00062 ,
                                     'H2O':0.001234,
                                     'H2O2':0.00254 ,
                                     'Ar': 0.9974968858350951},
                         initialTime=0,
                         finalTime=.0002,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

test_tube.run()

ksens = test_tube.kineticSensitivities
print(ksens.shape)


test_p2 = pr.Processor('MSI/data/chebyshev_data/FFCM1_custom_cheb_test.cti')

test_tube = st.shockTube(pressure=1.909,
                         temperature=1398,
                         observables=['OH','H2O'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'O2':0.00062 ,
                                     'H2O':0.001234,
                                     'H2O2':0.00254 ,
                                     'Ar': 0.9974968858350951},
                         initialTime=0,
                         finalTime=.0002,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant volume',
                         processor=test_p2,
                         save_timeHistories=1,
                         save_physSensHistories=1)


test_tube.run()

ksens2 = test_tube.kineticSensitivities