# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:33:46 2018

@author: Mark Barbet
"""

import cantera as ct
from .. import simulation as sim
from ...cti_core import cti_processor as ctp
import pandas as pd
import numpy as np
import time
import copy
import re

class free_flame(sim.Simulation):
    
    '''child class of sim.Simulaton.  Inherits all attributes and methods including __init__().  
    Also has internal init due to data requirements'''
    
    
    
    
    
    def __init__(self,pressure:float,temperature:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={},flame_width:float=1.0,
                 save_timeHistories:int=0,T_profile=pd.DataFrame(columns=['z','T']),soret=True,
                 tol_ss=[1.0e-5, 1.0e-13],tol_ts=[1.0e-4, 1.0e-10],loglevel=1,flametype='Flame Speed'):
        
        #sim.Simulation.__init__(self,pressure,temperature,observables,kineticSens,physicalSens,
        #                        conditions,processor,cti_path)
        #set up processor and initialize all variables  common to every simulation  
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)
        
        self.pressure=pressure
        self.temperature=temperature
        self.observables=observables
        self.kineticSens=kineticSens
        self.physicalSens=physicalSens
        self.conditions=conditions
        self.cti_path=cti_path
        self.thermalBoundary = thermalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        self.energycon='off'
        self.flame_width=flame_width
        self.timeHistory = None
        self.experimentalData = None
        self.tol_ss=tol_ss
        self.tol_ts=tol_ts
        self.soret=soret
        self.loglevel=loglevel
        self.flametype=flametype
        
        if save_timeHistories == 1:
            self.timeHistories=[]
            self.timeHistoryInterpToExperiment = None
            self.pressureAndTemperatureToExperiment = None
        else:
            self.timeHistories=None
        if save_physSensHistories == 1:
            self.physSensHistories = []
        self.setTPX()
        self.dk = 0.01
        self.solution=None
        
        
    def printVars(self):
        print()
        
    def settingFlameConditions(self):
        '''
        Determine the mechanical and thermal boundary conditions for a 
        shock tube.
        '''
        
        #assigning the thermal boundary variable
        if re.match('[aA]diabatic',self.thermalBoundary):
            energy = 'on'
        elif re.match('[iI]sothermal',self.thermalBoundary):
            energy = 'off'
        else:
            raise Exception('Please specify a thermal boundary condition, adiabatic or isothermal')
        #assigning the mehcanical boundary variable 
        if re.match('[Cc]onstant [Pp]ressure',self.mechanicalBoundary):
            mechBoundary = 'constant pressure'
        elif re.match('[Cc]onstant [Vv]olume',self.mechanicalBoundary):
            mechBoundary = 'constant volume'
        else:
            raise Exception('Please specifiy a mehcanical boundary condition, constant pressure or constant volume')
        #return the thermal and mechanical boundary of the shock tube 
        return energy,mechBoundary
    
    def sensitivity_adjustment(self,temp_del:float=0.0,
                               pres_del:float=0.0,
                               spec_pair:(str,float)=('',0.0)):
        
        #this is where we would make the dk fix
        if temp_del != 0.0:
            self.dk.append(temp_del)
        if pres_del != 0.0:       
            self.dk.append(pres_del) 
        if spec_pair[1] != 0.0:
            self.dk.append(spec_pair[1])
        
        
        kin_temp = self.kineticSens
        self.kineticSens = 0
        data = sim.Simulation.sensitivity_adjustment(self,temp_del,pres_del,spec_pair)
        self.kineticSens = kin_temp
        return data
    
    def run_single(self):
        
        gas=self.processor.solution       
        
        self.flame=ct.FreeFlame(gas,width=self.flame_width)
        
        self.flame.flame.set_steady_tolerances(default=self.tol_ss)   #Set steady state tolerances
        self.flame.flame.set_transient_tolerances(default=self.tol_ts) #Set transient tolerances
        print('Running simulation at T = '+str(round(gas.T,5))+', P = '+str(round(gas.P,5))+'\\Conditions: '+str(self.conditions))
        if re.match('[aA]diabatic',self.thermalBoundary):
            energycon = True
        self.flame.energy_enabled = energycon
        
        self.flame.transport_model = 'Multi'
        self.flame.set_max_jac_age(10,10)
        self.flame.solve(self.loglevel,refine_grid=False)
        self.flame.soret_enabled = self.soret
        
        
        
        
        self.flame.set_refine_criteria(ratio=2.0,slope=0.05,curve=0.1)
        self.flame.transport_model = 'Multi'
        self.flame.solve(self.loglevel,refine_grid=True)
        self.forward_rates=self.flame.forward_rate_constants
        self.net_rates=self.flame.net_rates_of_progress
        self.reverse_rates=self.flame.reverse_rate_constants
        
        
               
        
        
        if re.match('[Ff]lame [Ss]peed',self.flametype):
            columns=['T_in','P']+list(self.conditions.keys())+['u0']
            self.solution=pd.DataFrame(columns=columns)
            for i in range(len(columns)):
                if i==0:
                    self.solution['T_in']=[self.temperature]
                elif i==1:
                    self.solution['P']=[self.pressure]
                elif i>1 and i<len(columns)-1:
                    self.solution[columns[i]]=[self.conditions[columns[i]]]
                elif i==len(columns)-1:
                    self.solution['u0']=[self.flame.u[0]]
                    
        elif re.match('[Aa]diabatic [Ff]lame',self.flametype):
            self.solution=self.flame.X
            columnNames=gas.species_names
            tempdata=pd.DataFrame(columns=columnNames,data=self.flame.X)
            print(tempdata)
        if re.match('[Aa]diabatic [Ff]lame',self.flametype) and self.kineticSens==1 and bool(self.observables):
            self.dk=0.01
            #Calculate kinetic sensitivities
            sensIndex = [self.flame.grid.tolist(),gas.reaction_equations(),self.observables]
            
            S = np.zeros((len(self.flame.grid),gas.n_reactions,len(self.observables)))
            #print(solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1])
            #a=solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1]
            for m in range(gas.n_reactions):
                gas.set_multiplier(1.0)
                gas.set_multiplier(1+self.dk,m)
                self.flame.solve(loglevel=1,refine_grid=False)
                for i in np.arange(len(self.observables)):
                    for k in np.arange(len(self.flame.grid)):                    
                        S[k,m,i]=np.log10(self.solution[self.flame.flame.component_index(self.observables[i])-4,k])-np.log10(self.flame.X[self.flame.flame.component_index(self.observables[i])-4,k])
                        #print(solution.X[solution.flame.component_index(observables[i])-4,k])
                        #print(f.X[f.flame.component_index(observables[i])-4,k])
                        S[k,m,i]=np.divide(S[k,m,i],np.log10(self.dk))    
        
        
        elif self.kineticSens==1 and bool(self.observables)==False and not re.match('[Ff]lame [Ss]peed',self.flametype):
            raise Exception('Please supply a list of observables in order to run kinetic sensitivity analysis')
        #gas.set_multiplier(1.0)
        elif self.kineticSens==1 and re.match('[Ff]lame [Ss]peed',self.flametype):    
            self.fsens = self.flame.get_flame_speed_reaction_sensitivities()
        

        
        if self.kineticSens==1 and re.match('[Ff]lame [Ss]peed',self.flametype):
            #numpyMatrixsksens = [dfs[dataframe].values for dataframe in range(len(dfs))]
            #self.kineticSensitivities = np.dstack(numpyMatrixsksens)
            #print(np.shape(self.kineticSensitivities))
            #self.solution=data
            return (self.solution,self.fsens)
        elif self.kineticSens==1 and re.match('[Aa]diabatic [Ff]lame',self.flametype):
            dfs = [pd.DataFrame() for x in range(len(self.observables))]
            #print((pd.DataFrame(sens[0,:])).transpose())
            #test=pd.concat([pd.DataFrame(),pd.DataFrame(sens[0,:]).transpose()])
            #print(test)
            #for k in range(len(self.observables)):
                #dfs[k] = dfs[k].append(((pd.DataFrame(sens[k,:])).transpose()),ignore_index=True)
                #dfs[k]=pd.concat([dfs[k],pd.DataFrame(sens[k,:]).transpose()])
                #dfs[k]=pd.DataFrame(sens[k,:]).transpose()
            #print(dfs)  
            numpyMatrixsksens = [dfs[dataframe].values for dataframe in range(len(dfs))]
            self.kineticSensitivities = np.dstack(numpyMatrixsksens)
            #print(np.shape(self.kineticSensitivities))
            #self.solution=data
            return (self.solution,self.kineticSensitivities)
        elif self.kineticSens==0:
            #numpyMatrixsksens = [dfs[dataframe].values for dataframe in range(len(dfs))]
            
            return (self.solution,[])
        else:
            #self.solution=data
            return (self.solution,[])
        
       

class flamespeed_multi_condition(sim.Simulation):
        
    def __init__(self,pressures:float,temperatures:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,
                 processor:ctp.Processor=None, 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={},flame_width:float=1.0,
                 save_timeHistories:int=0,T_profile=pd.DataFrame(columns=['z','T']),soret=True,
                 tol_ss=[1.0e-5, 1.0e-13],tol_ts=[1.0e-4, 1.0e-10],loglevel=1,flametype='Flame Speed',cti_path=""):
        
#    sim.Simulation.__init__(self,pressure,temperature,observables,kineticSens,physicalSens,
#                                conditions,processor,cti_path)
        
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)        
        
        self.save_physSensHistories=save_physSensHistories
        self.temperatures=temperatures
        self.JSR_objects=[]
        self.pressures=pressures
        self.observables=observables
        self.kineticSens=kineticSens
        self.physicalSens=physicalSens
        self.conditions=conditions
        self.thermalBoundary=thermalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        self.energycon='off'
        self.flame_width=flame_width
        self.tol_ss=tol_ss
        self.tol_ts=tol_ts
        self.soret=soret
        self.loglevel=loglevel
        self.flametype=flametype
        self.save_timeHistories=save_timeHistories
               
        if save_timeHistories == 1:
            self.timeHistories=[]
            self.timeHistoryInterpToExperiment = None
            self.pressureAndTemperatureToExperiment = None
        else:
            self.timeHistories=None
        
        #self.setTPX()
        self.dk = [0]
       
            
    def run(self):
        
        
        solution=[]
        ksens=[]
        ksens_1stIter=False
        for i in range(len(self.temperatures)):
            for j in range(len(self.pressures)):
                for k in range(len(self.conditions)):
                    temp_flame=free_flame(pressure=self.pressures[j],
                                          temperature=self.temperatures[i],
                                          observables=self.observables,
                                          kineticSens=self.kineticSens,
                                          physicalSens=self.physicalSens,
                                          conditions=self.conditions[k],
                                          thermalBoundary=self.thermalBoundary,
                                          processor=self.processor,
                                          save_physSensHistories=self.save_physSensHistories,
                                          flame_width=self.flame_width,
                                          save_timeHistories=self.save_timeHistories,
                                          soret=self.soret,
                                          tol_ss=self.tol_ss,
                                          tol_ts=self.tol_ts,
                                          loglevel=self.loglevel,
                                          flametype=self.flametype)
            
                    a,b=temp_flame.run_single()
            
                    temp=[]
                    temp1=[]
                    temp=copy.deepcopy(a)
                    temp1=copy.deepcopy(b)
                    #print(a)
                    solution.append(temp)
                    if not ksens_1stIter and self.kineticSens==1:
                        ksens=temp1
                        ksens_1stIter=True
                    elif self.kineticSens==1 and ksens_1stIter:
                        ksens=np.vstack([ksens,temp1])
                    #print(ksens)
        solution=pd.concat(solution)
        #print(np.shape(ksens))
        if self.timeHistories != None:
            self.timeHistories.append(solution)
        self.kineticSensitivities=ksens
        return (solution,ksens)
        
        
    def sensitivity_adjustment(self,temp_del:float=0.0,
                               pres_del:float=0.0,
                               spec_pair:(str,float)=('',0.0),
                               res_del:float=0.0):
        
        #this is where we would make the dk fix
        if temp_del != 0.0:
            self.dk.append(temp_del)
        if pres_del != 0.0:       
            self.dk.append(pres_del) 
        if spec_pair[1] != 0.0:
            self.dk.append(spec_pair[1])
        
        temptemp=copy.deepcopy(self.temperatures)
        temppres=copy.deepcopy(self.pressure)
        tempcond=copy.deepcopy(self.conditions)
        kin_temp = self.kineticSens
        self.kineticSens = 0
        '''
          Passes the Perturbed observable to the setTPX function. Temperature and pressure 
        are passed and set directly species need to go through an additional step in the 
        setTPX function. 
        '''
        if spec_pair[0] != '':
            self.temperatures=np.array(self.temperatures)+temp_del*np.array(self.temperatures)
            self.pressure=self.pressure+pres_del*self.pressure
            xj=self.conditions[spec_pair[0]]
            delxj=spec_pair[1]*self.conditions[spec_pair[0]]
            #print(xj,delxj)
            self.conditions[spec_pair[0]]=np.divide(np.multiply(xj+delxj,1-xj),1-xj-delxj)
#           self.setTPX(self.temperature+self.temperature*temp_del,
#                   self.pressure+self.pressure*pres_del,
#                   {spec_pair[0]:self.conditions[spec_pair[0]]*spec_pair[1]})
          
           
        else:
           self.temperatures=np.array(self.temperatures)+temp_del*np.array(self.temperatures)
           self.pressure=self.pressure+pres_del*self.pressure
           self.residence_time=self.residence_time+res_del*self.residence_time
           
#           self.setTPX(self.temperature+self.temperature*temp_del,
#                       self.pressure+self.pressure*pres_del)
        
        data,trash = self.run() #Ignore trash, just temp storage for empty kinetic sens array
        #print(data)
        
        #data = sim.Simulation.sensitivity_adjustment(self,temp_del,pres_del,spec_pair)
        self.temperatures=temptemp
        self.pressure=temppres
        self.conditions=tempcond
        self.kineticSens = kin_temp
        
        
        return data
    
    def species_adjustment(self,spec_del:float=0.0):
        inert_species=['Ar','AR','HE','He','Kr','KR',
                       'Xe','XE','NE','Ne']
        
        '''
        Creates tuples of specie that need to be perturbed and the
        percent value by which to perturb its mole fraction 
        '''
        # gets the mole fraction and the species which are going to be 
        #perturbed in order to run a sensitivity calculation 
        data = []
        for x in self.conditions.keys():
            if x not in inert_species:
                data.append(self.sensitivity_adjustment(spec_pair=(x,spec_del)))

        return data
    
    def importExperimentalData(self,csvFileList):
        print('Importing jsr data the following csv files...') 
        print(csvFileList)
        experimentalData = [pd.read_csv(csv) for csv in csvFileList]
        experimentalData = [experimentalData[x].dropna(how='any') for x in range(len(experimentalData))]
        experimentalData = [experimentalData[x].apply(pd.to_numeric, errors = 'coerce').dropna() for x in range(len(experimentalData))]
        for x in range(len(experimentalData)):
            experimentalData[x] = experimentalData[x][~(experimentalData[x][experimentalData[x].columns[1]] < 0)]
        self.experimentalData = experimentalData
        return experimentalData
    
    def map_and_interp_ksens(self,temp_history=None):
        A = self.kineticSensitivities
        N = np.zeros(A.shape)
        Ea = np.zeros(A.shape)
        for i in range(0,A.shape[2]):
            sheetA = A[:,:,i] #sheet for specific observable
            for x,column in enumerate(sheetA.T):
                N[:,x,i]= np.multiply(column,np.log(self.timeHistories[0]['temperature'])) if temp_history is None else np.multiply(column,np.log(temp_history['temperature']))
                #not sure if this mapping is correct, check with burke and also update absorption mapping
                #to_mult_ea = np.divide(-1,np.multiply(1/ct.gas_constant,self.timeHistories[0]['temperature'])) if time_history is None else np.divide(-1,np.multiply(ct.gas_constant,time_history['temperature']))
                to_mult_ea = np.divide(-1,np.multiply(1,self.timeHistories[0]['temperature'])) if temp_history is None else np.divide(-1,np.multiply(1,temp_history['temperature']))
                Ea[:,x,i]= np.multiply(column,to_mult_ea)
        #print(np.shape(A))
        tempA=[]
        tempn=[]
        tempEa=[]
        for i in range(0,A.shape[2]):
            tempA.append(A[:,:,i])
            tempn.append(N[:,:,i])
            tempEa.append(Ea[:,:,i])
        A=tempA
        N=tempn
        Ea=tempEa
        return {'A':A,
                'N':N,
                'Ea':Ea}
    def sensitivityCalculation(self,originalValues,newValues,thingToFindSensitivtyOf,dk=.01):
        if isinstance(originalValues,pd.DataFrame) and isinstance(newValues,pd.DataFrame):
            
            #newValues.columns = thingToFindSensitivtyOf
            
            newValues = newValues.applymap(np.log)
            originalValues = originalValues.applymap(np.log)
            #tab
            
            sensitivity = (newValues.subtract(originalValues)/dk)
            return sensitivity
        else:
            print("Error: wrong datatype, both must be pandas data frames")
            return -1
        
        