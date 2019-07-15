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
<<<<<<< HEAD
=======
import copy
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0

class JSR_steadystate(sim.Simulation):
    
    '''child class of sim.Simulaton.  Inherits all attributes and methods including __init__().  
    Also has internal init due to data requirements'''
    
    
    
    
    
    def __init__(self,pressure:float,temperature:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={},residence_time:float=1.0,pvalveCoefficient:float=0.01,maxpRise:float=0.001,save_timeHistories:int=0):
        
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
        self.mechanicalBoundary = mechanicalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        self.pvalveCoefficient=pvalveCoefficient
        self.maxPrise=maxpRise
        self.energycon='off'
        self.residence_time=residence_time
        self.timeHistory = None
        self.experimentalData = None
        
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
        self.rtol=1e-6
        self.atol=1e-6
        self.solution=None
<<<<<<< HEAD
        
=======
        
        
    def set_geometry(self,volume=0.1):
        self.reactor_volume=volume
        
        
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
    def printVars(self):
        print()
        
    def settingJSRConditions(self):
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
    
<<<<<<< HEAD
    
    
    def run(self,):
        
        gas=processor.solution.P
        reactorPressure=gas
        pressureValveCoefficient=pvalveCoefficient
        maxPressureRiseAllowed=maxPrise
        
        
        #Build the system components for JSR
        fuelAirMixtureTank=ct.Reservoir(gas)
        exhaust=ct.Reservoir(gas)
        
        stirredReactor=ct.IdealGasReactor(gas,energy=energycon,volume=volume)
        massFlowController=ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                 downstream=stirredReactor,mdot=stirredReactor.mass/restime)
        pressureRegulator=ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)
        reactorNetwork=ct.ReactorNet([stirredReactor])
        
        if bool(self.moleFractionObservables) and self.kineticSens==1:
            for i in range(gas.n_reactions):
                stirredReactor.add_sensitivity_reaction(i)
            
        if self.kineticSens and bool(self.moleFractionObservables)==False:
            except:
                print('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        if self.physicalSens==1 and bool(self.moleFractionObservables)==False:
            except:
=======
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
        reactorPressure=gas.P
        self.reactorPressure=self.processor.solution.P
        pressureValveCoefficient=self.pvalveCoefficient
        maxPressureRiseAllowed=self.maxPrise
        
        print(maxPressureRiseAllowed,self.reactorPressure,pressureValveCoefficient)
        #Build the system components for JSR
        fuelAirMixtureTank=ct.Reservoir(self.processor.solution)
        exhaust=ct.Reservoir(self.processor.solution)
        
        stirredReactor=ct.IdealGasReactor(self.processor.solution,energy=self.energycon,volume=self.reactor_volume)
        massFlowController=ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                 downstream=stirredReactor,mdot=stirredReactor.mass/self.residence_time)
        pressureRegulator=ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)
        reactorNetwork=ct.ReactorNet([stirredReactor])
        
        if bool(self.observables) and self.kineticSens==1:
            for i in range(gas.n_reactions):
                stirredReactor.add_sensitivity_reaction(i)
            
        if self.kineticSens and bool(self.observables)==False:
            #except:
                print('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        if self.physicalSens==1 and bool(self.observables)==False:
            #except:
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
                print('Please supply a non-empty list of observables for sensitivity analysis or set physical_sens=0')
        
        # now compile a list of all variables for which we will store data
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        columnNames = ['pressure'] + columnNames

        # use the above list to create a DataFrame
        timeHistory = pd.DataFrame(columns=columnNames)

        # Start the stopwatch
        tic = time.time()
        
        
        #Establish a matrix to hold sensitivities for kinetic parameters, along with tolerances
<<<<<<< HEAD
        if self.kineticSens==1 and bool(self.moleFractionObservables):
            #senscolumnNames = ['Reaction']+observables     
            senscolumnNames = self.moleFractionObservables
            #sensArray = pd.DataFrame(columns=senscolumnNames)
            #senstempArray = np.zeros((gas.n_reactions,len(observables)))
            dfs = [pd.DataFrame() for x in range(len(self.moleFractionObservables))]
            tempArray = [np.zeros(self.processor.solution.n_reactions) for x in range(len(self.moleFractionObservables))]
=======
        if self.kineticSens==1 and bool(self.observables):
            #senscolumnNames = ['Reaction']+observables     
            senscolumnNames = self.observables
            #sensArray = pd.DataFrame(columns=senscolumnNames)
            #senstempArray = np.zeros((gas.n_reactions,len(observables)))
            dfs = [pd.DataFrame() for x in range(len(self.observables))]
            tempArray = [np.zeros(self.processor.solution.n_reactions) for x in range(len(self.observables))]
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
            
            reactorNetwork.rtol_sensitivity = self.rtol
            reactorNetwork.atol_sensitivity = self.atol        
        
        
        reactorNetwork.advance_to_steady_state()
        final_pressure=stirredReactor.thermo.P
<<<<<<< HEAD
        sens=reactorNetwork.sensitivites()
        if self.kineticSens==1 and bool(self.moleFractionObservables):
            for k in range(len(self.moleFractionObservables)):
                dfs[k] = dfs[k].append(((pd.DataFrame(sens[k,:])).transpose()),ignore_index=True)        
=======
        sens=reactorNetwork.sensitivities()
        #print(sens)
        if self.kineticSens==1 and bool(self.observables):
            #print((pd.DataFrame(sens[0,:])).transpose())
            #test=pd.concat([pd.DataFrame(),pd.DataFrame(sens[0,:]).transpose()])
            #print(test)
            for k in range(len(self.observables)):
                dfs[k] = dfs[k].append(((pd.DataFrame(sens[k,:])).transpose()),ignore_index=True)
                #dfs[k]=pd.concat([dfs[k],pd.DataFrame(sens[k,:]).transpose()])
                #dfs[k]=pd.DataFrame(sens[k,:]).transpose()
            #print(dfs)  
        toc = time.time()
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
        print('Simulation Took {:3.2f}s to compute'.format(toc-tic))
   
        columnNames = []
        #Store solution to a solution array
        #for l in np.arange(stirredReactor.n_vars):
            #columnNames.append(stirredReactor.component_name(l))
        columnNames=[stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        #state=stirredReactor.get_state()
        state=np.hstack([stirredReactor.mass, 
                   stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
        data=pd.DataFrame(state).transpose()
        data.columns=columnNames
        pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
<<<<<<< HEAD
        if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
            except:
=======
        if(abs(pressureDifferential/self.reactorPressure) > maxPressureRiseAllowed):
            #except:
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
                print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
        
        if self.kineticSens==1:
            numpyMatrixsksens = [dfs[dataframe].values for dataframe in range(len(dfs))]
            self.kineticSensitivities = np.dstack(numpyMatrixsksens)
<<<<<<< HEAD
            self.solution=data
            return (self.solution,self.kineticSensitivities)
        else:
            return self.solution
=======
            #print(np.shape(self.kineticSensitivities))
            self.solution=data
            return (self.solution,self.kineticSensitivities)
        else:
            self.solution=data
            return (self.solution,[])
        
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
        
        
        
        
<<<<<<< HEAD
        
        
        
        
        
=======
#class JSR_multiTemp_steadystate2(sim.Simulation):
#        
#    def __init__(self,volume:float,pressure:float,temperatures:list,observables:list,
#                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,mechanicalBoundary,
#                 processor:ctp.Processor=None,cti_path="", 
#                 save_physSensHistories=0,moleFractionObservables:list=[],
#                 absorbanceObservables:list=[],concentrationObservables:list=[],
#                 fullParsedYamlFile:dict={},residence_time:float=1.0,pvalveCoefficient:float=0.01,maxpRise:float=0.001):
#        
##    sim.Simulation.__init__(self,pressure,temperature,observables,kineticSens,physicalSens,
##                                conditions,processor,cti_path)
#        self.volume=volume
#        self.temperatures=temperatures
#        self.JSR_objects=[]
#        for i in range(len(self.temperatures)):
#            a=[]
#            a=copy.deepcopy(JSR_steadystate(pressure,self.temperatures[i],observables,
#                 kineticSens,physicalSens,conditions,thermalBoundary,mechanicalBoundary,
#                 processor,cti_path, 
#                 save_physSensHistories,moleFractionObservables,
#                 absorbanceObservables,concentrationObservables,
#                 fullParsedYamlFile,residence_time,pvalveCoefficient,maxpRise))
#            self.JSR_objects.append(a)
#        for i in self.JSR_objects:
#            print(i.temperature)
#            
#            
#    def run_multitemp(self):
#        
#        for i in range(len(self.JSR_objects)):
#            self.JSR_objects[i].set_geometry(self.volume)
#            
#            self.JSR_objects[i].run()
#        
#        
#        
class JSR_multiTemp_steadystate(sim.Simulation):
        
    def __init__(self,volume:float,pressure:float,temperatures:list,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],save_timeHistories:int=0,
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={},residence_time:float=1.0,pvalveCoefficient:float=0.01,maxpRise:float=0.001):
        
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
        self.volume=volume
        self.temperatures=temperatures
        self.JSR_objects=[]
        self.pressure=pressure
        self.observables=observables
        self.kineticSens=kineticSens
        self.physicalSens=physicalSens
        self.conditions=conditions
        self.thermalBoundary=thermalBoundary
        self.mechanicalBoundary = mechanicalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        self.pvalveCoefficient=pvalveCoefficient
        self.maxPrise=maxpRise
        self.energycon='off'
        self.residence_time=residence_time
        
               
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
            temp_jsr=JSR_steadystate(pressure=self.pressure,
                                     temperature=self.temperatures[i],
                                     observables=self.observables,
                                     kineticSens=self.kineticSens,
                                     physicalSens=self.physicalSens,
                                     conditions=self.conditions,
                                     thermalBoundary=self.thermalBoundary,
                                     mechanicalBoundary=self.mechanicalBoundary,
                                     processor=self.processor,
                                     save_physSensHistories=self.save_physSensHistories)
            temp_jsr.set_geometry(volume=self.volume)
            a,b=temp_jsr.run_single()
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
        return (solution,ksens)
        
        
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
            delxj=spec_pair[1]
            self.conditions[spec_pair[0]]=np.divide(np.multiply(xj+delxj,1-xj),1-xj-delxj)
#           self.setTPX(self.temperature+self.temperature*temp_del,
#                   self.pressure+self.pressure*pres_del,
#                   {spec_pair[0]:self.conditions[spec_pair[0]]*spec_pair[1]})
          
           
        else:
           self.temperatures=np.array(self.temperatures)+temp_del*np.array(self.temperatures)
           self.pressure=self.pressure+pres_del*self.pressure
           
#           self.setTPX(self.temperature+self.temperature*temp_del,
#                       self.pressure+self.pressure*pres_del)
        
        data,trash = self.run() #Ignore trash, just temp storage for empty kinetic sens array
        
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
        data = ''
        for x in self.conditions.keys():
            if x not in inert_species:
                data =  self.sensitivity_adjustment(spec_pair=(x,spec_del))

        return data
>>>>>>> f1481addf2130afeceaf2e09ed855c3bf1a0fdf0
        
        