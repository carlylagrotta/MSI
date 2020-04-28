import cantera as ct
from .. import simulation as sim
from ...cti_core import cti_processor as ctp
import pandas as pd
import numpy as np
import time
import copy
import re
import MSI.simulations.instruments.shock_tube as st
import time
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

class ignition_delay(sim.Simulation):
    
    
    def __init__(self,pressure:float,temperature:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,
                 mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={}, save_timeHistories:int=0,
                 log_file=True,log_name='log.txt',timeshift:float=0.0,initialTime:float=0.0,
                 finalTime:float=1.0,target:str='temperature',target_type:str='max derivative',n_processors:int=2):
        
        
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)
        self.n_processors=n_processors
        self.pressure=pressure
        self.temperature=temperature
        self.observables=observables
        self.kineticSens=kineticSens
        self.physicalSens=physicalSens
        self.conditions=conditions
        self.cti_path=cti_path
        self.thermalBoundary = thermalBoundary
        self.mechanicalBoundary=mechanicalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        #self.energycon='off'
        self.timeshift=timeshift
        self.timeHistory = None
        self.experimentalData = None
        self.initialTime=initialTime
        self.finalTime=finalTime
        self.log_name=log_name
        self.log_file=log_file
        #self.yaml_file=yaml_file
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
        self.target=target
        self.target_type=target_type
    
    
    def run_shocktube(self,args=[],temp_proc=None):
        
        if not args:
            shock_tube = st.shockTube(pressure =self.pressure,
                         temperature = self.temperature,
                         observables = self.observables,
                         kineticSens = self.kineticSens,
                         physicalSens = 0,
                         conditions = self.conditions,
                         initialTime = self.initialTime,
                         finalTime = self.finalTime,
                         thermalBoundary = self.thermalBoundary,
                         mechanicalBoundary = self.mechanicalBoundary,
                         processor = self.processor,
                         save_timeHistories = 1,
                         save_physSensHistories = 0,
                         moleFractionObservables = self.moleFractionObservables,
                         concentrationObservables = self.concentrationObservables,
                         fullParsedYamlFile = self.fullParsedYamlFile,
                         time_shift_value = self.timeshift)
        elif args:
            shock_tube = st.shockTube(pressure =args[0],
                         temperature = args[1],
                         observables = args[2],
                         kineticSens = args[3],
                         physicalSens = 0,
                         conditions = args[4],
                         initialTime = args[5],
                         finalTime = args[6],
                         thermalBoundary = args[7],
                         mechanicalBoundary = args[8],
                         processor = temp_proc,
                         save_timeHistories = 1,
                         save_physSensHistories = 0,
                         moleFractionObservables = args[9],
                         concentrationObservables = args[10],
                         fullParsedYamlFile = args[11],
                         time_shift_value = args[12])
        
        shock_tube.run()
        return shock_tube
    
    
    def run_single(self):
        
        tic=time.time()
        s=self.run_shocktube()
        
        self.timehistory=copy.deepcopy(s.timeHistory)
        delay=None
        if re.match('[Mm]ax [Dd]erivative',self.target_type):
            if re.match('[Tt]emperature',self.target):
                
                delay=self.ig_dTdt(self.timehistory)
                
            elif re.match('[Pp]ressure',self.target_type):
                
                delay=self.ig_dPdt(self.timehistory)
                
            elif not re.match('[Tt]emperature',self.target) and not re.match('[Pp]ressure',self.target_type):
                
                delay=self.ig_dXdt(self.timehistory,self.target)
                
        sens=[]        
        if self.kineticSens:   
                
                sens=self.BFM(delay)
                #sens=self.BFM_pool(delay,self.n_processors)
            
        
        toc=time.time()
        print('Simulation took '+str(round(toc-tic,9))+' seconds.')
        return delay,sens
            
        
        
                
        
    def ig_dTdt(self,data):
        
        delay=[]
        tt=data['time'].values
        TT=data['temperature'].values
        #PP=data['pressure'].values
        dTdt=np.zeros(np.array(tt).shape,np.float)
        dTdt[0:-1]=np.diff(TT)/np.diff(tt)
        dTdt[-1]=(TT[-1] - TT[-2])/(tt[-1] - tt[-2])
        delay=tt[np.argmax(dTdt)]
        return delay
        
    def ig_dPdt(self,data):
        
        delay=[]        
        tt=data['time'].values
        PP=data['pressure'].values
        dPdt=np.zeros(np.array(tt).shape,np.float)
        dPdt[0:-1]=np.diff(PP)/np.diff(tt)
        dPdt[-1]=(PP[-1] - PP[-2])/(tt[-1] - tt[-2])
        delay=tt[np.argmax(dPdt)]
        return delay
        
    def ig_dXdt(self,data,target):
        
        delay=[]
        tt=data['time'].values
        XX=data[target].values
        dXdt=np.zeros(np.array(tt).shape,np.float)
        dXdt[0:-1]=np.diff(XX)/np.diff(tt)
        dXdt[-1]=(XX[-1] - XX[-2])/(tt[-1] - tt[-2])
        delay=tt[np.argmax(dXdt)]
        return delay
    
    def sensitivityCalculation(self,originalValues,newValues,dk=.01):
        sensitivity=(np.log(newValues)-np.log(originalValues))/dk
                           
        return sensitivity
        
    def BFM(self,nominal):
         
         sens=np.zeros(self.processor.solution.n_reactions)
         if re.match('[Mm]ax [Dd]erivative',self.target_type):
            if re.match('[Tt]emperature',self.target):
                
                for i in range(self.processor.solution.n_reactions):
                    print('Solving kinetic sensitivity for reaction '+str(i+1))
                    self.processor.solution.set_multiplier(1+self.dk,i)
                    temp_history=copy.deepcopy(self.run_shocktube().timeHistory)                    
                    delay=self.ig_dTdt(temp_history)
                    sens[i]=self.sensitivityCalculation(nominal,delay,self.dk)
                    self.processor.solution.set_multiplier(1.0,i)
                    
            elif re.match('[Pp]ressure',self.target_type):
                
                for i in range(self.processor.solution.n_reactions):
                    print('Solving kinetic sensitivity for reaction '+str(i+1))
                    self.processor.solution.set_multiplier(1+self.dk,i)
                    temp_history=copy.deepcopy(self.run_shocktube().timeHistory)                    
                    delay=self.ig_dPdt(temp_history)
                    sens[i]=self.sensitivityCalculation(nominal,delay,self.dk)
                    self.processor.solution.set_multiplier(1.0,i)
                    
            elif not re.match('[Tt]emperature',self.target) and not re.match('[Pp]ressure',self.target_type):
                
                for i in range(self.processor.solution.n_reactions):
                    print('Solving kinetic sensitivity for reaction '+str(i+1))
                    self.processor.solution.set_multiplier(1+self.dk,i)
                    temp_history=copy.deepcopy(self.run_shocktube().timeHistory)                    
                    delay=self.ig_dXdt(temp_history,self.target)
                    sens[i]=self.sensitivityCalculation(nominal,delay,self.dk)
                    self.processor.solution.set_multiplier(1.0,i)
         return sens
        
    def BFM_pool(self,nominal,n_procs):
         info=[]
         sens=[]
         args=(self.pressure,self.temperature,self.observables,self.kineticSens,self.conditions,self.initialTime,self.finalTime,
               self.thermalBoundary,self.mechanicalBoundary,self.moleFractionObservables,self.concentrationObservables,
               self.fullParsedYamlFile,self.timeshift)
         for i in range(self.processor.solution.n_reactions):
             info.append([self.dk,self.processor.cti_path,nominal,i,self.target_type,args,self.target])
         def solver(info):
             if re.match('[Mm]ax [Dd]erivative',info[4]):
                 if re.match('[Tt]emperature',info[6]):
             
                     print('Solving kinetic sensitivity for reaction '+str(info[3]+1))
                     temp_processor=ctp.Processor(info[1])
                     temp_processor.solution.set_multiplier(1+info[0],info[3])
                     temp_history=copy.deepcopy(self.run_shocktube(info[5],temp_processor).timeHistory)                    
                     delay=self.ig_dTdt(temp_history)
                     sens=self.sensitivityCalculation(nominal,delay,info[0])
                     
                 elif re.match('[Pp]ressure',info[6]):
             
                     print('Solving kinetic sensitivity for reaction '+str(info[3]+1))
                     temp_processor=ctp.Processor(info[1])
                     temp_processor.solution.set_multiplier(1+info[0],info[3])
                     temp_history=copy.deepcopy(self.run_shocktube(info[5],temp_processor).timeHistory)                    
                     delay=self.ig_dPdt(temp_history)
                     sens=self.sensitivityCalculation(nominal,delay,info[0])
                     
                 elif not re.match('[Tt]emperature',info[6]) and not re.match('[Pp]ressure',info[6]):
             
                     print('Solving kinetic sensitivity for reaction '+str(info[3]+1))
                     temp_processor=ctp.Processor(info[1])
                     temp_processor.solution.set_multiplier(1+info[0],info[3])
                     temp_history=copy.deepcopy(self.run_shocktube(info[5],temp_processor).timeHistory)                    
                     delay=self.ig_dXdt(temp_history,info[6])
                     sens=self.sensitivityCalculation(nominal,delay,info[0])
             return sens 
                     
                     
                     
         pool = ThreadPool(n_procs) 
         #results = results+pool.map(solver,conditionsTups)
         sens=pool.map(solver,info)
         #sens=np.zeros(self.processor.solution.n_reactions)
         
         return sens
        
class ignition_delay_wrapper(sim.Simulation):
    
    def __init__(self,pressures,temperatures,observables:list,
                 kineticSens:int,physicalSens:int,conditions,thermalBoundary,
                 mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={}, save_timeHistories:int=0,
                 log_file=True,log_name='log.txt',timeshift:float=0.0,initialTime:float=0.0,
                 finalTime:float=1.0,target:str='temperature',target_type:str='max derivative',n_processors:int=2):
        
        
        
        
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)
        self.n_processors=n_processors
        self.pressures=pressures
        self.temperatures=temperatures
        self.observables=observables
        self.kineticSens=kineticSens
        self.physicalSens=physicalSens
        self.conditions=conditions
        self.cti_path=cti_path
        self.thermalBoundary = thermalBoundary
        self.mechanicalBoundary=mechanicalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        #self.energycon='off'
        self.timeshift=timeshift
        self.timeHistory = None
        self.experimentalData = None
        self.initialTime=initialTime
        self.finalTime=finalTime
        self.log_name=log_name
        self.log_file=log_file
        #self.yaml_file=yaml_file
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
        self.target=target
        self.target_type=target_type
        
        
        
    def run(self):
        
        
        
        solution=[]
        ksens=[]
        ksens_1stIter=False
        for i in range(len(self.temperatures)):
            for j in range(len(self.pressures)):
                for k in range(len(self.conditions)):
                    temp_ig=ignition_delay(self,pressure=self.pressures[j],
                                           temperature=self.temperatures[i],
                                           observables=self.observables,
                                           kineticSens=self.kineticSens,
                                           physicalSens=self.physicalSens,
                                           conditions=self.conditions[k],
                                           thermalBoundary=self.thermalBoundary,
                                           mechanicalBoundary=self.mechanicalBoundary,
                                           processor=self.processor,
                                           cti_path=self.cti_path, 
                                           save_physSensHistories=self.saave_physSensHistories,
                                           moleFractionObservables=self.moleFractionObservables,
                                           absorbanceObservables=self.absorbanceObservables,
                                           concentrationObservables=self.concentrationObservables,
                                           fullParsedYamlFile=self.fullParsedYamlFile, 
                                           save_timeHistories=self.save_timeHistories,
                                           log_file=self.log_file,
                                           log_name=self.log_name,
                                           timeshift=self.timeshift,
                                           initialTime=self.initialTime,
                                           finalTime=self.finalTime,
                                           target=self.target,
                                           target_type=self.target_type,
                                           n_processors=self.n_processors)
            
                    a,b=temp_ig.run_single()
            
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        