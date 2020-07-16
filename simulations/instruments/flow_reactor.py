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

class flow_reactor(sim.Simulation):
    
    
    def __init__(self,pressure:float,temperature:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,
                 mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 concentrationObservables:list=[],
                 fullParsedYamlFile:dict={}, save_timeHistories:int=0,
                 log_file=True,log_name='log.txt',timeshift:float=0.0,initialTime:float=0.0,
                 residenceTime:float=1.0):
        
        
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
        self.mechanicalBoundary=mechanicalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        #self.energycon='off'
        self.timeshift=timeshift
        self.timeHistory = None
        self.experimentalData = None
        self.initialTime=initialTime
        self.residenceTime = residenceTime
        self.finalTime = self.timeshift + self.residenceTime
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

    
    

    
    def run_shocktube(self,ksens_marker=1 ,psens_marker=1):
        if ksens_marker ==0 and psens_marker==0:
            shock_tube = st.shockTube(pressure =self.pressure,
                          temperature = self.temperature,
                          observables = self.observables,
                          kineticSens = 0,
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
            shock_tube.run()

            return shock_tube
        
        elif ksens_marker ==1 and psens_marker==0:
            shock_tube = st.shockTube(pressure =self.pressure,
                              temperature = self.temperature,
                              observables = self.observables,
                              kineticSens = 1,
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
            shock_tube.run()
            return shock_tube
        
        elif ksens_marker ==0 and psens_marker==1:
            shock_tube = st.shockTube(pressure =self.pressure,
                              temperature = self.temperature,
                              observables = self.observables,
                              kineticSens = 0,
                              physicalSens = 1,
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
            shock_tube.run()
            return shock_tube

        
        elif ksens_marker ==1 and psens_marker==1:
            shock_tube = st.shockTube(pressure =self.pressure,
                              temperature = self.temperature,
                              observables = self.observables,
                              kineticSens = 1,
                              physicalSens = 1,
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
            shock_tube.run()
            return shock_tube
    
    def run_single(self,ksens_marker=1,psens_marker=1):
        
               
               
        if self.kineticSens: 
            
            s = self.run_shocktube(ksens_marker=1,psens_marker=0)
            self.timehistory=copy.deepcopy(s.timeHistory)
            res_time_measurment=None
            res_time_measurment,index = self.get_res_time_data(self.timehistory,self.finalTime)   
            
            #print(s.kineticSensitivities.shape)
            #ksens = s.kineticSensitivities[-1,:,:]
            ksens,temp_arrays = self.get_ksens_at_res_time(s.kineticSensitivities,self.timehistory['time'],self.finalTime)
            #ksens = s.kineticSensitivities[index,:,:]
            
            #xdim = s.kineticSensitivities.shape[0]
            #ydim = s.kineticSensitivities.shape[1]
            #zdim = s.kineticSensitivities.shape[2]


            #ksens = ksens.reshape((1,ydim,zdim))
            
            
            self.kineticSensitivities = ksens
            
        elif ksens_marker==0 and psens_marker==1:
            
            s = self.run_shocktube(ksens_marker=0,psens_marker=1)
            
            self.timehistory=copy.deepcopy(s.timeHistory)
            
            res_time_measurment=None
            res_time_measurment,index = self.get_res_time_data(self.timehistory,self.finalTime)             
            
        else:
            s = self.run_shocktube(ksens_marker=0,psens_marker=0)
            self.timehistory=copy.deepcopy(s.timeHistory)
            res_time_measurment=None
            res_time_measurment,index = self.get_res_time_data(self.timehistory,self.finalTime) 
            

        if self.kineticSens:  
            
            return res_time_measurment,self.kineticSensitivities,self.timehistory,temp_arrays
        else:
            return res_time_measurment,[],None,None
            
        
        
    def get_ksens_at_res_time(self,ksens,time_array,res_time):
        ksens_array = []
        temp_arrays = []
        for sheet in range(ksens.shape[2]):
            temp = ksens[:,:,sheet]
            time=time_array.values
            time=time.reshape((time.shape[0],1))
            temp_with_time = np.hstack((time,temp))
            df =copy.deepcopy(temp_with_time)
            df = pd.DataFrame(temp_with_time)
            df=df.rename(columns = {0:'time'})
            temp_arrays.append(df)
            df.loc[-1, 'time'] = float(res_time)

            df = df.sort_values('time').reset_index(drop=True)
            
            df = df.interpolate()
            res_time_k_sens_data = df.iloc[(df['time']-res_time).abs().argsort()[:1]]
            res_time_k_sens_data = res_time_k_sens_data.reset_index(drop=True)
            res_time_k_sens_data = res_time_k_sens_data.drop(columns="time")
            res_time_k_sens_data = res_time_k_sens_data.to_numpy()
            
    
            res_time_k_sens_data = res_time_k_sens_data.reshape((res_time_k_sens_data.shape[0],res_time_k_sens_data.shape[1],1))
            ksens_array.append(res_time_k_sens_data)
            

        ksens_array = np.dstack((ksens_array))
        
        return ksens_array,temp_arrays


           
        
    def get_res_time_data(self,data,res_time):        
        #res_time_data = data.tail(1)
        #res_time_data = res_time_data.reset_index(drop=True)
        
        #print(res_time)
        
        #print(data.tail(1))
        #reset index
        df = copy.deepcopy(data)
        
        df.loc[-1, 'time'] = float(res_time)
        df = df.sort_values('time').reset_index(drop=True)
        df = df.interpolate()
        
        res_time_data = df.iloc[(df['time']-res_time).abs().argsort()[:1]]
        res_time_data = res_time_data.reset_index(drop=True)
        
        index = df.iloc[(df['time']-res_time).abs().argsort()[:1]].index.values[0]
        
        return res_time_data,index

    
    def sensitivityCalculation(self,originalValues,newValues,dk=.01):

        sensitivity=(np.log(newValues)-np.log(originalValues))/dk
                           
        return sensitivity
    

        
class flow_reactor_wrapper(sim.Simulation):
    
    def __init__(self,pressure:float,temperatures:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,
                 mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 concentrationObservables:list=[],
                 fullParsedYamlFile:dict={}, save_timeHistories:int=0,
                 timeshifts:list=[],initialTime:float=0.0,
                 residenceTimes:list=1.0):
        
        
        
        
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)
        self.pressure=pressure
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
        self.fullParsedYamlFile =  fullParsedYamlFile
        #self.energycon='off'
        self.timeshifts=timeshifts
        self.timeHistory = None
        self.experimentalData = None
        self.initialTime=initialTime
        self.residenceTimes = residenceTimes
        self.finalTimes = list(np.array(self.timeshifts) + np.array(self.residenceTimes))
        self.save_physSensHistories = save_physSensHistories
        self.save_timeHistories = save_timeHistories
        
        #self.yaml_file=yaml_file
        if save_timeHistories == 1:
            self.timeHistories=[]
            self.fullTimeHistories=[]
            self.temp_arrays=[]

        else:
            self.timeHistories=None
        if save_physSensHistories == 1:
            self.physSensHistories = []
        #self.setTPX()
        self.dk = [0]
        self.solution=None
        
        
        
    def run(self,ksens_marker=1,psens_marker=1):
        
        
        
        solution=[]
        ksens=[]
        ksens_1stIter=False
        #print(self.conditions)
        for i in range(len(self.temperatures)):
            temp_flow=flow_reactor(pressure=self.pressure,
                                   temperature=self.temperatures[i],
                                   observables=self.observables,
                                   kineticSens=self.kineticSens,
                                   physicalSens=self.physicalSens,
                                   conditions=self.conditions,
                                   thermalBoundary=self.thermalBoundary,
                                   mechanicalBoundary=self.mechanicalBoundary,
                                   processor=self.processor,
                                   cti_path=self.cti_path, 
                                   save_physSensHistories=self.save_physSensHistories,
                                   moleFractionObservables=self.moleFractionObservables,
                                   concentrationObservables=self.concentrationObservables,
                                   fullParsedYamlFile=self.fullParsedYamlFile, 
                                   save_timeHistories=self.save_timeHistories,
                                   timeshift=self.timeshifts[i],
                                   initialTime=self.initialTime,
                                   residenceTime=self.residenceTimes[i])
            
            #res_time_data,k_sens=temp_flow.run_single(ksens=self.kineticSens,psens=self.physicalSens)
            
            
            res_time_data,k_sens,fullTimeHistory,temp_array=temp_flow.run_single(ksens_marker=ksens_marker,psens_marker=psens_marker)
            if self.kineticSens==1:
                self.fullTimeHistories.append(fullTimeHistory)
                self.temp_arrays.append(temp_array)
                
            
            
            temp=[]
            temp1=[]
            temp=copy.deepcopy(res_time_data)
            #print(temp)
            temp1=copy.deepcopy(k_sens)
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
        #print(self.timeHistories)
        #print(solution)

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
            #self.pressures=np.array(self.pressures)+pres_del*np.array(self.pressures)
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
           #self.pressures=np.array(self.pressure)+pres_del*np.array(self.pressure)
           self.pressure=self.pressure+pres_del*self.pressure
           #self.residence_time=self.residence_time+res_del*self.residence_time
           
#           self.setTPX(self.temperature+self.temperature*temp_del,
#                       self.pressure+self.pressure*pres_del)
        
        data,trash = self.run(ksens_marker=0,psens_marker=1) #Ignore trash, just temp storage for empty kinetic sens array
        #print(data)
        
        #data = sim.Simulation.sensitivity_adjustment(self,temp_del,pres_del,spec_pair)
        self.temperatures=temptemp
        self.pressures=temppres
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
        print('Importing flow reactor data the following csv files...') 
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
        #print(self.kineticSensitivities.shape)
        #print(np.shape(A))
        N = np.zeros(A.shape)
        Ea = np.zeros(A.shape)
        for i in range(0,A.shape[2]):
            sheetA = A[:,:,i] #sheet for specific observable

            #print(sheetA)
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
    def sensitivityCalculation(self,originalValues,newValues,dk=.01):
        #print(originalValues)
        #print(newValues)
        if isinstance(originalValues,pd.DataFrame) and isinstance(newValues,pd.DataFrame) or isinstance(originalValues,pd.Series) and isinstance(newValues,pd.Series):
            if isinstance(originalValues,pd.Series) or isinstance(newValues,pd.Series):
                originalValues=originalValues.to_frame()
                newValues=newValues.to_frame()
            #newValues.columns = thingToFindSensitivtyOf
            
            newValues = newValues.applymap(np.log)
            originalValues = originalValues.applymap(np.log)
            #tab
            
            sensitivity = (newValues.subtract(originalValues)/dk)
            return sensitivity
        else:
            print("Error: wrong datatype, both must be pandas data frames")
            return -1
    

        
            
    def get_res_time_data(self,data,res_time):        
        df = copy.deepcopy(data)
        
        df.loc[-1, 'time'] = float(res_time)
        df = df.sort_values('time').reset_index(drop=True)
        df = df.interpolate()
        
        res_time_data = df.iloc[(df['time']-res_time).abs().argsort()[:1]]
        res_time_data = res_time_data.reset_index(drop=True)
        
        index = df.iloc[(df['time']-res_time).abs().argsort()[:1]].index.values[0]
        
        return res_time_data

    
    def calculate_time_shift_sensitivity(self,simulation,timeHistory,dk,finalTime):
        
        lst_obs = simulation.moleFractionObservables + simulation.concentrationObservables
        lst_obs = [i for i in lst_obs if i] 
        mean_times_of_experiments = []
            
        one_percent_of_average = 1e-8
                        
        
        original_time = timeHistory['time']
        new_time = original_time + one_percent_of_average
            

        #interpolate to the orignal time 

        interpolated_against_original_time = []
        for i,obs in enumerate(lst_obs):
            interpolated_original_observable_against_original_time = np.interp(original_time,new_time,timeHistory[lst_obs[i]])
            s1 = pd.Series(interpolated_original_observable_against_original_time,name=lst_obs[i])
            interpolated_against_original_time.append(s1)
        
        
        observables_interpolated_against_original_time_df = pd.concat(interpolated_against_original_time,axis=1)
        
        #calculate sensitivity
        
        calculated_sensitivity = []
        for i,obs in enumerate(lst_obs):
                      
           sens = (observables_interpolated_against_original_time_df[obs].apply(np.log) - timeHistory[obs].apply(np.log))/one_percent_of_average
           s1 = pd.Series(sens,name=lst_obs[i])
           calculated_sensitivity.append(s1)
            
        calculated_sensitivity.append(original_time)
        calculated_sensitivity_df = pd.concat(calculated_sensitivity,axis=1)
        
        

        time_shift_sensitivity = calculated_sensitivity_df
        #how to call this from the other class?
        time_shift_sensitivity = self.get_res_time_data(time_shift_sensitivity,finalTime)
        time_shift_sensitivity = time_shift_sensitivity.drop(columns="time")

        #self.time_shift_sensitivity = time_shift_sensitivity
        average_time=1
        self.average_time = average_time
        return time_shift_sensitivity       
        
        
        
        
        
        
        
        
        