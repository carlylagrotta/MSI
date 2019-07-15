from ..cti_core import cti_processor as ctp
import copy

class Simulation(object):
    pasc_to_atm = 101325
    def __init__(self,pressure:float,temperature:float,observables:list,kineticSens:int,physicalSens:int
            ,conditions:dict,processor:ctp.Processor=None,cti_path = ""):
        '''
        Input:
            - pressure = float, pressure in [atm]
            - temperature = float, temperature in [K]
            - observables = list, species which sensitivity analysis is performed for
            - kineticSen = integer, 0 for off, 1 for on 
            - physicalSens = integer, 0 for off, 1 for on 
            - processor = ctp.Processor
            - cti_path = string, path to cti file, will construct an internal processor
            
        '''
        #set up processor and initialize all variables  common to every simulation  
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)
        self.pressure = pressure
        self.temperature = temperature
        self.observables = observables
        self.kineticSens = kineticSens
        self.physicalSens = physicalSens
        self.conditions = conditions
        self.dk = []        
    def setTPX(self,temperature:float=-1,pressure:float=-1,conditions_perturb:dict={},reset_value={}):
        '''
        Set solution object for a simulation
        '''
        #set the temperature, pressure and species mole fractions for the simulation
        if temperature== -1:
            temperature = self.temperature
        if pressure == -1:
            pressure = self.pressure
        if conditions_perturb == {}:
            new_conditions = self.conditions
        else:
            #make copy of the original mole fractions so they can be changed 
            #for sensitivity analysis but the original ones are saved 
            conditions_copy = copy.deepcopy(self.conditions)
            for x in conditions_perturb.keys():
                if x != '':
                    
                    #print(conditions_perturb[x],'these are the perturbed conditions')
                    #conditions_copy[x] = ((conditions_copy[x]+conditions_perturb[x])*(1-conditions_copy[x]))/((1-conditions_copy[x])-(conditions_copy[x]+conditions_perturb[x]))
                    conditions_copy[x] = ((conditions_copy[x]+conditions_perturb[x])*(1-conditions_copy[x]))/(1 - (conditions_copy[x]+conditions_perturb[x]))
            new_conditions = conditions_copy
            #print(new_conditions,'these are the new conditions')
        
        self.processor.solution.TPX=temperature,pressure*self.pasc_to_atm, new_conditions
        #print(self.processor.solution.TPX)
        
           
    #always overwritten since each simulation is very different
    def run(self):
        print("Error: Simulation class itself does not implement the run method, please run a child class")
        '''
        Run simulation
        '''
        # run the simulation 
    def sensitivity_adjustment(self,temp_del:float=0.0,
                               pres_del:float=0.0,
                               spec_pair:(str,float)=('',0.0)):
        '''
          Passes the Perturbed observable to the setTPX function. Temperature and pressure 
        are passed and set directly species need to go through an additional step in the 
        setTPX function. 
        '''
        if spec_pair[0] != '':
            
           # calculates the value to run the sensitivity calculation at 
           #for physical observables 
           self.setTPX(self.temperature+self.temperature*temp_del,
                   self.pressure+self.pressure*pres_del,
                   {spec_pair[0]:self.conditions[spec_pair[0]]*spec_pair[1]})
          
           
        else:
           self.setTPX(self.temperature+self.temperature*temp_del,
                       self.pressure+self.pressure*pres_del)
        
        data = self.run()

            
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

