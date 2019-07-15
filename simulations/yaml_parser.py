import yaml 
import shutil
import numpy as np
import copy

# subpackage for reading yaml files that describe simulations and absorbance data
class Parser(object):
    def __init__(self,original_experimental_conditions=None):
        self.original_experimental_conditions =  original_experimental_conditions
       

    #config is a dict containing the yaml information
    def load_to_obj(self, path:str = ''):
        with open(path) as f:
            config = yaml.load(f)
        return config
    
    def parse_shock_tube_obj(self,loaded_exp:dict={}, loaded_absorption:dict={}):
        simulation_type = loaded_exp['apparatus']['kind']
        pressure = loaded_exp['common-properties']['pressure']['value']
        temperature = loaded_exp['common-properties']['temperature']['value']
        mole_fractions = [((concentration['mole-fraction'])) for concentration in loaded_exp['common-properties']['composition']]
        mole_fractions = [float(elm) for elm in mole_fractions]
        species_names = [(species['species']) for species in loaded_exp['common-properties']['composition']]
        conditions = dict(zip(species_names,mole_fractions))
        thermal_boundary = loaded_exp['common-properties']['assumptions']['thermal-boundary']
        mechanical_boundary = loaded_exp['common-properties']['assumptions']['mechanical-boundary']
        
        mole_fraction_observables = [point['targets'][0]['name'] for point in loaded_exp['datapoints']['mole-fraction']]
        species_uncertainties = [uncert['relative-uncertainty'] for uncert in loaded_exp['common-properties']['composition']]
        species_uncertainties = [float(elm) for elm in species_uncertainties]
        species_uncertainties = dict(zip(species_names,species_uncertainties))
        
        
        concentration_observables = [datapoint['targets'][0]['name'] for datapoint in loaded_exp['datapoints']['concentration']]            
        observables = [x for x in (mole_fraction_observables + concentration_observables) if x is not None]
        
        initial_time = loaded_exp['common-properties']['time']['initial-time']['value']
        #eventually going to get this from a csv file 
        final_time = loaded_exp['common-properties']['time']['final-time']['value']
    
   
        mole_fraction_csv_files = [csvfile['csvfile'] for csvfile in loaded_exp['datapoints']['mole-fraction']]
        concentration_csv_files = [csvfile['csvfile'] for csvfile in loaded_exp['datapoints']['concentration']]
        path_length = loaded_exp['apparatus']['inner-diameter']['value']
        csv_files = [x for x in (mole_fraction_csv_files + concentration_csv_files) if x is not None]


        #importing unceratinty values 
        temp_relative_uncertainty = loaded_exp['common-properties']['temperature']['relative-uncertainty']
        temp_relative_uncertainty = float(temp_relative_uncertainty)
        pressure_relative_uncertainty = loaded_exp['common-properties']['pressure']['relative-uncertainty']
        pressure_relative_uncertainty = float(pressure_relative_uncertainty)
        time_shift_uncertainty = loaded_exp['common-properties']['time-shift']['absolute-uncertainty']['value']
        concentration_absolute_uncertainty = [point['targets'][0]['absolute-uncertainty'] for point in loaded_exp['datapoints']['concentration']]
        concentration_relative_uncertainity = [point['targets'][0]['relative-uncertainty'] for point in loaded_exp['datapoints']['concentration']]

        mole_fraction_absolute_uncertainty = [point['targets'][0]['absolute-uncertainty'] for point in loaded_exp['datapoints']['mole-fraction']]

        mole_fraction_relative_uncertainty = [point['targets'][0]['relative-uncertainty'] for point in loaded_exp['datapoints']['mole-fraction']]        

        if loaded_absorption == {}:
            return{
               'pressure':pressure,
               'temperature':temperature,
               'conditions':conditions,
               'speciesUncertaintys':species_uncertainties,
               'thermalBoundary':thermal_boundary,
               'mechanicalBoundary':mechanical_boundary,
               'moleFractionObservables':mole_fraction_observables,
               'concentrationObservables': concentration_observables,               
               'observables':observables,
               'initialTime':initial_time,
               'finalTime':final_time,
               'speciesNames':species_names,
               'pathLength':path_length,
               'MoleFractions':mole_fractions,
               'moleFractionCsvFiles':mole_fraction_csv_files,
               'concentrationCsvFiles':concentration_csv_files,
               'tempRelativeUncertainty':temp_relative_uncertainty,
               'pressureRelativeUncertainty': pressure_relative_uncertainty,
               'timeShiftUncertainty':time_shift_uncertainty,
               'concentrationAbsoluteUncertainty':concentration_absolute_uncertainty,
               'concentrationRelativeUncertainity':concentration_relative_uncertainity,
               'moleFractionAbsoluteUncertainty':mole_fraction_absolute_uncertainty,
               'moleFractionRelativeUncertainty':mole_fraction_relative_uncertainty,
               'csvFiles': csv_files,
               'simulationType':  simulation_type
           }
        
        else: #absorbtion file given
            absorbance_absolute_uncertainty = [point['absolute-uncertainty'] for point in loaded_exp['datapoints']['absorbance']]
            absorbance_relative_uncertainty = [point['relative-uncertainty'] for point in loaded_exp['datapoints']['absorbance']]
            #importing absorbance uncertainty 

            absorbance_csv_files = [csvfile['csvfile'] for csvfile in loaded_exp['datapoints']['absorbance']]
            absorbance_csv_wavelengths = [csvfile['wavelength']['value'] for csvfile in loaded_exp['datapoints']['absorbance']]
            absorption_observables = [species['species'] for species in loaded_absorption['Absorption-coefficients']]

            observables = [x for x in (mole_fraction_observables + concentration_observables + absorption_observables) if x is not None]


            uncertainty_parameter_ones = [[] for i in range(len(loaded_absorption['Absorption-coefficients']))]
            for uncertainty in range(len(loaded_absorption['Absorption-coefficients'])):
                temp = [wavelength['parameter-one']['absolute-uncertainty']['value'] for wavelength in loaded_absorption['Absorption-coefficients'][uncertainty]['wave-lengths']]
                uncertainty_parameter_ones[uncertainty] = temp
                
            uncertainty_parameter_twos = [[] for i in range(len(loaded_absorption['Absorption-coefficients']))]
            for uncertainty in range(len(loaded_absorption['Absorption-coefficients'])):
                temp = [wavelength['parameter-two']['absolute-uncertainty']['value'] for wavelength in loaded_absorption['Absorption-coefficients'][uncertainty]['wave-lengths']]
                uncertainty_parameter_twos[uncertainty] = temp        
                
                
            # add the function which will return the coupled paramters here 
            parameter_ones = []
            for p1 in range(len(loaded_absorption['Absorption-coefficients'])):
                temp = [wl['parameter-one']['value'] for wl in loaded_absorption['Absorption-coefficients'][p1]['wave-lengths']]
                parameter_ones.append(temp)

            parameter_twos = [] 
            for p2 in range(len(loaded_absorption['Absorption-coefficients'])):
                temp = [wl['parameter-two']['value'] for wl in loaded_absorption['Absorption-coefficients'][p2]['wave-lengths']]
                parameter_twos.append(temp)
                
            coupledCoefficients = [list(zip(parameter_ones[x],parameter_twos[x])) for x in range(len(parameter_ones))]
            functional_form = []
            for form in range(len(loaded_absorption['Absorption-coefficients'])):
                temp = [wl['functional-form'] for wl in loaded_absorption['Absorption-coefficients'][form]['wave-lengths']]
                functional_form.append(temp)

   
            return {
                   'pressure':pressure,
                   'temperature':temperature,
                   'conditions':conditions,
                   'thermalBoundary':thermal_boundary,
                   'mechanicalBoundary':mechanical_boundary,
                   'speciesNames': species_names,
                   'observables': observables,
                   'moleFractionObservables':mole_fraction_observables,
                   'concentrationObservables':concentration_observables, 
                   'absorbanceObservables':absorption_observables,
                   'initialTime': initial_time,
                   'finalTime':final_time,
                   'speciesNames': species_names,
                   'MoleFractions':mole_fractions,
                   'absorbanceCsvFiles': absorbance_csv_files,
                   'moleFractionCsvFiles':mole_fraction_csv_files,
                   'concentrationCsvFiles':concentration_csv_files,
                   'absorbanceCsvWavelengths': absorbance_csv_wavelengths,
                   'pathLength':path_length,
                   'tempRelativeUncertainty': temp_relative_uncertainty,
                   'pressureRelativeUncertainty': pressure_relative_uncertainty,
                   'speciesUncertaintys': species_uncertainties,
                   'timeShiftUncertainty': time_shift_uncertainty,
                   'concentrationAbsoluteUncertainty': concentration_absolute_uncertainty,
                   'concentrationRelativeUncertainity': concentration_relative_uncertainity,
                   'moleFractionAbsoluteUncertainty': mole_fraction_absolute_uncertainty,
                   'moleFractionRelativeUncertainty': mole_fraction_relative_uncertainty,
                   'absorbanceAbsoluteUncertainty': absorbance_absolute_uncertainty,
                   'absorbanceRelativeUncertainty': absorbance_relative_uncertainty,
                   'uncertaintyParameterOnes':uncertainty_parameter_ones,
                   'uncertaintyParameterTwos':uncertainty_parameter_twos,
                   'coupledCoefficients':coupledCoefficients,
                   'simulationType':  simulation_type,
                   'parameterOnes':parameter_ones,
                   'parameterTwos':parameter_twos,
                   'functionalForm':functional_form
                   }
            
    def load_yaml_list(self, yaml_list:list = []):
        list_of_yaml_objects = []
        for tup in yaml_list:
            temp = []
            for file in tup:
                temp.append(self.load_to_obj(file))
            list_of_yaml_objects.append(temp) 
        list_of_yaml_objects = [tuple(lst) for lst in list_of_yaml_objects ]               
        return list_of_yaml_objects
    
    def parsing_multiple_dictonaries(self,list_of_yaml_objects:list = [],loop_counter=0):
        experiment_dictonaries = []
        for tup in list_of_yaml_objects:
            if len(tup)>1:
                experiment_dictonaries.append(self.parse_shock_tube_obj(loaded_exp = tup[0],
                                                                        loaded_absorption = tup[1]))

            else:
                experiment_dictonaries.append(self.parse_shock_tube_obj(loaded_exp = tup[0]))
        if loop_counter == 0   :     
            self.original_experimental_conditions = experiment_dictonaries
            
        return experiment_dictonaries
    
    def assemble_dicts_for_master_equation(self,experiment_dictonaries:list=[],
                                           master_equation_reactions:list=[],
                                           additional_parameters:dict={}):
        temperatures = []
        pressures = []
        conditions = []
        master_equation_parameters = []
        for exp in experiment_dictonaries:
            temperatures.append(exp['temperature'])
            pressures.append(exp['pressure'])
            conditions.append(exp['conditions'])
            
        if bool(additional_parameters) == False:
            parameters = {'W':['Energy','Frequencies','SymmetryFactor'],
                          'B':['ImaginaryFrequency']}
            
        for reaction in range(len(master_equation_reactions)):
            temp_dict = {}
            for key in parameters.keys():
                for param in parameters[key]:
                    string = str(key+str(reaction)+'_'+param)
                    temp_dict[string] = [temperatures,pressures,conditions]
            master_equation_parameters.append(temp_dict)
        
        return master_equation_parameters
    
    
    def yaml_file_copy(self,fileName):
    
        tempName = fileName[0:(len(fileName)-5)]
        yamlExtention = fileName[(len(fileName)-5):]
        NewName = tempName +'_updated'+yamlExtention
        shutil.copy2(fileName, NewName) 
    
        return NewName
    
    def yaml_file_updates(self,file_name_list,
                          parsed_yaml_list,
                          experiment_dict_list,
                          physical_observables_updates_list,
                          loop_counter=0):
        
        #always pass in the updated file name list except for the first run of the code
        if loop_counter == 0:
            updated_file_name_list = []
        
        for yaml_file in range(len(file_name_list)):
            temp = []
            if loop_counter == 0: 
                new_file_name = self.yaml_file_copy(file_name_list[yaml_file][0])                 
                temp.append(new_file_name)
                updated_file_name_list.append(temp)
                if len(file_name_list[yaml_file])>1:
                    temp.append(file_name_list[yaml_file][1])
                
            else: 
                new_file_name = file_name_list[yaml_file][0]
    
            if experiment_dict_list[0]['simulation'].physicalSens ==1 :
                temp = self.original_experimental_conditions[yaml_file]['temperature']
                press = self.original_experimental_conditions[yaml_file]['pressure']
                mole_fractions = self.original_experimental_conditions[yaml_file]['MoleFractions']
                conditions = self.original_experimental_conditions[yaml_file]['conditions']
                
                print('__________________________________________________________________________')
                print('loop:',loop_counter)
                print(temp)
                print(press)
                print(conditions)
                print('__________________________________________________________________________')
                
                
                updatedTemp = np.exp(physical_observables_updates_list[yaml_file]['T_experiment_'+str(yaml_file)]) * temp
                updatedTemp = round(updatedTemp,9)
                updatedPress = np.exp(physical_observables_updates_list[yaml_file]['P_experiment_'+str(yaml_file)]) * press
                updatedPress = round(updatedPress,9)
                
                

                species_to_loop =  experiment_dict_list[yaml_file]['uncertainty']['species_relative_uncertainty']['species']
                dilluant = ['Ar','AR','ar','HE','He','he','Kr','KR','kr','Xe','XE','xe','NE','Ne','ne']
                updated_mole_fractions = {}
                count = 0
                for specie in species_to_loop:
                    if specie in dilluant:
                        continue
                    updated = np.exp(physical_observables_updates_list[yaml_file]['X_'+str(count)+'_experiment_'+str(yaml_file)])*conditions[specie]
                    updated = round(updated,9)
                    updated_mole_fractions[specie] = updated
                    count+=1

                for specie in species_to_loop:
                    if specie in dilluant:
                        updated_mole_fractions[specie] = conditions[specie]
                
                updated_mole_fraction_list = []
                for specie in species_to_loop:
                    updated_mole_fraction_list.append(updated_mole_fractions[specie])
 # starting to do file updates here 
               
                with open(new_file_name) as f:
                    config2 = yaml.safe_load(f)
                
                config2['common-properties']['pressure']['value']=float(updatedPress)
                config2['common-properties']['temperature']['value']=float(updatedTemp)
                    
                for i,moleFraction in enumerate(updated_mole_fraction_list):
                    config2['common-properties']['composition'][i]['mole-fraction']=float(moleFraction)
                    
                with open(new_file_name,'w') as f:
                    yaml.safe_dump(config2, f,default_flow_style=False)
                    
  # make a list of updated yaml files here that we return from this and then use these names       
        if loop_counter ==0:  
            return updated_file_name_list
        else:
            return file_name_list
    
    def absorption_file_updates(self,file_name_list,
                          parsed_yaml_list,
                          experiment_dict_list,
                          absorption_observables_updates_dict,
                          loop_counter=0):
        
        
        
        
        #change is happening somewhere after here 
        for yaml_file in range(len(file_name_list)):
            if len(file_name_list[yaml_file])<2:
                continue
            
            if loop_counter ==0:
                new_absorption_file_name = self.yaml_file_copy(file_name_list[yaml_file][1])
                file_name_list[yaml_file][1] = new_absorption_file_name
                
                
            else:
                new_absorption_file_name = file_name_list[yaml_file][1]
                
            
            coupledCoefficients = self.original_experimental_conditions[yaml_file]['coupledCoefficients']
            coupledCoefficentsUpdated = copy.deepcopy(coupledCoefficients)
            ########changes somewhere down there
            
            for species in range(len(coupledCoefficients)):
                for wavelength in range(len(coupledCoefficients[species])):
                    lst = list(coupledCoefficients[species][wavelength])
                    tup = tuple(lst)

                    temp = []
                    for i,values in enumerate(absorption_observables_updates_dict[tup]):
                        
                        temp.append(np.exp(values)*lst[i])

                    coupledCoefficentsUpdated[species][wavelength] = tuple(temp)
                    
            combinationOfNewParameters = list(map(list, list(zip(*(list(map(list, list(zip(*x)))) for x in coupledCoefficentsUpdated)))))
            parameterOnesUpdated = combinationOfNewParameters[0]
             

            for x in range(len(parameterOnesUpdated)):
                for y in range(len(parameterOnesUpdated[x])):
                    parameterOnesUpdated[x][y] = round(parameterOnesUpdated[x][y], 8)

        

            parameterTwosUpdated = combinationOfNewParameters[1]
            for x in range(len(parameterTwosUpdated)):
                for y in range(len(parameterTwosUpdated[x])):
                    parameterTwosUpdated[x][y] = round(parameterTwosUpdated[x][y], 8)
            

            with open(new_absorption_file_name)as f:
                config3 = yaml.safe_load(f)
            
            
            
            for parameterOne in range(len(config3['Absorption-coefficients'])):
                for wavelength in range(len(config3['Absorption-coefficients'][parameterOne]['wave-lengths'])):        
                    config3['Absorption-coefficients'][parameterOne]['wave-lengths'][wavelength]['parameter-one']['value'] = float(parameterOnesUpdated[parameterOne][0])

        
     

            for parameterTwo in range(len(config3['Absorption-coefficients'])):
                for wavelength in range(len(config3['Absorption-coefficients'][parameterTwo]['wave-lengths'])):
                    config3['Absorption-coefficients'][parameterTwo]['wave-lengths'][wavelength]['parameter-two']['value'] = float(parameterTwosUpdated[parameterTwo][0])
            
            with open(new_absorption_file_name,'w') as f:
                yaml.safe_dump(config3, f,default_flow_style=False)
        
        return file_name_list
    
    
                
                
        
        
    