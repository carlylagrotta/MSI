import MSI.simulations as sim
import re
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import MSI.simulations.instruments.shock_tube as st


#acts as front end to the rest of the system


# takes the data from one experiment and puts it in a dict of dicts
# that follows the format of the S matrix
# only need the last 3 elements of the interpolated absorbance
# absorbance is of form [interp_original,interp_abs_kinetic_sens,interp_abs_phys_sens,interp_abs_coef_sens]
# where each list element is a dict. keys are wavelengths, values are the sensitivities for that wavelength 
# psens should match interpolated_tp and species sens in size, but again is dict with wavelength keys
# index from 1, so if you have 3 experiments, their indices will be 1,2,3
class Optimization_Utility(object):
    def __init__(self):
        self.matrix = None
        
        
    def build_single_exp_dict(self,exp_index:int,
                              simulation:sim.instruments.shock_tube.shockTube,
                              interpolated_kinetic_sens:dict,
                              interpolated_tp_sens:list,
                              interpolated_species_sens:list,
                              interpolated_absorbance:list=[],
                              experimental_data:list =[],
                              absorbance_experimental_data:list=[],
                              time_history_interpolated_against_absorbance_experiment:dict={},
                              absorbance_calculated_from_model=None):
        exp_dict = {}
        exp_dict['index']              = exp_index
        exp_dict['simulation']         = simulation
       
        if interpolated_kinetic_sens==None:
            exp_dict['ksens']  = None
            exp_dict['temperature'] = None
            exp_dict['pressure'] = None
            exp_dict['species'] = None
        else:
            exp_dict['ksens']              = interpolated_kinetic_sens
            exp_dict['temperature']        = interpolated_tp_sens[0]
            exp_dict['pressure']           = interpolated_tp_sens[1]
            exp_dict['species']            = interpolated_species_sens
            
        exp_dict['observables']        = simulation.observables
        exp_dict['concentration_observables'] = simulation.concentrationObservables
        exp_dict['mole_fraction_observables'] = simulation.moleFractionObservables
        #needs to be in the order of mole fraction csv files + concentration csv files 
        exp_dict['experimental_data']  = experimental_data
        # start here 
        exp_dict['uncertainty']        = self.build_uncertainty_shock_tube_dict(exp_dict['simulation'].fullParsedYamlFile)
        #decide how we want to build uncertainty dict and if we want to pass in the parsed yaml file?
        
        if len(interpolated_absorbance) != 0:
            exp_dict['absorbance_model_data'] = interpolated_absorbance[0]
            exp_dict['absorbance_ksens']   = interpolated_absorbance[1]
            exp_dict['absorbance_psens']   = interpolated_absorbance[2]
            exp_dict['perturbed_coef']     = interpolated_absorbance[3]
            exp_dict['absorbance_observables'] = simulation.absorbanceObservables
            exp_dict['absorbance_experimental_data'] = absorbance_experimental_data
            exp_dict['absorbance_calculated_from_model'] = absorbance_calculated_from_model
            exp_dict['time_history_interpolated_against_abs'] = time_history_interpolated_against_absorbance_experiment
            
            
        return exp_dict
    
    
    def load_exp_from_file(self,yaml_exp_file_list = []):
        for file in yaml_exp_file_list:
            continue
        
    def build_uncertainty_shock_tube_dict(self,experiment_dictonarie:dict={}):
        uncertainty_dict = {}
        #need to make an exception to this if there is no absortpion in dict
        if 'coupledCoefficients' in experiment_dictonarie.keys():
            coupled_coefficients = experiment_dictonarie['coupledCoefficients']
            coupled_coefficients = [item for sublist in coupled_coefficients for item in sublist]
            uncertain_parameters_ones = experiment_dictonarie['uncertaintyParameterOnes']
            uncertain_parameter_twos = experiment_dictonarie['uncertaintyParameterTwos']
            zip_uncertain_paramters = list(zip(uncertain_parameters_ones,uncertain_parameter_twos))
            dict_of_coupled_unc_and_param = dict(zip(coupled_coefficients,zip_uncertain_paramters))
            
            
            uncertainty_dict['coupled_coef_and_uncertainty'] = dict_of_coupled_unc_and_param
            uncertainty_dict['absorbance_relative_uncertainty'] = experiment_dictonarie['absorbanceRelativeUncertainty']
            uncertainty_dict['absorbance_absolute_uncertainty'] = experiment_dictonarie['absorbanceAbsoluteUncertainty']
        #finish making this dictonary
        uncertainty_dict['temperature_relative_uncertainty'] = experiment_dictonarie['tempRelativeUncertainty']
        uncertainty_dict['pressure_relative_uncertainty'] = experiment_dictonarie['pressureRelativeUncertainty']
        uncertainty_dict['species_relative_uncertainty'] = {'dictonary_of_values':experiment_dictonarie['speciesUncertaintys'],
                        'species':experiment_dictonarie['speciesNames']}
        uncertainty_dict['mole_fraction_relative_uncertainty'] = experiment_dictonarie['moleFractionRelativeUncertainty']
        uncertainty_dict['mole_fraction_absolute_uncertainty'] = experiment_dictonarie['moleFractionAbsoluteUncertainty'] 
        uncertainty_dict['concentration_relative_uncertainty'] = experiment_dictonarie['concentrationRelativeUncertainity']
        uncertainty_dict['concentration_absolute_uncertainty'] = experiment_dictonarie['concentrationAbsoluteUncertainty']
        
        return uncertainty_dict
    
    
    
    
    def running_full_shock_tube(self,processor=None,
                                           experiment_dictonary:dict={},
                                           kineticSens = 1,
                                           physicalSens = 1,
                                           dk = .01,
                                           exp_number=1):
        shock_tube = st.shockTube(pressure = experiment_dictonary['pressure'],
                     temperature = experiment_dictonary['temperature'],
                     observables = experiment_dictonary['observables'],
                     kineticSens = kineticSens,
                     physicalSens = physicalSens,
                     conditions = experiment_dictonary['conditions'],
                     initialTime = experiment_dictonary['initialTime'],
                     finalTime = experiment_dictonary['finalTime'],
                     thermalBoundary = experiment_dictonary['thermalBoundary'],
                     mechanicalBoundary = experiment_dictonary['mechanicalBoundary'],
                     processor = processor,
                     save_timeHistories = 1,
                     save_physSensHistories = 1,
                     moleFractionObservables = experiment_dictonary['moleFractionObservables'],
                     concentrationObservables = experiment_dictonary['concentrationObservables'],
                     fullParsedYamlFile = experiment_dictonary)
        
        csv_paths = [x for x in  experiment_dictonary['moleFractionCsvFiles'] + experiment_dictonary['concentrationCsvFiles'] if x is not None]
        exp_data = shock_tube.importExperimentalData(csv_paths)
        
        shock_tube.run()
        ########################################################################

        #check this tomorrow 
        ################################################################################
        int_ksens_exp_mapped= shock_tube.map_and_interp_ksens()#ksens is wiped on rerun so int it before
        shock_tube.sensitivity_adjustment(temp_del = dk)
        shock_tube.sensitivity_adjustment(pres_del = dk)
        shock_tube.species_adjustment(dk)
        ############################################### check to make sure these aren't effected 
        int_tp_psen_against_experimental = shock_tube.interpolate_experimental([shock_tube.interpolate_physical_sensitivities(index=1),
                                                                           shock_tube.interpolate_physical_sensitivities(index=2)])
    
        int_spec_psen_against_experimental = shock_tube.interpolate_experimental(pre_interpolated=shock_tube.interpolate_species_sensitivities())
    ###############saving the shock tube experimental interpolated time history     
        single_data = shock_tube.interpolate_experimental(single=shock_tube.timeHistories[0])
        shock_tube.savingInterpTimeHistoryAgainstExp(single_data)
        #tab starting here tomorrow
        shock_tube.interpolatePressureandTempToExperiment(shock_tube,exp_data)
    ###############  ###############  
        experiment = self.build_single_exp_dict(exp_number,
                                           shock_tube,
                                           int_ksens_exp_mapped,
                                           int_tp_psen_against_experimental,
                                           int_spec_psen_against_experimental,
                                           experimental_data = exp_data)
        
        #write test case and check if we can get as far as just returnign the experiment
        return experiment
    
    def running_full_shock_tube_absorption(self,processor=None,
                                           experiment_dictonary:dict={},
                                           absorbance_yaml_file_path = '',
                                           kineticSens = 1,
                                           physicalSens = 1,
                                           dk = .01,
                                           exp_number=1):
        shock_tube = st.shockTube(pressure = experiment_dictonary['pressure'],
                     temperature = experiment_dictonary['temperature'],
                     observables = experiment_dictonary['observables'],
                     kineticSens = kineticSens,
                     physicalSens = physicalSens,
                     conditions = experiment_dictonary['conditions'],
                     initialTime = experiment_dictonary['initialTime'],
                     finalTime = experiment_dictonary['finalTime'],
                     thermalBoundary = experiment_dictonary['thermalBoundary'],
                     mechanicalBoundary = experiment_dictonary['mechanicalBoundary'],
                     processor = processor,
                     save_timeHistories = 1,
                     save_physSensHistories = 1,
                     moleFractionObservables = experiment_dictonary['moleFractionObservables'],
                     absorbanceObservables = experiment_dictonary['absorbanceObservables'],
                     concentrationObservables = experiment_dictonary['concentrationObservables'],
                     fullParsedYamlFile = experiment_dictonary)
    
        
        csv_paths = [x for x in  experiment_dictonary['moleFractionCsvFiles'] + experiment_dictonary['concentrationCsvFiles'] if x is not None]
        
        exp_data = shock_tube.importExperimentalData(csv_paths)
        shock_tube.run()
        #this might be in the wrong spot 
        int_ksens_exp_mapped= shock_tube.map_and_interp_ksens()
    
    
        abs_instance = csp.Absorb()
        parser = yp.Parser()
        abs_loaded = parser.load_to_obj(absorbance_yaml_file_path)
        abs_data = abs_instance.superimpose_shock_tube(shock_tube,abs_loaded,experiment_dictonary['pathLength'],
                                                       kinetic_sens=kineticSens)    
        
        
        perturbed_coef = abs_instance.perturb_abs_coef(dk,
                                              shock_tube,
                                              abs_loaded,
                                              experiment_dictonary['pathLength'],
                                              summed_data = abs_data[0]) 
        
        
        
        shock_tube.sensitivity_adjustment(temp_del = dk)
        shock_tube.sensitivity_adjustment(pres_del = dk)
        shock_tube.species_adjustment(dk)
        int_tp_psen_against_experimental = shock_tube.interpolate_experimental([shock_tube.interpolate_physical_sensitivities(index=1),
                                                                                 shock_tube.interpolate_physical_sensitivities(index=2)])
        
        int_spec_psen_against_experimental = shock_tube.interpolate_experimental(pre_interpolated=shock_tube.interpolate_species_sensitivities())
        
        abs_phys_sens = abs_instance.absorb_phys_sensitivities(shock_tube,abs_data[0],abs_loaded,
                                                               experiment_dictonary['pathLength'],
                                                               dk = dk)
       


        loaded_experimental_data_absorbance = abs_instance.import_experimental_data(experiment_dictonary['absorbanceCsvFiles'])
        
        interp_abs_exp= abs_instance.interpolate_experimental(shock_tube,loaded_experimental_data_absorbance,
                                                            original_summed_absorption=abs_data[0],
                                                            abs_kinetic_sens = abs_data[1],
                                                            abs_phys_sens = abs_phys_sens,
                                                            abs_coef_sens = perturbed_coef)

        time_history_interp_against_experiment_dict = abs_instance.interpolate_experimental(shock_tube,
                                                                                            loaded_experimental_data_absorbance,
                                                                                            time_history = shock_tube.timeHistories[0])
     #################################################################################   
        single_data = shock_tube.interpolate_experimental(single=shock_tube.timeHistories[0])
        shock_tube.savingInterpTimeHistoryAgainstExp(single_data)
        shock_tube.interpolatePressureandTempToExperiment(shock_tube,exp_data)
    ####################################################################################    
        
        experiment = self.build_single_exp_dict(exp_number,
                                                shock_tube,
                                      int_ksens_exp_mapped,
                                      int_tp_psen_against_experimental,
                                      int_spec_psen_against_experimental,
                                      interpolated_absorbance=interp_abs_exp,
                                      experimental_data = exp_data,
                                      absorbance_experimental_data = loaded_experimental_data_absorbance,
                                      time_history_interpolated_against_absorbance_experiment = time_history_interp_against_experiment_dict,
                                      absorbance_calculated_from_model = abs_data[0])
        
        return experiment
    
    
    
    def running_shock_tube_absorption_only(self,processor=None,
                                           experiment_dictonary:dict={},
                                           absorbance_yaml_file_path = '',
                                           kineticSens = 1,
                                           physicalSens = 1,
                                           dk = .01,
                                           exp_number=1):
        shock_tube = st.shockTube(pressure = experiment_dictonary['pressure'],
                     temperature = experiment_dictonary['temperature'],
                     observables = experiment_dictonary['observables'],
                     kineticSens = kineticSens,
                     physicalSens = physicalSens,
                     conditions = experiment_dictonary['conditions'],
                     initialTime = experiment_dictonary['initialTime'],
                     finalTime = experiment_dictonary['finalTime'],
                     thermalBoundary = experiment_dictonary['thermalBoundary'],
                     mechanicalBoundary = experiment_dictonary['mechanicalBoundary'],
                     processor = processor,
                     save_timeHistories = 1,
                     save_physSensHistories = 1,
                     moleFractionObservables = experiment_dictonary['moleFractionObservables'],
                     absorbanceObservables = experiment_dictonary['absorbanceObservables'],
                     concentrationObservables = experiment_dictonary['concentrationObservables'],
                     fullParsedYamlFile = experiment_dictonary)

    
        shock_tube.run()
        abs_instance = csp.Absorb()
        parser = yp.Parser()
        abs_loaded = parser.load_to_obj(absorbance_yaml_file_path)
        abs_data = abs_instance.superimpose_shock_tube(shock_tube,abs_loaded,experiment_dictonary['pathLength'],
                                                       kinetic_sens=kineticSens)
        
        #print('first go')
        
        
        perturbed_coef = abs_instance.perturb_abs_coef(dk,
                                              shock_tube,
                                              abs_loaded,
                                              experiment_dictonary['pathLength'],
                                              summed_data = abs_data[0]) 
        
        #print('second go')
       
        #print(perturbed_coef)
        
        
        shock_tube.sensitivity_adjustment(temp_del = dk)
        shock_tube.sensitivity_adjustment(pres_del = dk)
        shock_tube.species_adjustment(dk)        
        abs_phys_sens = abs_instance.absorb_phys_sensitivities(shock_tube,abs_data[0],abs_loaded,
                                                               experiment_dictonary['pathLength'],
                                                               dk = dk)
       

        
        loaded_experimental_data_absorbance = abs_instance.import_experimental_data(experiment_dictonary['absorbanceCsvFiles'])
        
        interp_abs_exp= abs_instance.interpolate_experimental(shock_tube,loaded_experimental_data_absorbance,
                                                            original_summed_absorption=abs_data[0],
                                                            abs_kinetic_sens = abs_data[1],
                                                            abs_phys_sens = abs_phys_sens,
                                                            abs_coef_sens = perturbed_coef)
        
        
        time_history_interp_against_experiment_dict = abs_instance.interpolate_experimental(shock_tube,
                                                                                            loaded_experimental_data_absorbance,
                                                                                            time_history = shock_tube.timeHistories[0])
        experiment = self.build_single_exp_dict(exp_number,
                                                shock_tube,
                                      None,
                                      None,
                                      None,
                                      interpolated_absorbance=interp_abs_exp,
                                      absorbance_experimental_data = loaded_experimental_data_absorbance,
                                      time_history_interpolated_against_absorbance_experiment = time_history_interp_against_experiment_dict,
                                      absorbance_calculated_from_model = abs_data[0] )
        
        return experiment    
    
    
    
    
    
    def looping_over_parsed_yaml_files(self,list_of_parsed_yamls,list_of_yaml_paths,processor=None,kineticSens=1,physicalSens=1,dk=.01):
        experiment_list = []
        for i,yamlDict in enumerate(list_of_parsed_yamls):
         
            simulation_type = yamlDict['simulationType']
            
            if re.match('[Ss]hock [Tt]ube',simulation_type):
                simulation_type = 'shock tube'

                if simulation_type == 'shock tube':
                    if 'absorbanceObservables' not in yamlDict.keys():
                        experiment = self.running_full_shock_tube(processor=processor,
                                           experiment_dictonary=yamlDict,
                                           kineticSens = kineticSens,
                                           physicalSens = physicalSens,
                                           dk = dk,
                                           exp_number=i)
                        experiment_list.append(experiment)
                        ####FINISH writing this function and start writing main function tomorrow 
                    elif 'absorbanceObservables' in yamlDict.keys() and yamlDict['moleFractionObservables'][0] == None and yamlDict['concentrationObservables'][0]==None:
                        path = list_of_yaml_paths[i][1]
                        print(path)
                        experiment = self.running_shock_tube_absorption_only(processor=processor,
                                                                             experiment_dictonary = yamlDict,
                                                                             absorbance_yaml_file_path = path,
                                                                             kineticSens = kineticSens,
                                                                             physicalSens = physicalSens,
                                                                             dk = dk,
                                                                             exp_number=i)
                        experiment_list.append(experiment)
                    
                    else:
                        path = list_of_yaml_paths[i][1]
                        experiment = self.running_full_shock_tube_absorption(processor=processor,
                                           experiment_dictonary=yamlDict,
                                           absorbance_yaml_file_path = path,
                                           kineticSens = kineticSens,
                                           physicalSens = physicalSens,
                                           dk = dk,
                                           exp_number=i)
                        experiment_list.append(experiment)
                                            
            else:
                print('We do not have this simulation installed yet')
            
        return experiment_list
    def saving_experimental_dict(self,list_off_experiment_dictonaires):
        uncertainty_list = []
        for i,exp in enumerate(list_off_experiment_dictonaires):
            if 'perturbed_coef' not in exp.keys():
                uncertainty_list.append({})
            else:
                uncertainty_list.append(exp['uncertainty'])
                            
        return uncertainty_list
    


   
    

    
    
