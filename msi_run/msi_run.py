# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:33:39 2020

@author: Skoron
"""
import sys
import fire


class multiscale_informatics:
    
    
    def __init__(self,input_options):
        self.input_options=input_options
        self.get_parameters()
        
    def get_parameters(self):
        self.get_working_directory()
        self.get_iterations()
        self.get_yaml_files()
        self.get_kinetic_models()
        self.get_reaction_uncertainties_list()
        self.get_master_equation_models()
        self.get_master_equation_uncertainties()
        self.get_targets()
        self.get_optional_plotting_targets()
        self.get_master_equation_sens()
        self.get_master_reactions()
        
    def get_working_directory(self):
        self.wdir=self.input_options[0].split('=')[1]
    
    def get_iterations(self):
        self.iterations=int(self.input_options[1].split('=')[1])
        
    def get_yaml_files(self):
        yaml_bool=False
        self.yaml_files=[]
        for i,string in enumerate(self.input_options):
            if 'end_yaml' in string:
                break
            if yaml_bool:
                temp_entry=string.lstrip('[')
                temp_entry=temp_entry.rstrip(']')
                entry=temp_entry.split(',')
                #print(entry)
                if entry[1]=='':
                    entry=[entry[0]]
                self.yaml_files.append(entry)               
                
                
            if 'begin_yaml_list' in string:
                yaml_bool=True
            
            
    def get_kinetic_models(self):
        model_bool=False
        self.models=[]
        for i,string in enumerate(self.input_options):
            if 'end_model_list' in string:
                break
            if model_bool:
                self.models.append(string)
            if 'begin_model_list' in string:
                model_bool=True
            
            
    def get_reaction_uncertainties_list(self):
        uncertainties_bool=False
        self.reaction_uncertainties=[]
        for i,string in enumerate(self.input_options):
            if 'end_reaction_uncertainty_list' in string:
                break
            if uncertainties_bool:
                self.reaction_uncertainties.append(string)
            if 'begin_reaction_uncertainty_list' in string:
                uncertainties_bool=True
    def get_master_equation_models(self):
        master_bool=False
        self.master_equation_models=[]
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_model_list' in string:
                break
            if master_bool:
                self.master_equation_models.append(string)
            if 'begin_master_equation_model_list' in string:
                master_bool=True
                
                
    def get_master_equation_uncertainties(self):
        master_bool=False
        self.master_uncertainties=[]
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_uncertainties' in string:
                break
            if master_bool:
                self.master_uncertainties.append(string)
            if 'begin_master_equation_uncertainties' in string:
                master_bool=True
                
    def get_targets(self):
        targets_bool=False
        self.targets=[]
        for i,string in enumerate(self.input_options):  
            if 'end_rate_constant_targets' in string:
                break
            if targets_bool:
                self.targets.append(string)
            if 'begin_rate_constant_targets' in string:
                targets_bool=True
        
    def get_optional_plotting_targets(self):
        targets_bool=False
        self.optional_targets=[]
        for i,string in enumerate(self.input_options):  
            if 'end_optional_plotting_targets' in string:
                break
            if targets_bool:
                self.optional_targets.append(string)
            if 'begin_optional_plotting_targets' in string:
                targets_bool=True
    
    def get_master_equation_sens(self):
        sens_bool=False
        self.master_sens=[]
        for i,string in enumerate(self.input_options):  
            if 'end_master_equation_sensitivities' in string:
                break
            if sens_bool:
                self.master_sens.append(string)
            if 'begin_master_equation_sensitivities' in string:
                sens_bool=True
            
            
    def get_master_reactions(self):
        rxn_list_len=len(self.models)
        set_names=[]
        reactions_bool=False
        set_bools=[False]*rxn_list_len
        self.master_equation_reactions=[]
        for i in range(rxn_list_len):
            set_names.append('set'+str(i+1))
        for i,string in enumerate(self.input_options):
            if 'end_master_equation_reactions' in string:
                break
            if reactions_bool:
                for j,string1 in enumerate(set_bools):
                    #print(j)
                    index_top=self.input_options.index('set'+str(j+1))
                    index_bottom=self.input_options.index('end set'+str(j+1))
                    self.master_equation_reactions.append(self.input_options[index_top+1:index_bottom])
                break
                
                self.master_sens.append(string)
            if 'begin_master_equation_reactions' in string:
                reactions_bool=True
        


def parser(input_file):
    #print(#input_file)
    
    
    with open(input_file,'r') as f:
        input_lines=f.readlines()
        
    for i in range(len(input_lines)):
        #rint('char='+str(input_lines[i]))
        input_lines[i]=input_lines[i].lstrip()
        input_lines[i]=input_lines[i].rstrip('\n')
        input_lines[i]=input_lines[i].rstrip()
        if input_lines[i]:
            if input_lines[i][0]=='#':
               input_lines[i]=''
        
        
    while '' in input_lines:
        #print('a')
        input_lines.remove('')
    
      
    return input_lines
    
        



def main(input_file=''):

    if input_file=='':
        print('Please run program with defined input file using --input_file=FILEPATH')

    elif input_file !='':
        input_options=parser(input_file)
        simulation=multiscale_informatics(input_options)
        
    return simulation

if __name__ == '__main__':
    a=fire.Fire(main)