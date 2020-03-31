#Sigma Adjuster
import pandas as pd
import numpy as np
import yaml
def sigma_adjuster(sigma):
    file_range = 28
    for file in range(file_range+1):
        df = pd.read_csv('masten_oh_'+str(file)+'.csv')
        shape = df.shape[0]
        new_sigma = np.array([sigma]*shape)
        new_sigma = new_sigma.reshape((new_sigma.shape[0],1))
        new_sigma = pd.DataFrame(new_sigma,columns=['Relative_Uncertainty'])
        df['Relative_Uncertainty'] = new_sigma['Relative_Uncertainty']
        df.to_csv('masten_oh_'+str(file)+'.csv', index = False)

sigma_adjuster(.25)
#sigma_adjuster(1)
#sigma_adjuster(.1)

def sigma_adjuster_2(sigma):
    file_range = 10
    for file in range(file_range+1):
        df = pd.read_csv('masten_oh_'+str(file)+'.csv')
        shape = df.shape[0]
        new_sigma = np.array([sigma]*shape)
        new_sigma = new_sigma.reshape((new_sigma.shape[0],1))
        new_sigma = pd.DataFrame(new_sigma,columns=['Relative_Uncertainty'])
        df['Relative_Uncertainty'] = new_sigma['Relative_Uncertainty']
        df.to_csv('masten_oh_'+str(file)+'.csv', index = False)

#sigma_adjuster(.15)
        
        
def temp_uncertainty_adjust(sigma):
    def load_to_obj(path:str = ''):
        with open(path) as f:
            config = yaml.load(f,Loader=yaml.FullLoader)
        return config
    
    file_range = 29
    for file in range(file_range):
        config = load_to_obj('Masten_'+str(file)+'.yaml')
        config['common-properties']['temperature']['relative-uncertainty'] = sigma
        
        with open('Masten_'+str(file)+'.yaml','w') as f:        
            yaml.safe_dump(config, f,default_flow_style=False)
          #was .01  
#temp_uncertainty_adjust(.01)
#temp_uncertainty_adjust(.03)
