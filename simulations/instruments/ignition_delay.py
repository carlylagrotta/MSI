import cantera as ct
from .. import simulation as sim
from ...cti_core import cti_processor as ctp
import pandas as pd
import numpy as np
import time
import copy
import re




class ignition_delay(sim.simulation):
    
    
    def __init__(self,pressure:float,temperature:float,observables:list,target:str,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={}, save_timeHistories:int=0,
                 log_file=True,log_name='log.txt'):
        
        
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
        self.log_name=log_name
        self.log_file=log_file
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
        