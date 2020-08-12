import cantera as ct
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
gas = ct.Solution('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/variable_pressure_shock_tube/jetsurf2.cti')
#gas.TPX = T0, P0, 'IC8H18:0.016503292, N2:0.785871069, CO2:0.197625639'
gas.TPX = 1350, 101325*1.4, 'NC7H16:.002, O2:0.022, Ar:.976'
pressure_trace = pd.read_csv('/Users/carlylagrotta/Dropbox/Columbia/MSI/data/variable_pressure_shock_tube/pressure_trace.csv')
print(pressure_trace.shape)
S0 = gas.entropy_mass
rho0 = gas.density
print(gas.SV)
Pressure_array = pressure_trace['Pressure'].to_numpy()*101325
plt.figure()
plt.plot(pressure_trace['Time'],pressure_trace['Pressure'].to_numpy()*101325)
#V0 = gas.
# for i in range(pressure_trace.shape[0]):
#     P = pressure_trace['Pressure'][i] * 101325
#     gas.SP = S0, P
#     T_new = gas.T
#     V_new = V0 * rho0 / gas.density

    # print / save P,T_new,V_new


class VolumeFromPressure(object):
    r"""Create a volume trace given a pressure trace.
    Using Cantera to evaluate the thermodynamic properties, compute a
    volume trace from a pressure trace.
    Parameters
    ----------
    pressure : `numpy.ndarray`
        1-D array containing the reactor pressure
    v_initial : `float`
        Initial volume of the experiment, in m**3
    T_initial : `float`, optional
        Initial temperature of the experiment, in Kelvin. Optional for
        Cantera versions greater than 2.2.0.
    chem_file : `str`, optional
        Filename of the chemistry file to be used
    Attributes
    ----------
    volume : `numpy.ndarray`
        The volume trace
    Notes
    -----
    The volume is computed according to the formula
    .. math:: v_i = v_{initial}*\rho_{initial}/\rho_i
    where the index :math:`i` indicates the current point. The state
    is set at each point by setting the pressure from the input array
    and the entropy to be constant. The volume is computed by the
    isentropic relationship described above.
    """
    def __init__(self, pressure, v_initial, T_initial, chem_file='species.cti', cti_source=None):
        if cti_source is None:
            gas = ct.Solution(chem_file)
        else:
            gas = ct.Solution(source=cti_source)
        gas.TP = T_initial, pressure[0]
        initial_entropy = gas.entropy_mass
        initial_density = gas.density
        self.volume = np.zeros((len(pressure)))
        for i, p in enumerate(pressure):
            gas.SP = initial_entropy, p
            self.volume[i] = v_initial*initial_density/gas.density

    def __repr__(self):
        return 'VolumeFromPressure(volume={self.volume!r})'.format(self=self)
instance = VolumeFromPressure(Pressure_array,.19,1350,chem_file='/Users/carlylagrotta/Dropbox/Columbia/MSI/data/variable_pressure_shock_tube/jetsurf2.cti')
