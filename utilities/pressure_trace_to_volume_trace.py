import cantera as ct
gas = ct.Solution('gri30.cti')
#gas.TPX = T0, P0, 'IC8H18:0.016503292, N2:0.785871069, CO2:0.197625639'
gas.TPX = 1000, 101325, 'IC8H18:0.016503292, N2:0.785871069, CO2:0.197625639'

S0 = gas.entropy_mass
rho0 = gas.density

for line in pressure_trace:
    P = float(line[1])
    gas.SP = S0, P
    T_new = gas.T
    V_new = V0 * rho0 / gas.density

    # print / save P,T_new,V_new
