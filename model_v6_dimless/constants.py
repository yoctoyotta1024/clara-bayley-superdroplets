### file for values of global constants 
### and useful helper functions for
### quick small calculations

import numpy as np
import matplotlib.pyplot as plt
#from typhon.physics import e_eq_water_mk as saturation_pressure

plt.rcParams.update({'font.size': 12})




######################## constants of problem ########################
### NOTE: Any variation in these quantities are ignored (e.g. temperature depedence)
G            = 9.80665           # acceleration of gravity [m/s^2]
R            = 8.314462618       # Universal Molar Gas Constant [J/mol/K]
Mr_water     = 0.01801528        # molecular mass of water
Mr_dry       = 0.0289647         # molecular mass (per 1 mole) of dry air (M_r)       
mr_ratio     = Mr_water/Mr_dry


Rgas_dry     = R/Mr_dry          # specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
Rgas_v       = R/Mr_water        # specific gas constant for water
Latent_v     = 2437300           # specific latent heat of vapourisation of water [J/Kg]  (at 300K)
Dynvisc      = 18.45*1e-6        # dynamic viscosity of air [Pa s] (at 300K)

Cp_dry       = 1004.9            # specific heat capacity (dry) air at const. pressure [J/Kg/K] (taken at 300K)   = 1.400*cv_dry 
Cp_v         = 1864              # specific heat capacity of water vapour [J/Kg/K] (at 300K)
C_l          = 4180              # specific heat capacity of liquid water[J/Kg/K] (at 300K)

Rho_dry      = 1.177             # density of dry air [Kg/m^3] (at 300K)
Rho_l        = 996.57            # density of liquid water condensing [kg/m^3] (at 300K)

Rho_sol        = 2200                             # density of dry areosol [Kg m^-3]
Mr_sol         = 0.058443                         # Mr of solute (dry areosol) [Kg/mol]
ionic          = 2                                # degree ionic dissociation (van't Hoff factor)


######################## Dimensional Parameters of Problem #################
W0  = 0.5           # [m/s]
T0 = 4000           # [s]

P0 = 100000         # [Pa]
Temp0 = 273.14      # [K]
Rho0 = P0/(Rgas_dry*Temp0)        # [Kg]

Cp0 = Cp_dry

R0 = 1e-6                      # [m]
N0 = 1e6                       # [m^-3]
Rho0 = P0/(Rgas_dry*Temp0)     # [Kg/m^3]
Mr0 = Mr_dry                   # [Kg/mol]
############################################################################





######################## de-dimensionalise problem ########################
cp_dry = Cp_dry / Cp0
cp_v = Cp_v / Cp0
c_l = C_l / Cp0
latent_v = Latent_v/(Temp0*Cp0)
rgas_dry = Rgas_dry/Cp0

rho_dry = Rho_dry/Rho0
rho_l = Rho_l/Rho0

rho_sol = Rho_sol/Rho0
mr_sol = Mr_sol/Mr0 
############################################################################



def dry_pot_temp(temp, p, qv):
    ''' calculate potential temperature
    assuming moist (unsaturated) air with
    vapour content qv '''

    Cp = Cp_dry * (1+qv*Cp_v/Cp_dry)/(1+qv)
    Rgas = Rgas_dry *(1+qv*Rgas_v/Rgas_dry)/(1+qv)
    
    theta = temp*(p[0]/p)**(Rgas/Cp)
    
    return theta






def moist_static_energy(z, temp, qv=0):
    ''' calculate the moist static energy /m^3
        (not assuming dry air cp) '''
    
    return  G*z + Latent_v*qv + Cp_dry*temp
