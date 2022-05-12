### file for values of global constants 
### and useful helper functions for
### quick small calculations

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import g, R
#from typhon.physics import e_eq_water_mk as saturation_pressure

plt.rcParams.update({'font.size': 12})


### some constants
mr_water     = 0.01801528        # molecular mass of water
mr_dry       = 0.0289647         # molecular mass (per 1 mole) of dry air (M_r)                 
mr_ratio     = mr_water/mr_dry
rgas_dry     = R/mr_dry          # specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
rgas_v       = R/mr_water        # specific gas constant for water
cp_dry       = 1004.9            # specific heat capacity (dry) air at const. pressure [J/Kg/K] (taken at 300K)   = 1.400*cv_dry 
cp_v         = 1864              # specific heat capacity of water vapour [J/Kg/K] (at 300K)
c_l          = 4180              # specific heat capacity of liquid water[J/Kg/K] (at 300K)
rho_dry      = 1.177             # density of dry air [Kg/m^3] (at 300K)
latent_v     = 2437300           # specific latent heat of vapourisation of water [J/Kg]  (at 300K)
rho_l        = 996.57            # density of liquid water condensing [kg/m^3] (at 300K)
dynvisc      = 18.45*1e-6        # dynamic viscosity of air [Pa s] (at 300K)
### NOTE! Any variation in above quantities are ignored (e.g. temperature depedence)





def saturation_pressure(T):
    ''' Calculate the equilibrium vapor pressure
    of water over liquid water. Equation taken from
    typhon.physics.thermodynamics.e_eq_water_mk '''

    if np.any(T <= 0):
        err = 'T must be larger than 0K.'
        raise ValueError(err)

    lnpsat = (54.842763      # ln(saturation vapor pressure) [Pa]
         - 6763.22 / T
         - 4.21 * np.log(T)
         + 0.000367 * T
         + np.tanh(0.0415 * (T - 218.8))
         * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025 * T))

    return np.exp(lnpsat)







def pv2qv(pv, p):
    ''' calculate mass mixing ratio
    qv = m_v/m_dry =rho_v/rho_dry
    given vapour pressure pv = p_v/p_tot '''
    
    qv = mr_ratio * pv/(p-pv)
    
    return qv



def qv2pv(qv, p):
    ''' calculate vapour pressure
    pv = p_v/p_tot
    given mass mixing ratio 
    qv = m_v/m_dry = rho_v/rho_dry'''
    
    pv = p*qv/(mr_ratio+qv)
    
    return pv


def cp_moist(qv, qc):
    ''' calculate effecitve specific heat 
    capacity of parcel with water vapour 
    qv=m_v/m_dry and liquid water qc= m_l/m_dry
    (notice m_dry not m_tot). Enthalpy change of
    parcel per unit (total) mass = dh. 
    (1+qv+qc)dh = cp_eff*dtemp + latent_v*dqv'''

    return cp_dry + cp_v*qv + c_l*qc 
 

def dry_pot_temp(temp, p, qv):
    ''' calculate potential temperature
    assuming moist (unsaturated) air with
    vapour content qv '''
    
    cp = cp_dry * (1+qv*cp_v/cp_dry)/(1+qv)
    rgas = rgas_dry *(1+qv*rgas_v/rgas_dry)/(1+qv)
    
    theta = temp*(p[0]/p)**(rgas/cp)
    
    return theta
