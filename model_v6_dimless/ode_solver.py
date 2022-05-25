### FUNCTIONS FOR ODE SOLVER
import numpy as np
import matplotlib.pyplot as plt

from constants import *
from superdroplets_as_objects import *


######################## constants of problem ########################
### NOTE: Any variation in these quantities are ignored (e.g. temperature depedence)
Thermk273 = 0.023822892 - 7.11756e-5*273.15              #[eq.7.24 lohmann intro 2 clouds] [J/s/m/K]
Diffuse_v273 = 2.11e-5 * 101325 / (273.15**1.94)         #[eq.7.26 lohmann intro 2 clouds] [m^2/s]



######################## small calculations involved ########################
def pv2qv(pv, p):
    ''' calculate mass mixing ratio
    qv = m_v/m_dry = rho_v/rho_dry
    given vapour pressure pv = p_v/p_tot.'''
    
    qv = mr_ratio * pv/(p-pv)
    
    return qv


def qv2pv(qv, p):
    ''' calculate vapour pressure
    pv = p_v/p_tot
    given mass mixing ratio 
    qv = m_v/m_dry = rho_v/rho_dry.
    True pv = pv * P0'''
    
    pv = p*qv/(mr_ratio+qv)
    
    return pv


def cp_moist(qv, qc):
    ''' calculate effecitve specific heat 
    capacity of parcel with water vapour 
    qv=m_v/m_dry and liquid water qc= m_l/m_dry
    (notice m_dry not m_tot). Enthalpy change of
    parcel per unit (total) mass = dh. 
    (1+qv+qc)dh = cp_eff*dtemp + latent_v*dqv.
    true cp_moist = Cp0 * cp_moist.'''

    cp_moist = cp_dry + cp_v*qv + c_l*qc 
    
    cp = cp_dry * (1+qv*cp_v/cp_dry)/(1+qv)
    
    return cp_moist



def saturation_pressure(T, Temp0, P0):
    ''' Calculate the equilibrium vapor pressure 
    of water over liquid water ie. the
    saturation pressure (psat). Equation taken from
    typhon.physics.thermodynamics.e_eq_water_mk.
    Real temp = T*Temp0, dimensionless psat = psat/P0'''

    T = T*Temp0                               # real T [K]
    
    if np.any(T <= 0):
        err = 'T must be larger than 0K.'
        raise ValueError(err)

    lnpsat = (54.842763                    # ln(psat) [Pa]
         - 6763.22 / T
         - 4.21 * np.log(T)
         + 0.000367 * T
         + np.tanh(0.0415 * (T - 218.8))
         * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025 * T))

    return np.exp(lnpsat)/P0               # dimensionless psat


############################################################################





######################## ODE Solver Functions ########################
def dp_dt_profile(t, w, Z0, Temp0, P0):
    '''dimensionless pressure profile used for 
    pressure change of rising parcel over time 
    (dP/dt = dP/dz * dz/dt. Here, dp/dz is chosen 
    to be that for a dry hydrostatic atmosphere 
    with constant lapse rate. 
    true dP/dt = P0/T0 * dp '''

    tempg = 273.15/Temp0             #dimless tempg
    pg = 100000/P0                   #dimless pg
    zg = 0/Z0                        #dimless zg
    
    lpsrate = 0.0062/Temp0*Z0     #dimless lapse rate
    
    profile = 1 - lpsrate/tempg*(w*t-zg)
    profile = profile**(G/(Rgas_dry*0.0062)-1)

    const = Z0*G/(Rgas_dry*Temp0)
    dp = -const*pg/tempg*profile
    
    return dp


    

def dtemp_expansion(temp, p, dp, qv, qc, Cp0):
    ''' calculate dT/dt due to pressure 
    change dp for adiabatic process
    (no heat loss) of parcel that has 
    water vapour mass mixing ratio 
    (m_v/m_dry) = qv and liquid water mass
    mixing ratio (m_c/m_dry) = qc
    Assumes instantaneous volume change 
    of parcel to change pressure according
    to dp_dt_profile. 
    True dTemp/dt = dtemp * Temp0/t0 '''
    
    pv = qv2pv(qv, p) 
    rho_d = (p-pv)/temp                            # dimless density of dry parcel
    
    dtemp = rgas_dry/(rho_d*cp_moist(qv, qc)) * dp
     
    return dtemp





def diffusion_factors(temp, p, psat, Temp0, P0, R0, Rho0, T0):
    ''' calculate Fkl and Fdl heat and vapour
    diffusion factors in equation for 
    radial growth of droplets according to
    eq.7.27 lohmann, luond and mahrt intro 2 
    clouds textbook '''

    temp, p, psat = temp*Temp0, p*P0, psat*P0

    # diffusional growth equation factors
    Rhermk = 7.11756e-5*temp + Thermk273                     #[eq.7.24 lohmann intro 2 clouds]
    Diffuse_v = Diffuse_v273 / p * temp**1.94
    
    Fkl = (Latent_v/(Rgas_v*temp) -1)*Latent_v/(Rhermk*temp) 
    Fdl = Rgas_v*temp/(Diffuse_v*psat)
    
    F0 = T0/(Rho0*R0**2)
    fkl = Fkl/F0                                          # dimless diffusion factors
    fdl = Fdl/F0 
    
    return fkl, fdl


        
        
def diffusion_growth(drops, temp, p, qv, qc, t):
    ''' diffusion growth of droplets by condensation
    given temperature, pressure and qv.
    NOTE! shrinking by evaporation not includeded, nor
    are ventillation effects fv(r) (see Seifert & stevens 2010)'''
    
    r = np.asarray(same_attr_from_objs(drops, 'r'))
    eps = np.asarray(same_attr_from_objs(drops, 'eps'))
    
    # supersaturation of parcel
    psat = saturation_pressure(temp, Temp0, P0)
    s_ratio = qv2pv(qv, p)/psat                                      # ambient supersaturation ratio 
    a, b = kohfactors(drops, temp, Temp0, R0)
    
    # radial growth/shrink of droplets
    fkl, fdl = diffusion_factors(temp, p, psat, Temp0, P0, R0, Rho0, T0)                             
    #dr = (s_ratio-1 -a/r +b/r**3) / (rho_l * (fkl+fdl) * r)           # [eq.7.27 lohmann intro 2 clouds]
    dr = (s_ratio-1) / (rho_l * (fkl+fdl) * r)           # [eq.7.27 lohmann intro 2 clouds]
       
    dry_r = np.asarray(same_attr_from_objs(drops, 'dry_r'))
    if (r <= dry_r).any():                                       # if droplets are dry, do not shrink further
        dry_drs = np.where(r<=dry_r, dr, 0.0)
        dr = np.where(dry_drs<=0.0, 0.0, dr)
        r = np.where(r <= dry_r, dry_r, r)
    
    # resultant change to temp, qv and qc of parcel
    dm = 4*np.pi*rho_l*(r**2)*dr                             # mass of water condensed onto each drop
    tot_drhov = np.sum(dm*eps)*R0**3*N0                      # change in density of water in volume
    dqc = tot_drhov/rho_dry
    dqv = -dqc
    dtemp_c = -latent_v/(cp_moist(qv, qc)) * dqv
    
    
    return dr, dtemp_c, dqv, dqc








def rising_parcel(t, y, w, drops, cond, printt=False):
    ''' differential equations for rising parcel 
    (thermo)dynamics with condensation onto droplets
    if supersaturation > 0.
    Parcel rises, expands adiabatically and then
    condenses water if qv > q_sat. Rate of condensation
    determined by growth rate of droplets (supersaturation
    is possible if growth slower than rise in relative
    humidity)'''
    
    print('t = {:.3f}'.format(t), end='\r')
    
    p, temp, qv, qc = y[:4]
    
    ### dynamics of droplets (position evolution)
    same_attr_from_objs(drops, 'r', y[4:len(drops)+4])                       # update radius of superdroplets
    dxyz_drops = velocities(drops, Rho0, R0, W0)
    
    ### thermodynamics of parcel and droplets
    #tau = 50*60/T0                                                              
    #dz = w*np.sin(2*np.pi*t/tau)
    dz = w                                                                   # parcel motion (upwards)                                              
    dp = dp_dt_profile(t, dz, W0*T0, Temp0, P0)                                # expansion of parcel according to dp_dz profile
    dtemp = dtemp_expansion(temp, p, dp, qv, qc, Cp0)                             # Temp change solely due to expansion (no condensation)
   
    # condensation (by diffusional growth of droplets)
    if cond:
        dr, dtemp_c, dqv, dqc = diffusion_growth(drops, 
                                temp, p, qv, qc, t)
    else:
        dr = np.zeros(len(drops))
        dtemp_c, dqv, dqc = 0, 0, 0

    
    dtemp+=dtemp_c
    dy = [dp, dtemp, dqv, dqc]
    dy += list(dr)
    dy += dxyz_drops

    #dy.extend(dr)
    #dy.extend(dxyz_drops)

    return dy