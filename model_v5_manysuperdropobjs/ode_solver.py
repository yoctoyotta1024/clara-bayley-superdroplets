### FUNCTIONS FOR ODE SOLVER
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import g, R
#from typhon.physics import e_eq_water_mk as saturation_pressure

from constants import *
from superdroplets_as_objects import *
thermk273 = 0.023822892 - 7.11756e-5*273.15             #[eq.7.24 lohmann intro 2 clouds]
diffuse_v273 = 2.11e-5 * 101325 / (273.15**1.94)         #[eq.7.26 lohmann intro 2 clouds]



### ODE input dp/dt = dz*dp_dz_profile
def dp_dz_profile(z, pprofile=False):
    '''pressure profile used for pressure change 
    of parcel over time (dP/dt = dP/dz * dz/dt.
    Here, dp/dz is chosen to be that for a dry
    hydrostatic atmosphere with constant lapse rate '''

    pg, zg = 100000, 0
    lps_rate, tempg = 0.0062, 273.15
    profile = 1 - lps_rate/tempg*(z-zg)
    temp_surr = tempg*profile
    p_surr = pg * profile**(g/(rgas_dry*lps_rate))
    
    if pprofile:
        return p_surr
    else:
        dp_dz = -g*p_surr/(rgas_dry*temp_surr)      
        return dp_dz


    
    
### ODE solver functions
def dtemp_expansion(temp, p, dp, qv, qc):
    ''' calculate dT/dt due to pressure 
    change dp for adiabatic process
    (no heat loss) of parcel that has 
    water vapour mass mixing ratio 
    (m_v/m_dry) = qv and liquid water mass
    mixing ratio (m_c/m_dry) = qc
    Assumes instantaneous volume change 
    of parcel to change pressure according
    to dp_dz profile. '''
    
    pv = qv2pv(qv, p) 
    rho_d = (p-pv)/(rgas_dry*temp)                                    # dry density of parcel (~ total density)
    
    dtemp = 1/cp_moist(qv, qc) * (dp/rho_d)
     
    return dtemp





def diffusion_factors(temp, p, psat):
    ''' calculate Fkl and Fdl heat and vapour
    diffusion factors in equation for 
    radial growth of droplets according to
    eq.7.27 lohmann, luond and mahrt intro 2 
    clouds textbook '''

    # diffusional growth equation factors
    thermk = 7.11756e-5*temp + thermk273                     #[eq.7.24 lohmann intro 2 clouds]
    diffuse_v = diffuse_v273 / p * temp**1.94
    
    fkl = (latent_v/(rgas_v*temp) -1)*latent_v/(thermk*temp) 
    fdl = rgas_v*temp/(diffuse_v*psat)

    
    return fkl, fdl


        
        
def diffusion_growth(drops, temp, p, qv, qc, t):
    ''' diffusion growth of droplets by condensation
    given temperature, pressure and qv.
    NOTE! shrinking by evaporation not includeded, nor
    are ventillation effects fv(r) (see Seifert & stevens 2010)'''
    
    r = np.asarray(same_attr_from_objs(drops, 'r'))
    eps = np.asarray(same_attr_from_objs(drops, 'eps'))
    
    # supersaturation of parcel
    psat = saturation_pressure(temp)
    s_ratio = qv2pv(qv, p)/psat                                      # ambient supersaturation ratio 
    a, b = kohfactors(drops, temp)
    
    # radial growth/shrink of droplets
    fkl, fdl = diffusion_factors(temp, p, psat)                             
    #dr = (s_ratio-1 -a/r +b/r**3) / (rho_l * (fkl+fdl) * r)           # [eq.7.27 lohmann intro 2 clouds]
    dr = (s_ratio-1) / (rho_l * (fkl+fdl) * r)           # [eq.7.27 lohmann intro 2 clouds]
    
    dry_r = np.asarray(same_attr_from_objs(drops, 'dry_r'))
    if (r < dry_r).any():                                       # if droplets are dry, do not shrink further
        dry_drs = np.where(r<dry_r, dr, 0.0)
        dr = np.where(dry_drs<0.0, 0.0, dr)
        r = np.where(r < dry_r, dry_r, r)
    
    # resultant change to temp, qv and qc of parcel
    dm = rho_l*4*np.pi*(r**2)*dr                                # mass of water condensed onto each drop
    tot_drhov = np.sum(dm*eps)                                  # change in density of water in volume
    dqc = tot_drhov/rho_dry
    dqv = -dqc
    dtemp_c = -latent_v/cp_moist(qv, qc) * dqv
    
    
    
    return dr, dtemp_c, dqv, dqc








def rising_parcel(t, y, w, drops, pclsize, cond, printt=False):
    ''' differential equations for rising parcel 
    (thermo)dynamics with condensation onto droplets
    if supersaturation > 0.
    Parcel rises, expands adiabatically and then
    condenses water if qv > q_sat. Rate of condensation
    determined by growth rate of droplets (supersaturation
    is possible if growth slower than rise in relative
    humidity)'''
    
    print('t = {:.1f}s'.format(t), end='\r')
    
    z, temp, p, qv, qc = y[:5]
    
    ### dynamics of droplets (position evolution)
    same_attr_from_objs(drops, 'r', y[5:len(drops)+5])                       # update radius of superdroplets
    dxyz_drops = velocities(drops)
    
    ### dynamics of parcel and droplets
    #tau = 50*60                                                                 # tau_w = 50min
    #dz = w*np.sin(2*np.pi*t/tau)
    dz = w                                                                       # parcel motion (upwards)     
    
    
    ### thermodynamics of parcel and droplets
    dp = dp_dz_profile(z) * dz                                                   # expansion of parcel according to dp_dz profile
    dtemp = dtemp_expansion(temp, p, dp, qv, qc)                                 # Temp change solely due to expansion (no condensation)
   
    # condensation (by diffusional growth of droplets)
    if cond:
        dr, dtemp_c, dqv, dqc = diffusion_growth(drops, 
                                temp, p, qv, qc, t)
    else:
        dr = np.zeros(len(drops))
        dtemp_c, dqv, dqc = 0, 0, 0

    
    dtemp+=dtemp_c
    dy = [dz, dtemp, dp, dqv, dqc]
    dy += list(dr)
    dy += dxyz_drops

    #dy.extend(dr)
    #dy.extend(dxyz_drops)


    return dy