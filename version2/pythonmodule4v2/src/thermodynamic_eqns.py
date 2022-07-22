import numpy as np

def vapour_pressure(p, qv, Mr_ratio):
    
    pv = qv*p/(Mr_ratio + qv) 
    
    return pv




def saturation_pressure(T):
  ''' Calculate the equilibrium vapor pressure 
  of water over liquid water ie. the
  saturation pressure (psat [Pa]). Equation taken from
  typhon.physics.thermodynamics.e_eq_water_mk.'''
  
  if np.any(T <= 0):
      err = 'T must be larger than 0K.'
      raise ValueError(err)

  lnpsat = (54.842763                    # ln(psat) [Pa]
        - 6763.22 / T
        - 4.21 * np.log(T)
        + 0.000367 * T
        + np.tanh(0.0415 * (T - 218.8))
        * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025 * T))

  return np.exp(lnpsat)               # psat [Pa]



def dry_pot_temp(Temp, P, qv, CONSTS, MCONSTS):
    ''' calculate potential Temperature [K]
    assuming moist (unsaturated) air with
    vapour content qv '''

    Cpdry = CONSTS["CP_DRY"]
    Cpv = CONSTS["CP_V"], 
    Rgasdry = MCONSTS["RGAS_DRY"]
    Rgasv = MCONSTS["RGAS_V"]

    Cp =  Cpdry * (1+qv*Cpv/Cpdry)/(1+qv)
    Rgas = Rgasdry *(1+qv*Rgasv/Rgasdry)/(1+qv)
    
    Theta = Temp*(P[0]/P)**(Rgas/Cp)
    
    return Theta



def dry_adiabat(p, temp, CONSTS, MCONSTS):

    gamma = (MCONSTS["RGAS_DRY"]/CONSTS["CP_DRY"]) 
    dry_adia = temp[0]*(p/p[0])**gamma            # dry adiabatic temp

    dry_adia_theta = dry_adia*(p[0]/p)**gamma     # dry adiabatic theta (=const)


    return dry_adia, dry_adia_theta



def relative_humidity(p, temp, qv, MCONSTS):

    pv = vapour_pressure(p, qv, MCONSTS["Mr_ratio"])
    psat = saturation_pressure(temp)
    relh = pv/psat
    
    qsat = MCONSTS["Mr_ratio"] * psat/(p-pv) 
    s = qv/qsat - 1  


    return relh, s



def moist_static_energy(Z, Temp, qv, CONSTS):
    ''' calculate the moist static energy J/m^3
        (not assuming dry air cp) '''
    
    G = CONSTS["G"]
    Latent_v = CONSTS["LATENT_V"]
    Cp_dry = CONSTS["CP_DRY"]

    return  G*Z + Latent_v*qv + Cp_dry*Temp




