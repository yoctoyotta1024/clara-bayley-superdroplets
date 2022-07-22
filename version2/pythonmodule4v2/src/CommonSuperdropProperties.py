import numpy as np


######## functions for superdroplet attributes  ########

class CommonSuperdropProperties():
  '''Parent class for all superdroplets. Each 
    Superdrop instance has these properties'''
    
   
  def __init__(self, nsupers, VOL, RHO_L, RHO_SOL, MR_SOL, IONIC):

    # Common attributes of all superdroplets
    self.nsupers = nsupers
    self.VOL = VOL
    self.RHO_L = RHO_L                                # density of liquid in droplets (=density of water at 300K) [Kg/m^3]
    
    # droplet solute properties
    self.RHO_SOL = RHO_SOL                           # density of (dry) solute [Kg/m^3]
    self.MR_SOL = MR_SOL                             # Mr of solute [g/mol]
    self.IONIC = IONIC                               # degree ionic dissociation (van't Hoff factor)
  

  def print_properties(self):
    print("\n---- Superdrop Properties -----")
    print("nsupers =", self.nsupers)
    print("parcel volume =", self.VOL, "m^3")
    print("RHO_L =", self.RHO_L, "Kg/m^3")
    print("RHO_SOL =", self.RHO_SOL, "Kg/m^3")
    print("MR_SOL =", self.MR_SOL, "Kg/mol")
    print("IONIC =", self.IONIC)
    print("-------------------------------\n")
 


  def rhoeff(self, r, m_sol):
    ''' calculates effective density [Kgm^-3] of 
  droplet such that mass_droplet, m = 4/3*pi*r^3 * rhoeff
  taking into account mass of liquid and mass of
  solute assuming solute occupies volume it
  would given its (dry) density, RHO_SOL. '''

    solfactor = 3*m_sol/(4.0*np.pi*(r**3))
    rhoeff = self.RHO_L + solfactor*(1-self.RHO_L/self.RHO_SOL)

    return rhoeff


  def vol(self, r):
    ''' volume of droplet [m^3] '''
  
    return 4.0/3 * np.pi * r**3
      
      
  def m(self, r, m_sol):
    ''' total mass of droplet (water + (dry) areosol),
     m =  4/3*pi*rho_l**3 + m_sol(1-rho_l/rho_sol) 
    ie. m = 4/3*pi*rhoeff*R**3 '''
        
    m = m_sol*(1-self.RHO_L/self.RHO_SOL)
    m += 4/3.0*np.pi*(r**3)*self.RHO_L                                  

    return m


  def m_w(self, r, m_sol):
    ''' mass of only water in droplet '''
        
    v_w = 4/3.0*np.pi*(r**3)
    v_w += -m_sol/self.RHO_SOL
    
    return self.RHO_L*v_w