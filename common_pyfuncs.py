import numpy as np


######## functions for converting c++ values into python  ########

def read_cpp_into_floats(filename):
  """make dictionary of value: float for 
  doubles and ints in c++ file. Also make
  dictionary of notfloats for values that
  couldn't be converted"""
  constants = {}
  notfloats = {}
  with open(filename) as file:
      rlines=[]
      filelines = file.readlines()
      for line in filelines:
        ind1 = line.find("=")
        ind2 = line.find(";")
        if ind2 != -1 and ind1 != -1:
          line = line[:ind2+1]
          if "double" in line or "int" in line: 
            rlines.append(line)
      
      for line in rlines:
        if line[:13] == "const double ": x=13

        elif line[:10] == "const int ": x=10
        elif line[:4] == "int ": x=4
        elif line[:7] == "double ": x=7
        ind1 = line.find("=")
        ind2 = line.find(";")
        indx = line.find(" ", x) 
        symbol = line[x:indx]
        try: 
          float(line[ind1+2:ind2])
          constants[symbol] = float(line[ind1+2:ind2])
        except ValueError:
          notfloats[symbol] = line[ind1+2:ind2]
  
  print("---- Constants read from ", filename, "-----")
  for c in constants:
    print(c, "=", constants[c])
  print("---------------------------------------------")
  print("---- Not floats read from ", filename, "-----")
  for st in notfloats:
    print(st, "=", notfloats[st])
  print("---------------------------------------------")

  return constants, notfloats


###################################################################


######## functions for reading .csv files into python
### and converting into variables with dimensions ########

def get_soldata(sol_filename, TIME0, P0, TEMP0):
    
    #### Load data from .csv file ###
    with open(sol_filename) as file_name:
        t, p, temp, qv, qc = np.loadtxt(file_name, delimiter=",", comments="/*", unpack=True)

    print("--- Raw Data Shapes ---")
    print("variables: t, p, temp, qv, qc")
    print(t.shape, p.shape, temp.shape, qv.shape, qc.shape)

    print("--- Non Dimensional Max/Mins of Data ---")
    print("time:", np.amin(t), np.amax(t))
    print("p:", np.amin(p), np.amax(p))
    print("temp:", np.amin(temp), np.max(temp))
    print("(qv, qc)", (np.amin(qv), np.amin(qc)), (np.amax(qv), np.amax(qc)))
    
    
    return t*TIME0, p*P0, temp*TEMP0, qv, qc      ### (Re-)dimensionalise and return


def get_SDdata(SDsol_filename, nsupers, R0, RHO0):

    ### Load data from SD .csv file ###
    with open(SDsol_filename) as file_name:
        drops = np.loadtxt(file_name, delimiter=",", comments="/*", unpack=False)
    eps = drops[:,0:nsupers]
    r = drops[:,nsupers:2*nsupers]
    m_sol = drops[:,2*nsupers:]


    print("--- Raw SD Data Shapes ---")
    print("variables: eps, r, m_sol")
    print(eps.shape, r.shape, m_sol.shape)

    print("--- Non Dimensional Max/Mins of Data ---")
    print("droplet eps:", np.amin(eps), np.amax(eps))
    print("droplet r:", np.amin(r), np.amax(r))
    print("droplet m_sol:", np.amin(m_sol), np.amax(m_sol))
    

    return eps, r*R0, m_sol*RHO0*R0**3           ### (Re-)dimensionalise and return





######## functions for plotting droplet lnR distributions  ########


def logr_distribution(rspan, nbins, data, wghts, 
           ax=False, step=False, lab=None, c='k', 
                 xlab='ln(r /\u03BCm)', ylab=None, 
                        perlnR=False, smooth=False):
  ''' plot data distribution against logr (using np.histogram to
  get frequency of a particular value of data that falls in 
  each lnr -> lnr + dlnr  bin). Apply guassian kernel smoothing if wanted '''

  # create lnr bins (linearly spaced in lnr)
  hedgs = np.linspace(np.log(rspan[0]), np.log(rspan[1]), nbins+1)             # edges to lnr bins
  hwdths = hedgs[1:]- hedgs[:-1]                               # lnr bin widths
  hcens = (hedgs[1:]+hedgs[:-1])/2                             # lnr bin centres

  # get number frequency in each bin
  hist, hedgs = np.histogram(np.log(data), bins=hedgs, 
                  weights=wghts, density=None)
  if perlnR == True: # get frequency / bin width
      hist = hist/hwdths


  if smooth:
    hist, hcens = gaussian_kernel_smoothing(hist, hcens, smooth)

  if ax:
      if step:
          ax.step(hcens, hist, label=lab, where='mid', color=c)
      else:
          ax.bar(hcens, hist, hwdths, label=lab, color=c)

      ax.legend()
      ax.set_ylabel(ylab)
      if xlab:
          ax.set_xlabel(xlab)
          ax.xaxis.tick_top() 
          ax.xaxis.set_label_position('top')
      else:
          ax.set_xticks([])

  return hist, hedgs




def linear_twinax(ax, lnr, eps):
    ''' linear x axis for lognormal 
    radius distribution plots'''
     
    axb = ax.twiny()
   
    axb.plot(np.e**lnr*1e6, eps, alpha=0)
    axb.set_xscale('log')
    axb.set_xlabel('radius, r /\u03BCm)')

    axb.xaxis.tick_bottom() 
    axb.xaxis.set_label_position('bottom') 

    return axb



def gaussian_kernel_smoothing(hist, hcens, sig):

    for h in range(len(hist)):
        kernel = 1/(np.sqrt(2*np.pi)*sig)*np.exp(-(hcens - hcens[h])**2/(2*sig**2))
        kernel = kernel/np.sum(kernel)
        hist[h] = np.sum(hist*kernel)  

    hist = np.where(hist<1e-16, 0, hist)  

    return hist, hcens



###################################################################
 








######## functions for superdroplet attributes  ########


class Common2AllSuperdrops():
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
    
    print("---- Superdrop Properties -----")
    print("nsupers =", self.nsupers)
    print("parcel volume =", self.VOL, "m^3")
    print("RHO_L =", self.RHO_L, "Kg/m^3")
    print("RHO_SOL =", self.RHO_SOL, "Kg/m^3")
    print("MR_SOL =", self.MR_SOL, "Kg/mol")
    print("IONIC =", self.IONIC)
    print("-------------------------------")
 


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