import numpy as np
import matplotlib.pyplot as plt
import csv

import src.convert_cxx_2_pyfloats as cxx2py
import src.volexponential as volexp
import src.lognormal as lognorm
import src.plot_lnr_distrib as pltlnr
import src.generic_axfuncs as axfuncs

plt.rcParams.update({'font.size': 14})

##############################################
### droplet properties for initialisation ### 
INITDROPSCSV = "./build/dimlessSDinit.csv"

### settings for creating distribution from lognormal in R space distribution
use_lognormal    = False
if use_lognormal:
  rspan            = [1e-8, 1e-5]                  # initial range of droplet radii [m]
  # mus             = [0.075e-6]                     # [m] geometric mean droplet radius
  # sigs            = [1.5]                    
  # n_as            = [1e9]                          # [m^-3] total no. concentration of droplets          
  mus              = [0.02e-6, 0.2e-6, 3.5e-6]               
  sigs             = [1.55, 2.3, 2]                    
  n_as             = [1e6, 0.3e6, 0.025e6]   


### settings for distirbution from exponential in droplet volume
use_volexponential = True
if use_volexponential:
  r_a             = 30.531e-6                   # [m] radius of peak to volume dist
  n_a             = 2**(23)                     # [m^-3] total no. concentration of real droplets          
  rho_expo        = 1

if use_volexponential == use_lognormal:
  raise ValueError("same boolean for use_volexponential and use_lognormal")
#############################################

##############################################
###  CONSTANTS!! MAKE SURE THESE ARE SAME AS IN  ###
###  constants.hpp and init.hpp FOR ODE SOLVER!!  ###

### read in constants from .hpp files
CONSTS, notfloats = cxx2py.read_cpp_into_floats("./claras_SDconstants.hpp")
INITS, notfloats2 = cxx2py.read_cpp_into_floats("./claras_SDinit.hpp") 
INITS["nsupers"] = int(INITS["NSUPERS"])

RGAS_DRY   = CONSTS["RGAS_UNIV"]/CONSTS["MR_DRY"]  # specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
P0         = CONSTS["P0"]                          # pressure [Pa]
TEMP0      = CONSTS["TEMP0"]                       # temperature [K]
RHO0       = P0/(RGAS_DRY*TEMP0)                   # density [Kg/m^3]
R0         = CONSTS["R0"]                          # droplet radius lengthscale [m]
#VOL0         = CONSTS["VOL0"]                      # droplet multiplicity [m^-3]
Rho_sol    =  CONSTS["RHO_SOL"]/RHO0               # dimensionless solute density []
nsupers    = int(INITS["nsupers"])                 # no. of distinct superdrops (different initial radii (evenly spaced between ln(rspan))
VOL = INITS["DROPVOL"]                         #parcel volume [m]

print("---- Additional Constants Derived -----")
print("RGAS_DRY", "=", RGAS_DRY)
print("RHO0", "=", RHO0)
print("Rho_sol", "=", Rho_sol)
print("nsupers", "=", nsupers)
print("---------------------------------------------")

##############################################






#############################################
###       create initial superdroplet      ###
###  conditions from chosen distribiution  ###

if use_lognormal:
  ### each superdroplet is different radius with 
  ### epsilon determined from lognormal distribution
  for d in range(len(mus)):
    r0, eps1 = lognorm.dimless_lnr_dist(rspan, nsupers, n_as[d], mus[d], sigs[d])[0:2]

    if d==0:
      eps = eps1
    else:                        
      eps = eps+eps1    
  m_sol = Rho_sol*4/3*np.pi*r0**3                   # assuming initially dry areosol droplets (dimless mass_solute)



elif use_volexponential:
  ### each superdroplet has same epsilon. Radii 
  ### superdroplets is random sample of exponential
  ### volume (or mass) distribution
  r0, eps = volexp.dimless_randomexpo_dist(nsupers, 
                  VOL, n_a, r_a, rho=rho_expo, R0=R0)
  m_sol = Rho_sol*4/3*np.pi*r0**3                   # assuming initially dry areosol droplets (dimless mass_solute)

print(r0[0], m_sol[0])
print(Rho_sol, RHO0)
print(Rho_sol*4/3*np.pi)
#############################################


#############################################
###       write initial superdroplet      ###
###        conditions to .csv file        ###
print("Writing inital droplet distribution to: ./init_superdroplets.csv")

with open('./'+INITDROPSCSV, 'w', encoding='UTF8') as f:
  writer = csv.writer(f)
  writer.writerow(["/* Initial Dimensionless Superdroplets Data"])
  writer.writerow([" VOL = "+str(VOL)+"[m^3]", "  R0 = "+str(R0)+"[m]", 
                                "  RHO_SOL = "+str(Rho_sol)+"[Kg m^-3]"])
  space = "           "                              
  writer.writerow(["eps", space+"r", space+"m_sol */"])
  for i in range(nsupers):
    row = ["{:.14f}".format(eps[i]), "{:.14f}".format(r0[i]), "{:.14f}".format(m_sol[i])] 
    writer.writerow(row)

##############################################




##############################################
###        plot initial distribution       ###
###              as histogram              ###


print('Plotting Initial Droplet No. Concentration Distribution:\n')
fig, ax = plt.subplots(figsize=(9,6))


if use_lognormal:
  nbins = nsupers
elif use_volexponential:
  nbins = 50
  rspan = [np.amin(r0*R0), np.amax(r0*R0)]
  #wghts = [1]*nsupers
wghts = eps/VOL
axfuncs.linear_twinax(ax, np.log(r0*R0), wghts)
hist, hedgs = pltlnr.logr_distribution(rspan, nbins, r0*R0, wghts, ax=ax, 
  ylab="No. Conc. Real Droplets /m$^{-3}$", lab="initial superdroplets", 
  c='C0', perlnR=False, smooth=False)
       
if use_lognormal:
  # plot lognormal curves with x10 more bins for reference
  for d in range(len(mus)):
      pltr, plteps, pltwdths = lognorm.dimless_lnr_dist(rspan,                     
                      nsupers*10, n_as[d], mus[d], sigs[d])[0:3]
      ax.plot(np.log(pltr)+np.log(R0), plteps*10/VOL, 
                      color='C'+str(d+1), label='lognormal '+str(d)) 

ax.legend()
plt.tight_layout()
print("Total No. real droplets = ", np.sum(eps), "m^-3")
print("--> Conc. real droplets = ", np.sum(hist), "m^-3")
print("No. superdroplets = ", nsupers)
print("Average multiplicity = ", np.sum(hist)*VOL/nsupers)
plt.savefig("./build/initial_superdroplet_dist.png", dpi=400, bbox_inches="tight")
print("Initial Distribution saved as: ./initial_superdroplet_dist.png")
plt.show()



if use_volexponential:
  print('Plotting Initial Mass Density Distribution:\n')

  fig, ax = plt.subplots(figsize=(9,6))

  nbins = 100
  rspan = [np.amin(r0*R0), np.amax(r0*R0)]
  smoothsig = 0.62*nsupers**(-1/5)
  
  mass = R0**3*RHO0*m_sol*(1-CONSTS["RHO_L"]/CONSTS["RHO_SOL"])
  mass += 4/3.0*np.pi*((r0*R0)**3)*CONSTS["RHO_L"]                                 
  wghts = eps*mass/VOL

  ax2 = axfuncs.linear_twinax(ax, np.log(r0*R0), wghts)
  hist, hedgs = pltlnr.logr_distribution(rspan, nbins, r0*R0, wghts, ax=ax, 
    ylab="g(lnR) /Kg m$^{-3}$ / unit lnR", lab="initial superdroplets", 
    c='C0', perlnR=True, smooth=False)
        
  ax.legend()
  plt.tight_layout()
  plt.show()
  #############################################
