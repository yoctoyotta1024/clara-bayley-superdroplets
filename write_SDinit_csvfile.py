import numpy as np
import matplotlib.pyplot as plt
import csv

from common_pyfuncs import read_cpp_into_floats, logr_distribution, linear_twinax

plt.rcParams.update({'font.size': 14})


##############################################
### droplet properties for initialisation ### 
INITDROPSCSV = "dimlessSDinit.csv"

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


if use_volexponential == use_lognormal:
  raise ValueError("same boolean for use_volexponential and use_lognormal")
#############################################


##############################################
###  CONSTANTS!! MAKE SURE THESE ARE SAME AS IN  ###
###  constants.hpp and init.hpp FOR ODE SOLVER!!  ###

### read in constants from .hpp files
CONSTS, notfloats = read_cpp_into_floats("./constants.hpp")
INITS, notfloats2 = read_cpp_into_floats("./init.hpp") 
INITS["nsupers"] = int(INITS["nsupers"])

RGAS_DRY   = CONSTS["RGAS_UNIV"]/CONSTS["MR_DRY"]  # specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
P0         = CONSTS["P0"]                          # pressure [Pa]
TEMP0      = CONSTS["TEMP0"]                       # temperature [K]
RHO0       = P0/(RGAS_DRY*TEMP0)                   # density [Kg/m^3]
R0         = CONSTS["R0"]                          # droplet radius lengthscale [m]
#VOL0         = CONSTS["VOL0"]                      # droplet multiplicity [m^-3]
Rho_sol    =  CONSTS["RHO_SOL"]/RHO0               # dimensionless solute density []
nsupers    = int(INITS["nsupers"])                 # no. of distinct superdrops (different initial radii (evenly spaced between ln(rspan))
VOL = INITS["iVOL"]                         #parcel volume [m]

print("---- Additional Constants Derived -----")
print("RGAS_DRY", "=", RGAS_DRY)
print("RHO0", "=", RHO0)
print("Rho_sol", "=", Rho_sol)
print("nsupers", "=", nsupers)
print("---------------------------------------------")

##############################################



    
#############################################
###      functions for getting droplet    ###
###          radii distribution           ###

def dimless_lnr_dist(rspan, nbins, n_a, mu, sig):
  ''' returns 'nbins' no. of samples from 
  lognormal distribution with even spacing
  in lnr space between ln(r[0]) and ln(rs[1])'''
  
  edgs = np.linspace(np.log(rspan[0]), np.log(rspan[1]), nbins+1)        # edges to linearly spaced lnr bins
  wdths = edgs[1:]- edgs[:-1]                                      # lnr bin widths
  lnr = (edgs[1:]+edgs[:-1])/2                                     # lnr bin centres
  lnnorm = lnnormal_dist(np.e**lnr, n_a, mu, sig)                  # lognormal values

  ### make radii and no. concentrations dimensionless
  lnr = lnr - np.log(R0)                                           
  edgs = edgs - np.log(R0)
  wdths = edgs[1:]- edgs[:-1]                                     # lnr bin widths
  eps = np.asarray([int(n) for n in lnnorm*wdths*VOL])     # dimless multiplicities from lognorm dist
  

  return np.e**(lnr), eps, wdths, edgs



def lnnormal_dist(r, n_a, mu, sig):
  ''' return number concentration (n [m^-3])
  of particles with radius r in bin of unit length 
  on logarithmic scale based on monomodal 
  lognormal distribution. ie. No particles per cm"^-3
  in bin at radius r of width delta(lnr), n =
  dn_dlnr * delta(lnr)'''
  
  sigtilda = np.log(sig)
  mutilda = np.log(mu)
  norm = n_a/(np.sqrt(2*np.pi)*sigtilda)
  exponent = -(np.log(r)-mutilda)**2/(2*sigtilda**2)
  
  dn_dlnr = norm*np.exp(exponent)                                 # eq.5.8 [lohmann intro 2 clouds]
  
  
  return dn_dlnr



def dimless_randomexpo_dist(nsupers, VOL, n_a, r_a, rho=1):
  ''' return nsupers superdroplets with eps = n_a*vol/nsupers
  and radius taken from random rample of exponential volume
  (or mass) distribution '''
  
  eps = np.full(nsupers, n_a*VOL/nsupers)

  vol_a = 3/4*np.pi*rho*r_a**3
  vols = np.random.exponential(vol_a, nsupers)
  r = (4*vols/(3*np.pi*rho))**(1/3)
 
  return r/R0, eps
#############################################




#############################################
###       create initial superdroplet      ###
###  conditions from chosen distribiution  ###

if use_lognormal:
  ### each superdroplet is different radius with 
  ### epsilon determined from lognormal distribution
  for d in range(len(mus)):
    r0, eps1 = dimless_lnr_dist(rspan, nsupers, n_as[d], mus[d], sigs[d])[0:2]

    if d==0:
      eps = eps1
    else:                        
      eps = eps+eps1    
  m_sol = Rho_sol*4/3*np.pi*r0**3                   # assuming initially dry areosol droplets (dimless mass_solute)



elif use_volexponential:
  ### each superdroplet has same epsilon. Radii 
  ### superdroplets is random sample of exponential
  ### volume (or mass) distribution
  r0, eps = dimless_randomexpo_dist(nsupers, VOL, n_a, r_a)
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
linear_twinax(ax, np.log(r0*R0), wghts)
hist, hedgs = logr_distribution(rspan, nbins, r0*R0, wghts, ax=ax, 
  ylab="No. Conc. Real Droplets /m$^{-3}$", lab="initial superdroplets", 
  c='C0', perlnR=False, smooth=False)
       
if use_lognormal:
  # plot lognormal curves with x10 more bins for reference
  for d in range(len(mus)):
      pltr, plteps, pltwdths = dimless_lnr_dist(rspan,                     
                      nsupers*10, n_as[d], mus[d], sigs[d])[0:3]
      ax.plot(np.log(pltr)+np.log(R0), plteps*10/VOL, 
                      color='C'+str(d+1), label='lognormal '+str(d)) 

ax.legend()
plt.tight_layout()
print("Total No. real droplets = ", np.sum(eps), "m^-3")
print("--> Conc. real droplets = ", np.sum(hist), "m^-3")
print("No. superdroplets = ", nsupers)
print("Average multiplicity = ", np.sum(hist)*VOL/nsupers)
plt.savefig("initial_superdroplet_dist.png", dpi=400, bbox_inches="tight")
print("Initial Distribution saved as: ./initial_superdroplet_dist.png")
plt.show()



if use_volexponential:
  print('Plotting Initial Mass Density Distribution:\n')

  fig, ax = plt.subplots(figsize=(9,6))

  nbins = 50
  rspan = [np.amin(r0*R0), np.amax(r0*R0)]
  mass = 4/3*np.pi*(r0*R0)**3 * CONSTS["RHO_L"]
  wghts = eps*mass/VOL
  smoothsig = 0.62*nsupers**(-1/5)

  ax2 = linear_twinax(ax, np.log(r0*R0), wghts)
  hist, hedgs = logr_distribution(rspan, nbins, r0*R0, wghts, ax=ax, 
    ylab="g(lnR) /Kg m$^{-3}$ / unit lnR", lab="initial superdroplets", 
    c='C0', perlnR=True, smooth=False)
        
  ax.legend()
  plt.tight_layout()
  plt.show()
  #############################################