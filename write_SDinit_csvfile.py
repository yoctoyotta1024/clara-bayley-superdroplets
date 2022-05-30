import numpy as np
import matplotlib.pyplot as plt
import csv

plt.rcParams.update({'font.size': 14})

##############################################
###  CONSTANTS!! MAKE SURE THESE ARE SAME  ###
###  AS IN constants.hpp FOR ODE SOLVER!!  ###
RGAS_DRY   = 8.314462618/0.0289647                # specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
P0         = 100000                               # pressure [Pa]
TEMP0      = 273.15                               # temperature [K]
RHO0       = P0/(RGAS_DRY*TEMP0)                  # density [Kg/m^3]
R0         = 1e-6                                 # droplet radius lengthscale [m]
N0         = 1e6                                  # droplet multiplicity [m^-3]

Rho_sol    = 2200/RHO0                            # dimensionless solute density []
##############################################


##############################################
### droplet properties for initialisation ### 
INITDROPSCSV = "dimlessSDinit.csv"
nsupers          = 20                              # no. of distinct superdrops (different initial radii (evenly spaced between ln(rspan))

rspan            = [1e-8, 1e-5]                    # initial range of droplet radii [m]
#mus             = [0.075e-6]                      # [m] geometric mean droplet radius
#sigs            = [1.5]                    
#n_as            = [1e9]                           # [m^-3] total no. concentration of droplets          
mus              = [0.02e-6, 0.2e-6, 3.5e-6]               
sigs             = [1.55, 2.3, 2]                    
n_as             = [1e6, 0.3e6, 0.025e6]                  
#############################################


    
#############################################
###      functions for getting droplet    ###
###          radii distribution           ###

def dimless_lnr_dist(rspan, nbins, n_a, mu, sig):
  ''' returns 'nbins' no. of samples from 
  lognormal distribution with even spacing
  in lnr space between ln(r[0]) and ln(rs[1])'''
  
  edgs = np.linspace(np.log(rspan[0]), np.log(rspan[1]), nbins+1)        # edges to lnr bins
  wdths = edgs[1:]- edgs[:-1]                                      # lnr bin widths
  lnr = (edgs[1:]+edgs[:-1])/2                                     # lnr bin centres
  lnnorm = lnnormal_dist(np.e**lnr, n_a, mu, sig)                  # lognormal values

  ### make radii and no. concentrations dimensionless
  lnr = lnr - np.log(R0)                                           
  edgs = edgs - np.log(R0)
  wdths = edgs[1:]- edgs[:-1]                                      # lnr bin widths
  eps = np.asarray([int(n) for n in lnnorm*wdths])/N0     # dimless multiplicities from lognorm dist
  

  return lnr, eps, wdths, edgs



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

#############################################


#############################################
###     functions for plotting droplet    ###
###          radii distribution           ###

def plot_histogram(data, freq, span, nbins, 
          lab=None, c='k', plotdiff=False, ax=None):
  
  if ax==None:
    fig, ax = plt.subplots(figsize=(6,4))
  
  hedgs = np.linspace(np.log(span[0]), np.log(span[1]), nbins+1)        # edges to lnr bins
  hwdths = hedgs[1:]- hedgs[:-1]                                      # lnr bin widths
  hcens = (hedgs[1:]+hedgs[:-1])/2                                     # lnr bin centres
  hist, hedgs = np.histogram(np.log(data), bins=hedgs, weights=freq)
    
  ax.bar(hcens+np.log(1e6), hist, hwdths, label=lab, color=c)
  ax.set_xlabel('ln(r /\u03BCm)')
  ax.set_ylabel('Droplet Conc. [m$^{-3}$]')
  ax.legend()

  return hist, hedgs
  
def linear_twinax(ax, lnr, eps):
    ''' linear x axis for lognormal 
    radius distribution plots'''
     
    ax1 = ax.twiny()
   
    ax1.plot(np.e**lnr*1e6, eps, alpha=0)
    ax1.set_xscale('log')
    ax1.set_xlabel('radius, r /\u03BCm)')
    
    ax.set_xlabel('ln(r /\u03BCm)')
    ax.set_ylabel('No. concentration of droplets [m$^{-3}$]')
  
    return ax, ax1

#############################################
 





#############################################
###       write initial superdroplet      ###
###        conditions to .csv file        ###

for d in range(len(mus)):
    
  mu, sig, n_a = mus[d], sigs[d], n_as[d]
  lnr1, eps1 = dimless_lnr_dist(rspan, nsupers, n_a, mu, sig)[0:2]
  r0 = np.e**(lnr1)
    
  if d==0:
    eps = eps1
  else:                        
    eps = eps+eps1    
m_sol = Rho_sol*4/3*np.pi*r0**3                   # assuming initially dry areosol droplets (dimless mass_solute)


print("Writing inital droplet distribution to: ./init_superdroplets.csv")

with open('./'+INITDROPSCSV, 'w', encoding='UTF8') as f:
  writer = csv.writer(f)
  writer.writerow(["/* Initial Dimensionless Superdroplets Data"])
  writer.writerow([" N0 = "+str(N0)+"[m^-3]", "  R0 = "+str(R0)+"[m]", 
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
plot_histogram(R0*r0, eps*N0, rspan, nsupers, 
              lab='superdroplets', c='C0', ax=ax)          

# plot lognormal curves with x10 more bins for reference
for d in range(len(mus)):
    pltlnr, plteps, pltwdths = dimless_lnr_dist(rspan,                     
                    nsupers*10, n_as[d], mus[d], sigs[d])[0:3]
    ax.plot(pltlnr+np.log(R0*1e6), N0*plteps*10, 
                     color='C'+str(d+1), label='lognormal '+str(d)) 

linear_twinax(ax, np.log(R0*r0), eps)
ax.legend()
plt.tight_layout()
plt.savefig("initial_superdroplet_dist.png", dpi=400, bbox_inches="tight")
print("Initial Distribution saved as: ./initial_superdroplet_dist.png")
plt.show()

#############################################