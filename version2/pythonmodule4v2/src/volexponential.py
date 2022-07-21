import numpy as np

#############################################
###      functions for getting droplet    ###
###          radii distribution           ###
###     from randmly sampling volume      ###
###           exponential distribution    ###



def dimless_randomexpo_dist(nsupers, VOL, n_a, r_a, rho=1, R0=1):
  ''' return nsupers superdroplets with eps = n_a*vol/nsupers
  and radius taken from random rample of exponential volume
  (or mass) distribution '''
  
  eps = np.full(nsupers, n_a*VOL/nsupers)

  vol_a = 3/4*np.pi*rho*r_a**3
  vols = np.random.exponential(vol_a, nsupers)
  r = (4*vols/(3*np.pi*rho))**(1/3)
 
  return r/R0, eps