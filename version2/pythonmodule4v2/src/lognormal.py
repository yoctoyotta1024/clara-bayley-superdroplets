import numpy as np

#############################################
###      functions for getting droplet    ###
###          radii distribution           ###
###     from lognormal distribution(s)    ###

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




#############################################
