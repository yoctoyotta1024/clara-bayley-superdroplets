import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

from constants import *
terminal_const = -2*g/(9*dynvisc)

#plt.rcParams.update({'font.size': 12})


class Superdroplets():
    '''Class for all superdroplets, 
    eg. for object drops = Superdroplets(r0) where r0
    is array of drop radii, drops is object with
    radii. multiplicty and other attributes (eg. mass)
    of all superdroplets'''
    
   
    def __init__(self, nsupers, r, eps, m_sol, mr_sol, ionic, rho_sol, rho_l):
     
        self.nsupers = nsupers
        if self.nsupers != len(r):
            print('warning nsupers != number of different superdroplet radii')
            print(nsupers, len(r))
    
        # Common attributes of all superdroplets
        self.rho_l = rho_l                                # density of liquid in droplets (=density of water at 300K) [Kg/m^3]
       
        # droplet solute properties
        self.mr_sol = mr_sol                             # Mr of solute [g/mol]
        self.ionic = ionic                               # degree ionic dissociation (van't Hoff factor)
        self.rho_sol = rho_sol                           # density of dry droplets  (solute density)
        
        # attributes unique to each sperdroplet object
        self.eps = eps                                   # array of multiplicity of droplets
        self.r = r                                       # array of droplet radii
        self.m_sol = m_sol                               # mass of solute dissovled (for each droplet in superdrop) [kg]
        self.b = 43e-6*self.m_sol*self.ionic/self.mr_sol      # kohler b factor: 43e-6 m^3/mol * mass_solute/ Mr_solute [eq.6.22]
        self.dry_r = (3*self.m_sol/(4*np.pi*self.rho_sol))**(1/3)   # dry radius of droplets
        self.addsol = (self.rho_sol - self.rho_l)*self.dry_r**3    # contribution to density due to solute = addsol/self.r^3
        
        self.r0 = r                                      # array of droplet initial radii
        v_w0 = 4/3*np.pi*(self.r0**3 - self.dry_r**3)
        self.m_w0 = self.rho_l*v_w0                      # initial droplet water mass
        self.m0 = self.m_w0 + self.m_sol                 # initial droplet masses
       

    @property
    def m(self):
        ''' total mass of droplet
            (water + dry areosol) '''
        
        rho_eff = self.rho_l + self.addsol/(self.r**3)
             
        return  rho_eff*4/3*np.pi*(self.r)**3
   

    @property
    def m_w(self):
        ''' mass of water in droplet '''
        
        v_w = 4/3*np.pi*(self.r**3 - self.dry_r**3)
        
        return self.rho_l*v_w                           # droplet water mass
        
        
    @property
    def vol(self):
        ''' droplet volume (inc. solute) '''
        
        return 4/3*np.pi*(self.r)**3                     # droplet volume
    
    
    
    
    def velocity(self):
        ''' returns flattened array of 
        velocities of each particle 
        vflat = [vx1, vy1, vz1, 
        vx2, vy2, vz2, vx3, vy3 ...]'''
        
        #rho_eff = self.rho_l + self.addsol/(self.r**3)
        #terminal_v = terminal_const*rho_eff*(self.r)**2
    
        rho_rsqrd = self.rho_l*self.r**2 + self.addsol/self.r                  # r^2 * rho_effective (inc. solute mass)
        terminal_v = terminal_const*rho_rsqrd
        
        #v = np.zeros((self.nsupers, 3))
        #v[:,0] = dx
        #v[:,1] = dy
        #v[:,2] = terminal_v           # z velocity of each particle = terminal velocity
        #return v.flatten()
        
        vflat = np.zeros(self.nsupers*3)
        vflat[2::3] = terminal_v           # z velocity of each particle = terminal velocity
        #vflat[::3] = dx
        #vflat[1:3] = dy 
        
        return vflat
    
       
        
    
    
    def init_positions(self, pclsiz):
        ''' initialise droplet positions 
        in relative coordinates (real
        x of particle = x of parcel centre + x).
        pclsize is [Lx/2,Ly/2,Lz/2] dimensions of parcel'''
        
        dimless = np.random.rand(self.nsupers,3)-0.5                 # dimensionless random (x,y,x) of each particle (-1<d<1)
        
        self.coords0 = np.asarray(pclsiz)[None,:]*dimless            # all superdroplet droplet positions [droplet, [x,y,z]]
        
        return self.coords0
    
    

    def kohler_factors(self, temp):
        ''' calculate b in kelvin factor (1-b/r^3)
        and a in raoult factor (exp^(a/r)) to
        account for curvature and effect of solute
        on radial growth of droplet respectively.
        Using eq.6.24 and eq.6.22 of lohmann, luond
        and mahrt intro 2 clouds textbook'''
        
        a = 3.3e-7/temp                                #2*surface_tens (=0.0756N/m)/(rho_l*rgas_v*temp) [eq.6.24]
            
        return a, self.b
        
        
    def active(self, s_ratio, temp):
        ''' calculate b in kelvin factor (1-b/r^3)
        and a in raoult factor (exp^(a/r)) to
        account for curvature and effect of solute
        on radial growth of droplet respectively.
        Using eq.6.24 and eq.6.22 of lohmann, luond
        and mahrt intro 2 clouds textbook'''
        
        
        a, b = self.kohler_factors(temp)
       
        r_act = np.sqrt(3*b/a)
        active = np.where(self.r >= r_act, True, False)
        
        return active
    
    
    def distribution_stats(self, r=None, geo=True):
        '''calculate mean, standard deviation,
        skewness and kurtsis (tailedness) of
        distribution. 
        (geometric by default) 
        NOTE! all metrics (mean, deviaiton etc.) 
        all divided by n not n-1.'''
        
        if type(r) == type(None): r = self.r
        ntot = np.sum(self.eps)
        
        if r.ndim > 1:
            sl = np.s_[None,:]
            eps = np.asarray(self.eps)[np.s_[:, None]]
        else:
            sl = np.s_[None]
            eps = np.asarray(self.eps)
        
        if geo:
            mean = sts.mstats.gmean(r, axis=0, weights=self.eps)
            sqrds = np.sum(eps*((np.log(r/mean[sl]))**2), axis=0)
            lnstd = np.sqrt(sqrds/ntot)
            std = np.exp(lnstd)                                         # first standard deviation (geometric)
            
            m3 = np.sum(eps*((np.log(r/mean[sl]))**3), axis=0)/ntot     # 3rd and 4th moments (geometric)
            m4 = np.sum(eps*(r-mean[sl])**4, axis=0)/ntot
            skew = m3/np.exp(lnstd**3)                                  # 3rd and 4th standardised (geo)moments
            kurt = m4/np.exp(lnstd**4)

        else:
            mean = np.average(r, axis=0, weights=self.eps)     # mean of distribution
            sqrds = np.sum(eps*(r-mean[sl])**2, axis=0)       
            std = np.sqrt(sqrds/ntot)                           # std dev of distribution
            
            m3 = np.sum(eps*(r-mean[sl])**3, axis=0)/ntot        # 3rd and 4th arithmetic moments of disribution
            m4 = np.sum(eps*(r-mean[sl])**4, axis=0)/ntot
            skew = m3/std**3                                     #3rd and 4th standardised moments of distribition
            kurt = m4/std**4

        
        return mean, std, skew, kurt
    
    
    
    
    
    
### functions for plotting droplet distributions

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


def get_lnr_dist(rs, nbins, n_a, mu, sig):
    ''' returns 'nbins' no. of samples from 
    lognormal distribution with even spacing
    in lnr space between ln(r[0]) and ln(rs[1])'''
    
    edgs = np.linspace(np.log(rs[0]), np.log(rs[1]), nbins+1)          # edges to lnr bins
    wdths = edgs[1:]- edgs[:-1]                                      # lnr bin widths
    lnr = (edgs[1:]+edgs[:-1])/2                                     # lnr bin centres
    lnnorm = lnnormal_dist(np.e**lnr, n_a, mu, sig)                  # lognormal values
    
    return lnr, lnnorm, wdths, edgs


def raxes4dists(ax, lnr, eps, edgs):
    ''' axes labels, ticks etc. for lognormal
    radius distribution plots'''
     
    ax1 = ax.twiny()
    linear_wdths = np.e**edgs[1:] - np.e**edgs[:-1]
    ax1.bar(np.e**lnr*1e6, eps, linear_wdths*1e6, alpha=0)
    ax1.set_xscale('log')
    ax1.set_xlabel('radius, r /\u03BCm)')
    
    ax.set_xlabel('ln(r /\u03BCm)')
    ax.set_ylabel('No. concentration of droplets [m$^{-3}$]')
  
    return ax, ax1



def plot_lognormal(eps, lnr, lnnorm, wdths, edgs, pltref=True):
    ''' plots lognormal distribution * bin widths = 
    no. concentration of particles with given radius '''
        
    fig, ax = plt.subplots()
    
    ax.bar(lnr+np.log(1e6), eps, wdths, label='samples')
    
    raxes4dists(ax, lnr, eps, edgs)
    
    if pltref:
        # plot x100 more bins lognormal for reference
        pltlnr, pltnorm, pltwdths = get_lnr_dist(rspan, 
                                nbins*100, n_a, mu, sig)[0:3]
        ax.plot(pltlnr+np.log(1e6), pltnorm*pltwdths*100, 
                    color='k', label='lognormal distribution')
    
    ax.legend()
    plt.tight_layout()
    
    plt.show()

    
    
def radii2dist(data, eps, rs, nbins):
    ''' returns histogram of droplet radii (r) in 
    interval rs[0] t rs[1] in 'nbins' with even 
    spacing in lnr space'''
    
    edgs = np.linspace(np.log(rs[0]), np.log(rs[1]), nbins+1)          # edges to lnr bins
    wdths = edgs[1:]- edgs[:-1]                                      # lnr bin widths
    lnr = (edgs[1:]+edgs[:-1])/2                                     # lnr bin centres
    
    hist, hedgs = np.histogram(np.log(data), bins=edgs, weights=eps)
    
    return hist, wdths, lnr, hedgs



def plot_histogram(ax, data, eps, rspan, nbins, lab=None, c='k', plotdiff=False):

    hist, hwdths, hcens, hedgs = radii2dist(data, eps, rspan, nbins)
        
    ax.bar(hcens+np.log(1e6), hist, hwdths, label=lab, color=c)
    
    ax.set_xlabel('ln(r /\u03BCm)')
    ax.set_ylabel('Droplet Conc. [m$^{-3}$]')
    #ax, ax1 = raxes4dists(ax, hcens, eps, hedgs)
    ax.legend()

    return hist, hwdths, hcens, hedgs
    
    
