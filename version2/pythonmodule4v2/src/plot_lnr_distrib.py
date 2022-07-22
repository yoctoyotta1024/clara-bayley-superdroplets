import numpy as np


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


def gaussian_kernel_smoothing(hist, hcens, sig):

    for h in range(len(hist)):
        kernel = 1/(np.sqrt(2*np.pi)*sig)*np.exp(-(hcens - hcens[h])**2/(2*sig**2))
        kernel = kernel/np.sum(kernel)
        hist[h] = np.sum(hist*kernel)  

    hist = np.where(hist<1e-16, 0, hist)  

    return hist, hcens



###################################################################
 







