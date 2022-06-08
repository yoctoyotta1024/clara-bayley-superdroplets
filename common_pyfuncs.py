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
 
