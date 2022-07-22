import numpy as np
import matplotlib.pyplot as plt


def radii_growth_legend(r0, i, nsupers):

  label = None
  if i == 0:
      label = '{:.2g}\u03BCm'.format(r0*1e6)
      label = 'min r0 = '+label
  elif i == nsupers-1:
      label = '{:.2g}\u03BCm'.format(r0*1e6)
      label = 'max r0 = '+label

  return label



def plot_ndrops_individual_growths(t, r, ndrops2plot, is_show=True):
  ''' plots of droplet radii growth '''
  
  fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
  ax1, ax2 = axs

  for i in range(ndrops2plot):
      label = radii_growth_legend(r[0,i], i, ndrops2plot)
      ax1.plot(t, r[:,i]*1e6, linewidth=0.8, label=label)
      ax2.plot(t, r[:,i]/r[0,i], linewidth=0.8, label=label)

  for ax in axs:
    ax.set_xlabel('time /s')
    ax.set_yscale('log')
    ax.legend(fontsize=13)

  ax1.set_ylabel('droplet radius /\u03BCm')
  ax2.set_ylabel('droplet radius / initial radius')

  fig.tight_layout()
  
  if is_show:
    plt.show()

  return fig, axs