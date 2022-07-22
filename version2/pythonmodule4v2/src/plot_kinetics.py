import numpy as np
import matplotlib.pyplot as plt

import src.generic_axfuncs as axfuncs



def plot_kinetics_against_time(t, p, temp, qv, qc, 
      relh, s, dry_adia, is_show=True):

  fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,6))
  axs = axs.flatten()

  axfuncs.axplt(axs[0], t, p/100, ylab='Pressure /mbar', c=0)
  
  axfuncs.axplt(axs[1], t, dry_adia, lab='dry adiabat', c=2, l='--')
  axfuncs.axplt(axs[1], t, temp, ylab='T /K', lab='temp', c=3)

  axfuncs.axplt(axs[2], t, qv+qc, lab='total', c=7, l='--')
  axfuncs.axplt(axs[2], t, qv, lab='vapour', c=5)
  axfuncs.axplt(axs[2], t, qc, xlab="time /s", 
      ylab='Water Content', lab='liquid', c=6)

  axfuncs.axplt(axs[3], t, relh, 
      xlab="time /s", ylab='relative humidity', c=4)
  axfuncs.axplt(axs[3].twinx(), t, s, 
      ylab="supersaturation", c=4)

  for ax in [axs[1], axs[2]]:
    ax.legend(fontsize=11, loc='upper right')
  
  axs[0].set_yscale("log")
  
  fig.tight_layout()
  
  if is_show:
    plt.show()

  return fig, axs





def plot_kinetics_against_pressure(z, p, temp, qv, qc, 
      theta, pv, psat, relh, s, dry_adia, dry_adia_theta, is_show=True):

  fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15,6))
  axs = axs.flatten()

  axfuncs.axplt(axs[0],  z/1000, p/100, 'z=w*t /km', 'P /mbar', c=0)
  axfuncs.axplt(axs[1], dry_adia, p/100, lab='dry adiabat', c=2, l='--')
  axfuncs.axplt(axs[1], temp, p/100, 'T /K', lab='temp', c=3)
  axfuncs.axplt(axs[2], dry_adia_theta, p/100, lab='dry adiabat', c=2, l='--')
  axfuncs.axplt(axs[2], theta, p/100, '\u03B8 /K', lab='\u03B8', c=1)

  axfuncs.axplt(axs[4], qv+qc, p/100, lab='total', c=7, l='--')
  axfuncs.axplt(axs[4], qc, p/100, lab='liquid', c=6)
  axfuncs.axplt(axs[4], qv, p/100, 'Water Content', lab='vapour', c=5)

  axfuncs.axplt(axs[3], psat/100, p/100, lab='Saturation P', c=2, l='--')
  axfuncs.axplt(axs[3], pv/100, p/100, 
          'Water Vapour Pressure /mbar','Vapour P', c=4)

  axfuncs.axplt(axs[5], s, p/100, 'Supersaturation' , c=4)
  axfuncs.axplt(axs[5].twiny(), relh, p/100, 'Relative Humidity', c=4)

  for ax in axs:
    ax.set_yscale('log')
    ax.invert_yaxis()
  
  for ax in [axs[2], axs[3], axs[4]]:
      ax.legend()

  fig.tight_layout()

  if is_show:
    plt.show()

  return fig, axs





def plot_moist_static_energy(t, mse, is_show=True):

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4))

  axfuncs.axplt(ax, t, mse, xlab ="time /s", 
    ylab='moist static energy / J/Kgm$^{-3}$', c="k")
  axb = ax.twinx()
  axfuncs.axplt(axb, t, ((mse-mse[0])/mse[0]*100), 
    ylab='% change in moist static energy', c="k")

  fig.tight_layout()

  if is_show:
    plt.show()

  return fig, ax