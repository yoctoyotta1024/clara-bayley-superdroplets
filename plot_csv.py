import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

### Characterstic time, velocity, temperate etc. 
### used to make ODEs dimensionless (see constants.hpp)
W0         = 0.5;                         # characteristic velocity [m/s]
TIME0      = 4000;                        # timescale [s]
P0         = 100000;                      # pressure [Pa]
TEMP0      = 273.15;                      # temperature [K]
R0         = 1e-6;                       # droplet radius lengthscale [m]
N0         = 1e6;                        # droplet multiplicity [m^-3]
Mr_water     = 0.01801528                 # molecular mass of water
Mr_dry       = 0.0289647                  # molecular mass (per 1 mole) of dry air (M_r)       
mr_ratio     = Mr_water/Mr_dry
nsupers    = 20

### Input parameters to ODE (see init.hpp)
w           = 0.5;                        # vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
z_init      = 0;
temp_init   = 273.15;                     # initial parcel temperature [T]
p_init      = 100000;                     # initial pressure [Pa]
relh_init   = 60;                         # initial relative humidity (%)
qc_init     = 0;                          # initial liquid water content []




def saturation_pressure(T):
    ''' Calculate the equilibrium vapor pressure 
    of water over liquid water ie. the
    saturation pressure (psat). Equation taken from
    typhon.physics.thermodynamics.e_eq_water_mk.
    Real temp = T*Temp0, dimensionless psat = psat/P0'''

    T = T*TEMP0                               # real T [K]
    
    if np.any(T <= 0):
        err = 'T must be larger than 0K.'
        raise ValueError(err)

    lnpsat = (54.842763                    # ln(psat) [Pa]
         - 6763.22 / T
         - 4.21 * np.log(T)
         + 0.000367 * T
         + np.tanh(0.0415 * (T - 218.8))
         * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025 * T))

    return np.exp(lnpsat)/P0               # dimensionless psat




#### Load data from .csv files ###
with open("sundials2_sol.csv") as file_name:
    t, p, temp, qv, qc = np.loadtxt(file_name, delimiter=",", comments="/*", unpack=True)

with open("sundials2_SDsol.csv") as file_name:
    r = np.loadtxt(file_name, delimiter=",", comments="/*", unpack=False)




#### RE-Dimensionalise Solution ###

temp0, p0, qc0 = temp_init, p_init, qc_init
pv0 = relh_init/100 * saturation_pressure(temp0)
qv0 = mr_ratio * pv0/(p-pv0)

t = TIME0 * t
z = t*w + z_init
p, temp = P0*p, TEMP0*temp
r = R0*r.T
#theta = dry_pot_temp(temp, p, qv)
#xyz = sol.y[nsupers+4:nsupers*4+4]
#xyz = np.reshape(xyz, [nsupers,3, len(time)])*(W0*T0)

print("--- Data Shapes ---")
print(t.shape, temp.shape, p.shape, qv.shape, qc.shape)
print(r.shape)
print("\n--- Non Dimensional Max/Mins of Data ---")
print("time:", np.amax(t)/TIME0, np.amin(t)/TIME0)
print("p:", np.amax(p)/P0, np.amin(p)/P0)
print("temp:", np.amax(temp)/TEMP0, np.amin(temp)/TEMP0)
print("qv, qc", (np.amax(qv), np.amax(qc)), (np.amin(qv), np.amin(qc)))
print("droplet r:", np.amax(r)/R0, np.amin(r)/R0)
#print("droplet z coord:", np.amax(xyz[:,2,:]), np.amin(xyz[:,2,:]))



### plots of p, temp, qv and qc evolution
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 8))
[ax1, ax2, ax3, ax4] = [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]
ax1b = ax1.twinx()
ax2b = ax2.twiny()
ax1.plot(t, temp, label="T", color='r')
ax1b.plot(t, p/100, label="P", color='b')
ax2.plot(temp, z/1000, label="T", color='r')
ax2b.plot(p/100, z/1000, label="P", color='b')
ax1.set_xlabel("Time, t /s")
ax1.set_ylabel("Temperature, T /K")
ax1b.set_ylabel("Pressure, P /hPa")
ax2.set_xlabel("Temperature, T /K")
ax2b.set_xlabel("Pressure, P /hPa")
ax2.set_ylabel("Height, z /Km")
ax1.legend(loc="upper left")
ax1b.legend(loc="upper right")
ax2.legend(loc="lower left")
ax2b.legend(loc="upper left")

ax3.plot(t, qv, label="q$_{v}$", color='purple')
ax3.plot(t, qc, label="q$_{c}$", color='brown')
ax4.plot(qv, z/1000, label="q$_{v}$", color='purple')
ax4.plot(qc, z/1000, label="q$_{c}$", color='brown')
ax3.set_xlabel("Time, t /s")
ax3.set_ylabel("Water Content")
ax4.set_xlabel("Water Content")
ax4.set_ylabel("Height, z /Km")
ax3.legend(loc="lower right")
ax4.legend(loc="lower right")

plt.tight_layout()
plt.savefig("solution.png", dpi=400, bbox_inches="tight")
plt.show()




### plots of droplet radii growth
fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
labs = ['{:.2g}\u03BCm'.format(r*1e6) for r in r[:,0]]
for i in range(nsupers):
    label = None
    if i == 0:
        label = 'min r0 = '+labs[i]
    elif i == nsupers-1:
        label = 'max r0 = '+labs[i]
    ax1.plot(t, r[i]*1e6, linewidth=0.8, label=label)
ax1.set_xlabel('time /s')
ax1.set_ylabel('droplet radius /\u03BCm')
ax1.set_yscale('log')
ax1.legend(fontsize=13)

for i in range(nsupers):    
    label = None
    if i == 0:
        label = 'min r0 = '+labs[i]
    elif i == len(r)-1:
        label = 'max2 r0 = '+labs[i]
    ax2.plot(t, r[i]/r[i,0], linewidth=0.8, label=label) #,color=cols[i])
ax2.set_xlabel('time /s')
ax2.set_ylabel('droplet radius / initial radius')
ax2.set_yscale('log')
ax2.legend(fontsize=13)

plt.tight_layout()
plt.savefig("SDsolution.png", dpi=400, bbox_inches="tight")
plt.show()