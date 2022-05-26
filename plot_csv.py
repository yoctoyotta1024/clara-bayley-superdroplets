import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

### Characterstic time, velocity, temperate etc. 
### used to make ODEs dimensionless (see constants.hpp)
W0         = 0.5;                         # characteristic velocity [m/s]
TIME0      = 4000;                        # timescale [s]
P0         = 100000;                      # pressure [Pa]
TEMP0      = 273.15;                      # temperature [K]
#R0         = 1e-6;                        # droplet radius lengthscale [m]
#N0         = 1e6;                         # droplet multiplicity [m^-3]


### Input parameters to ODE (see init.hpp)
w           = 0.5;                         # vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
z_init      = 0;
temp_init   = 273.15;                     # initial parcel temperature [T]
p_init      = 100000;                     # initial pressure [Pa]
relh_init   = 60;                         # initial relative humidity (%)
qc_init     = 0;                          # initial liquid water content []



with open("sundials2_sol.csv") as file_name:
    t, p, temp = np.loadtxt(file_name, delimiter=",", comments="#", unpack=True)

#### RE-Dimensionalise Solution ###
temp0, p0 = temp_init, p_init

t = TIME0 * t
z = t*w + z_init
p, temp = P0*p, TEMP0*temp
qv, qc = 0, 0
#theta = dry_pot_temp(temp, p, qv)
#r = R0*sol.y[4:nsupers+4]
#xyz = sol.y[nsupers+4:nsupers*4+4]
#xyz = np.reshape(xyz, [nsupers,3, len(time)])*(W0*T0)

print("--- Data Shapes ---")
print(temp.shape, p.shape, temp.shape)
#print(temp.shape, p.shape, qv.shape, qc.shape)
#print(r.shape, xyz.shape)
print("\n--- Non Dimensional Max/Mins of Data ---")
print("time:", np.amax(t)/TIME0, np.amin(t)/TIME0)
print("p:", np.amax(p)/P0, np.amin(p)/P0)
print("temp:", np.amax(temp)/TEMP0, np.amin(temp)/TEMP0)
print("qv, qc", (np.amax(qv), np.amax(qc)), (np.amin(qv), np.amin(qc)))
#print("droplet r:", np.amax(r)/R0, np.amin(r)/R0)
#print("droplet z coord:", np.amax(xyz[:,2,:]), np.amin(xyz[:,2,:]))

fig, [axa, axb] = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
ax2a = axa.twinx()
axa.plot(t, temp, label="T", color='r')
ax2a.plot(t, p/100, label="p", color='b')
axa.set_xlabel("Time, t /s")
axa.set_ylabel("Temperature, T /K")
ax2a.set_ylabel("Pressure, P /hPa")
axa.legend(loc="upper left")
ax2a.legend(loc="upper right")

ax2b = axb.twiny()
axb.plot(temp, z/1000, label="T", color='r')
ax2b.plot(p/100, z/1000, label="p", color='b')
axb.set_xlabel("Temperature, T /K")
ax2b.set_xlabel("Pressure, P /hPa")
axb.set_ylabel("Height, z /Km")
axb.legend(loc="lower right")
ax2b.legend(loc="upper right")


plt.tight_layout()
plt.savefig("solution.png", dpi=400, bbox_inches="tight")
plt.show()


