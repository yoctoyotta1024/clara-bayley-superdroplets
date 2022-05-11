#!/usr/bin/env python
# coding: utf-8

# In[1]:


import  matplotlib.pyplot as plt
import numpy as np
import typhon as ty
from scipy.constants import g, R
from scipy.integrate import solve_ivp


# In[2]:


### some constants
dry_mr = 0.0289647         # molecular mass (per 1 mole) of dry air (M_r)                 
water_mr = 0.01801528      # molecular mass of water
rgas_dry = R/dry_mr        # specific gas constant for dry air [J/Kg/K]      <-- used in ideal gas equation for hydrostsic rather than moist??
rgas_v = R/water_mr        # specific gas constant for water
cp_dry = 1003.5            # specific heat capacity dry air [J/Kg/K]
latent_v = 2264.705        # specific latent heat of vapourisation of water [J/Kg]


# In[3]:


tempg = 297.2           # temperature at ground level [K]
lps_rate =  0.0062      # lapse rate above zb [K/km]
pb = 95000              # pressure at cloud base [Pa]
zb = 500                # height of cloud base [m]


# In[4]:


z = np.arange(0,2020,20)      # z coordinates of domain


# In[5]:


###  initial T,p and theta profiles

def temp_pressure_profile(z):
    
    if z<=zb:
        temp = tempg
        p = pb * np.exp( -g/(rgas_dry*tempg) * (z-zb))
    else:
        temp = tempg - lps_rate*(z-zb)
        p = pb * (1 - lps_rate/tempg*(z-zb))**(g/(rgas_dry*lps_rate))
        
    return temp, p
    


# In[6]:


#temp = tempg - lps_rate*(z-zb)
#temp = np.where(z<=zb, tempg, temp)
#p = pb * (1 - lps_rate/tempg*(z-zb))**(g/(rgas_dry*lps_rate))
#p_below = pb * np.exp( -g/(rgas_dry*temp0) * (z-zb))
#p=np.where(z<=zb, p_below, p_above)

temp_list, p_list = [], []
for x in z:
    temp, p = temp_pressure_profile(x)
    temp_list.append(temp)
    p_list.append(p)
temp = np.asarray(temp_list)
p = np.asarray(p_list)

theta = temp*(p[0]/p)**(rgas_dry/cp_dry)


# In[7]:


fig, axs = plt.subplots(nrows=2, ncols=2)

axs[0,0].plot(temp, z/1000)
axs[0,0].set_xlabel('T /K')
axs[0,0].set_ylabel('z /km')

axs[1,0].plot(p/100, z/1000)
axs[1,0].set_xlabel('pressure /mbar')
axs[1,0].set_ylabel('z /km')

axs[0,1].plot(theta, z/1000)
axs[0,1].set_xlabel('\u03B8 /K')
axs[0,1].set_ylabel('z /km')

axs[1,1].plot(theta, p/100)
axs[1,1].set_xlabel('\u03B8 /K')
axs[1,1].set_ylabel('pressure /mbar')

axs[1,1].set_yscale('log')
axs[1,1].invert_yaxis()

plt.tight_layout()


# In[8]:


###  initial saturation pressure profile
p_sat = ty.physics.e_eq_water_mk(temp)

fig, ax = plt.subplots()
ax.plot(p_sat/100, z/1000)
ax.set_xlabel('initial saturation pressure /mbar')
ax.set_ylabel('z / km')
plt.show()


# In[9]:


### initial parcel
w = 2                 # vertical wind speed of parcel [m/s]
pcl0_z = z[0]         # initial z coordinate of parcel
pcl0_temp = temp[0]   # initial temp of parcel
pcl0_p = p[0]         # initial pressure of parcel
pcl0_relh = 50        # initial relative humidity of parcel (%)
pcl0_pc = 0           # initial liquid water content in parcel

pcl0_pv = pcl0_relh/100 * ty.physics.e_eq_water_mk(pcl0_temp)    # initial vapour pressure of water vapour in parcel
pcl_qv = water_mr/dry_mr * pcl0_pv/(pcl0_p-pcl0_pv)              # water mass mixing ratio (mass water/mass dry air)


# In[17]:


#t_eval = np.arange(0,1001,1)
t_span = [0,1000]
y0 = [pcl0_z, pcl0_temp, pcl0_pv, pcl0_pc]



def condensation_adjustment(z, delt, temp, p, pcl_qv, pcl_qc):
    ### if saturation pressure < water vapour pressure 
    ### instigate condensation process until p_sat = p_v
    
    p_sat = ty.physics.e_eq_water_mk(temp)
    
    if pcl_qv > p_sat:

        #[calculate adjustment]
        
        q_sat = p_sat/(p-p_sat)
        dqv = (pcl_qv-p_sat)/()
        
    else:
        
        dtemp, dqv, dqc = 0,0,0       # saturation not reached, no condensation adjustment required
        
    return dtemp, dqv, -dqv       
    
    
   

def fun(t, y, w, tpop=[]):
    
    z, pcl_temp, pcl_qv, pcl_qc = y
    
    temp, p = temp_pressure_profile(z) 
    pcl_temp, pcl_p = temp, p         #instanataneous heat and volume change to match surroundings
    
    ### thermodynamics of parcel
    if t == 0:                        # in  first instance skip evolution of dT, dq_v and dq_c 
        dtemp, dqv, dqc = 0,0,0 
    else:
        t_m1 = tpop.pop()
        delt = t-t_m1
        dtemp, dqv, dqc = condensation_adjustment(z, delt, temp, pcl_qv, pcl_qc)
    tpop.append(t)

    ### dynamics of parcel
    dy = w
    #tau = 200
    #dy = w*np.sin(2*np.pi*t/tau)
    
    return [dy, dtemp, dqv, -dqv]      # dqc = -dqv  (vapour<-->condensation)
    

sol = solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, args=[w])      # max_step not specified
time = sol.t
z = sol.y[0]

plt.plot(time,z)


# In[ ]:





# In[ ]:





# In[ ]:




