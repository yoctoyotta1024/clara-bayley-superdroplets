import numpy as np


######## functions for reading .csv files into python
### and converting into variables with dimensions ########

def get_soldata(sol_filename, TIME0, P0, TEMP0):
    
    #### Load data from .csv file ###
    with open(sol_filename) as file_name:
        t, p, temp, qv, qc = np.loadtxt(file_name, delimiter=",", comments="/*", unpack=True)

    print("--- Raw Data Shapes ---")
    print("variables: t, p, temp, qv, qc")
    print(t.shape, p.shape, temp.shape, qv.shape, qc.shape)

    print("--- Non Dimensional Max/Mins of Data ---")
    print("time:", np.amin(t), np.amax(t))
    print("p:", np.amin(p), np.amax(p))
    print("temp:", np.amin(temp), np.max(temp))
    print("(qv, qc)", (np.amin(qv), np.amin(qc)), (np.amax(qv), np.amax(qc)), "\n")
    
    
    return t*TIME0, p*P0, temp*TEMP0, qv, qc      ### (Re-)dimensionalise and return


def get_SDdata(SDsol_filename, nsupers, R0, RHO0):

    ### Load data from SD .csv file ###
    with open(SDsol_filename) as file_name:
        drops = np.loadtxt(file_name, delimiter=",", comments="/*", unpack=False)
    eps = drops[:,0:nsupers]
    r = drops[:,nsupers:2*nsupers]
    m_sol = drops[:,2*nsupers:]

    print("--- Raw SD Data Shapes ---")
    print("variables: eps, r, m_sol")
    print(eps.shape, r.shape, m_sol.shape)

    print("--- Non Dimensional Max/Mins of Data ---")
    print("droplet eps:", np.amin(eps), np.amax(eps))
    print("droplet r:", np.amin(r), np.amax(r))
    print("droplet m_sol:", np.amin(m_sol), np.amax(m_sol), "\n")
    

    return eps, r*R0, m_sol*RHO0*R0**3           ### (Re-)dimensionalise and return


