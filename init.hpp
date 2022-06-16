#ifndef INIT
#define INIT


#include "constants.hpp"
namespace dlc = dimless_constants;

/* File containing initial conditions & setup
for ODE solver */

/* initial parcel conditions */
double iVOL        = 1e6;                        // parcel volume [m^3]
double iW          = 0.5;                        // vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
double temp_init   = 273.15;                     // initial parcel temperature [T]
double p_init      = 100000;                     // initial pressure [Pa]
double relh_init   = 60;                         // initial relative humidity (%)
double qc_init     = 0;                          // initial liquid water content []
bool doExpand      = true;                      // enable condensation droplet growth
bool doCond        = true;                      // enable condensation droplet growth
bool doColl        = false;                       // enable collisions of superdroplets

/* integration params */
double tspan[2]    = {0, 4000};                   // time span of integration [s]
int nout           = 40;                          // No. time points to evaluate (save data at)
double rtol        = 1e-6;                        // relative tolerance (tol) for integration
double atols[4]    = {1e-6,1e-6,1e-6,1e-6};       // absolute tols for thermodynamics [P, T, qv, qc]


/* droplet init params. First create superdroplet eps, r0
and m_sol data using python "create_superdrop_init.py" */
//int nsupers        = 32768;                           // no. distinct superdrop objects in array
int nsupers        = 10;                           // no. distinct superdrop objects in array
double iRho_l      = dlc::Rho_l;
double iRho_sol    = dlc::Rho_sol;
double iMr_sol     = dlc::Mr_sol;
int iIONIC         = dlc::IONIC;                 

/* SD Model timestep parameters */                                       
double cond_tstep     = 0.0001;                    // time between SD condensation events  = ceil(coll/cond)*min(coll,cond) [s]
double coll_tstep     = 1;                    // time between SD collision events = ceil(coll/cond)*min(coll,cond) [s]
       





#endif // INIT