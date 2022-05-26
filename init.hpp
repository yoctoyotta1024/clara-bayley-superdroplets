#ifndef INIT
#define INIT

/* File containing initial conditions & setup
for ODE solver */

/* initial parcel conditions */
double w           = 0.5;                        // vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
double z_init      = 0;                          // initial parcel z coordinate [m]
double temp_init   = 273.15;                     // initial parcel temperature [T]
double p_init      = 100000;                     // initial pressure [Pa]
double relh_init   = 60;                         // initial relative humidity (%)
double qc_init     = 0;                          // initial liquid water content []
bool doCond        = false;                      // enable condensation droplet growth

/* integration params */
double tspan[2]    = {0, 4000};                   // time span of integration [s]
int nout           = 1000;                         // No. time points to evaluate (save data at)
double rtol        = 1e-6;                         // relative tolerance (tol) for integration
double atols[2]    = {1e-6, 1e-6};                // absolute tols for [parcel thermodynamics, droplet radii]

                           


                                                






#endif // INIT