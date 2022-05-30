#ifndef DIFFERENTIALS
#define DIFFERENTIALS

#include <cmath>
#include "constants.hpp"

namespace DC = dimmed_constants;
namespace dlc = dimless_constants;

/* Some Definitions */
/* User-defined vector accessor macro: Ith
  Ith(v,i) references the ith component of the vector v, where i is in
  the range [1..NEQ] (and NEQ is No. of equations ie. variables). 
*/
#define Ith(v,i)    NV_Ith_S(v,i-1)         // i-th vector component i=1..NEQ
#define ZERO  RCONST(0.0)


/* -------------------------------
 * Functions called by the solver
 *------------------------------- */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

//static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

// ---------------------------------------------------------------------- //



/* ---- constants for use in dp_dt function ----- */
realtype dp_dt_const = dlc::W0*dlc::TIME0*DC::G/(DC::RGAS_DRY*dlc::TEMP0); 
realtype lpsrate = 0.0062/dlc::TEMP0*dlc::W0*dlc::TIME0;
realtype tempg = 273.15/dlc::TEMP0;
realtype pg = 100000/dlc::P0;
realtype zg = 0/(dlc::W0*dlc::TIME0);
realtype gamma = DC::G/(DC::RGAS_DRY*0.0062)-1;
/* ----------------------------------------------- */



/* user data structure for passing 
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct {
  realtype w;
  Superdrop drop1; //, drop2, drop3;
} *UserData;


static void InitUserData(UserData data, realtype w, Superdrop testdrop)
{
  data->w = w;
  data->drop1 = testdrop;
  // data->drop2 = 2;
  // data->drop3 = 3;
}

/* ----------------------------------------------- */




/*
 * dp/dt differential equation 
 describing pressure evolution over time.
 */

static realtype dp_dt(realtype t, N_Vector ydot, realtype w)
{
  realtype dp, profile;

  profile = 1 - lpsrate/tempg*(w*t-zg);
  profile = pow(profile, gamma);

  dp = Ith(ydot,1) = -dp_dt_const*pg/tempg*profile;
  // dp = -dp_dt_const*pg/tempg*profile;

  return(dp);
}




/*
 * dtemp/dt differential equation describing 
temperature evolution solely due to pressure 
changes in parcel during for adiabatic 
process (no heat loss). Parcel has water vapour
mass mixing ratio (m_v/m_dry) = qv and liquid 
water mass mixing ratio (m_c/m_dry) = qc.
(Assumes instantaneous volume change 
of parcel to change pressure according
to dp_dt_profile.) 
True dTemp/dt = dtemp * Temp0/t0 
 */

static realtype dtemp_dt_adia(N_Vector ydot, realtype p, 
                  realtype temp, realtype qv, realtype qc)
{
  realtype dtemp, rho_d, cp_moist;
  
  rho_d = dlc::Mr_ratio/(dlc::Mr_ratio+qv)*p/temp;      // density of dry parcel (p_dry/temp)    
  
  cp_moist = dlc::Cp_dry + dlc::Cp_v*qv + dlc::C_l*qc;       // specific heat capacity of moist parcel 

  dtemp = Ith(ydot,2) = dlc::Rgas_dry/(rho_d*cp_moist) * Ith(ydot,1);
  // dtemp = dlc::Rgas_dry/(rho_d*cp_moist) * Ith(ydot,1);

  return(dtemp);
}






// /*
//  * dr/dt differential equation describing radius
//   evolution over time for each superdrop.
//  */

// static realtype dr_dt(realtype t, N_Vector ydot, realtype w)
// {
//   realtype dp, profile;

//   profile = 1 - lpsrate/tempg*(w*t-zg);
//   profile = pow(profile, gamma);

//   dp = Ith(ydot,1) = -dp_dt_const*pg/tempg*profile;
//   // dp = -dp_dt_const*pg/tempg*profile;

//   return(dp);
// }






/*
 * Simple function f(t,y) called by ODE solver to 
 solve differential equations over time.
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  
  UserData data;
  data = (UserData) user_data;
  
  realtype p, temp, qv, qc, w;
  p = Ith(y,1); 
  temp = Ith(y,2);
  qv = 0.0;  //Ith(y,3);
  qc = 0.0;  //Ith(y,4);
  w = data -> w;

 // realtype w = 0.5/dlc::W0;
  dp_dt(t, ydot, w);     

  dtemp_dt_adia(ydot, p, temp, qv, qc);    

  


  return(0);
}










// /*
//  * g routine. Compute functions g_i(t,y) for i = 0,1.
//  */

// static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
// {
//   realtype y1, y3;

//   y1 = Ith(y,1); y3 = Ith(y,3);
//   gout[0] = y1 - RCONST(0.0001);
//   gout[1] = y3 - RCONST(0.01);

//   return(0);
// }












#endif //DIFFERENTIALS
