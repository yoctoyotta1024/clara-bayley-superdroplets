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
#define SDloop(i,nsupers) for(int i=0; i<nsupers; i++)  //for loop over all superdroplets
#define ZERO  RCONST(0.0)


/* -------------------------------
 * Functions called by the solver
 *------------------------------- */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

//static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

// ---------------------------------------------------------------------- //



/* ---- constants for use in dp_dt function ----- */
const realtype dp_dt_const = dlc::W0*dlc::TIME0*DC::G/(DC::RGAS_DRY*dlc::TEMP0); 
const realtype lpsrate = 0.0062/dlc::TEMP0*dlc::W0*dlc::TIME0;
const realtype tempg = 273.15/dlc::TEMP0;
const realtype pg = 100000/dlc::P0;
const realtype zg = 0/(dlc::W0*dlc::TIME0);
const realtype gamma = DC::G/(DC::RGAS_DRY*0.0062)-1;
/* ----------------------------------------------- */

// /* ---- constants for use in diffusion growth ----- */
//const realtype dm_dt_const = 4*M_PI*dlc::Rho_l*dlc::EPS0*pow(dlc::R0, 3.0);
const realtype dm_dt_const = 4*M_PI*dlc::Rho_l*pow(dlc::R0, 3.0)/iVOL;
// /* ----------------------------------------------- */



/* user data structure for passing 
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct {
  realtype w;
  bool doCond;
  int nsupers;
  Superdrop* ptr;
} *UserData;


static void InitUserData(UserData data, realtype w,
 bool doCond, int nsupers, Superdrop* ptr)
{
  data->w = w;
  data->doCond= doCond;
  data->nsupers = nsupers;
  data->ptr = ptr;
}

/* ----------------------------------------------- */


/* -------------- helper functions --------------- */
/*
* Calculate mass mixing ratio
qv = m_v/m_dry = rho_v/rho_dry
given vapour pressure pv = p_v/p_tot.'''
*/
realtype pv2qv(realtype pv, realtype p)
{  
  return dlc::Mr_ratio * pv/(p-pv);
}


/*
* Calculate the equilibrium vapor pressure 
of water over liquid water ie. the
saturation pressure (psat). Equation taken from
typhon.physics.thermodynamics.e_eq_water_mk.
Real temp = T*Temp0, dimensionless psat = psat/P0
 */

realtype saturation_pressure(realtype temp)
{
  realtype T, lnpsat;
  T = temp*dlc::TEMP0;                               // real T [K]
  
  if(T <= 0)
  {
      cout << "psat ERROR: T must be larger than 0K. T = " <<endl;
  }

  lnpsat = (54.842763                               // ln(psat) [Pa]
        - 6763.22 / T
        - 4.21 * log(T)
        + 0.000367 * T
        + tanh(0.0415 * (T - 218.8))
        * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T));

  return exp(lnpsat)/dlc::P0;               // dimensionless psat
}


realtype cp_moist(realtype* qv, realtype* qc){

  return dlc::Cp_dry + dlc::Cp_v*(*qv) + dlc::C_l*(*qc);

}


/* Calculate dimensionless Fkl and Fdl 
* heat and vapour diffusion factors in 
equation for radial growth of droplets 
according to eq.7.27 lohmann, luond and 
mahrt intro 2 clouds textbook 
*/
int diffusion_factors(realtype* fkl, realtype* fdl, 
            realtype temp, realtype p, realtype psat)
{
  realtype Thermk, Diffuse_v;
  realtype Temp = temp*dlc::TEMP0;
  realtype P = p*dlc::P0;
  realtype Psat = psat*dlc::P0;
  
  Thermk = 7.11756e-5*pow(Temp, 2.0) + Temp*4.38127686e-3;            // [eq.7.24]
  Diffuse_v = (4.012182971e-5 / P * pow(Temp,1.94))/DC::RGAS_V;

  *fkl = (DC::LATENT_RGAS_V/Temp-1)*DC::LATENT_V/(Thermk*dlc::F0); 
  *fdl = Temp/(Diffuse_v*Psat)/dlc::F0;
    
  return 0;
}
/* ----------------------------------------------- */




/* -------------- ODE functions -------------- */


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
  realtype dtemp, rho_d, cp_m;
  
  rho_d = dlc::Mr_ratio/(dlc::Mr_ratio+qv)*p/temp;      // density of dry parcel (p_dry/temp)    
  
  cp_m = cp_moist(&qv, &qc);

  dtemp = Ith(ydot,2) = dlc::Rgas_dry/(rho_d*cp_m) * Ith(ydot,1);

  return(dtemp);
}






/*
 * dr/dt differential equation describing radius
  evolution over time for each superdrop due to diffusion
  and condensation of water vapour.
 */

static int condensation_droplet_growth(realtype delt, realtype* p, 
    realtype* temp, realtype* qv, realtype*  qc,
    realtype* dtemp, realtype* dqv, realtype* dqc,
    Superdrop* ptr, int nsupers)
{
  realtype dr, pv, psat, s_ratio, r, eps, a, b, fkl, fdl;
  realtype tot_drhov;

  psat = saturation_pressure(*temp);
  pv = (*p)*(*qv)/(dlc::Mr_ratio+(*qv));
  s_ratio = pv/psat;              // parcel supersaturation ratio 
 
  /* radial growth/shrink each droplet 
  [eq.7.27 lohmann intro 2 clouds] */
  fkl=0.0; fdl=0.0; tot_drhov=0.0;
  SDloop(i, nsupers)
  { 
    r = (ptr+i) -> r;
    eps = (ptr+i) -> eps;
    a = (ptr+i) -> akohler_factor(*temp);
    b = (ptr+i) -> bkohler_factor();
    diffusion_factors(&fkl, &fdl, *p, *temp, psat);
    dr = (s_ratio-1 -a/r + b/pow(r, 3.0)) / (dlc::Rho_l * (fkl+fdl) * r);
    //dr = (s_ratio-1) / (dlc::Rho_l * (fkl+fdl) * r);
    
    /*  if droplets are dry, do not shrink further */
    if (r<= (ptr+i) -> dry_r() && dr<=0.0){
      dr = 0.0;
    }

    (ptr+i) -> r += dr*delt;

    realtype dm = dm_dt_const*pow(r,2.0)*eps*dr*delt;     // dm/dt * delta t
    tot_drhov += dm;                               // drho_condensed_vapour/dt * delta t
  }
  
  // *dqc += tot_drhov/dlc::Rho_dry; 
  *dqc += 0;
  *dqv += -(*dqc); 
  *dtemp +=  (dlc::Latent_v/cp_moist(qv, qc))*(*dqc);

  return(0);
}






/*
 * Simple function f(t,y) called by ODE solver to 
 solve differential equations over time.
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  
  UserData data = (UserData) user_data;

  // bool doCond;
  // int nsupers;
  // Superdrop* ptr;
  // doCond = data -> doCond;
  // nsupers = data -> nsupers;
  // ptr = data -> ptr;
  
  if(doExpand){
    dp_dt(t, ydot, data -> w);     
    dtemp_dt_adia(ydot, Ith(y,1), Ith(y,2), Ith(y,3), Ith(y,4));    
  }

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
