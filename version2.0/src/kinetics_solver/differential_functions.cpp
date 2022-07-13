// Author: Clara Bayley
// File: differential_functions.cpp
/* ODE functions which are solved by
 CVODE ode solver to model evolution of the
 kinetics (p, temp, qv and qc) over time */

#include "differential_functions.hpp"

int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */
{

  UserData data = (UserData)user_data;

  if (data->doThermo)
  {
    double pdot;
    pdot = dp_dt(t, data->w);
    NV_Ith_S(ydot, 0) = pdot;
    NV_Ith_S(ydot, 1) = dtemp_dt_adia(pdot, y);
  }
  else
  {
    NV_Ith_S(ydot, 0) = 0;
    NV_Ith_S(ydot, 1) = 0;
  }
 
  NV_Ith_S(ydot, 2) = 0;
  NV_Ith_S(ydot, 3) = 0;

  return 0;
}

static double dp_dt(const double &t, const double &w)
/* dp/dt differential equation (dimensionless)
  describing pressure evolution over time.
  note: true dP/dt = dp/dt * P0/TIME0 */
{
  double pdot, profile;

  static const double dp_dt_const = dlc::W0 * dlc::TIME0 * DC::G / (DC::RGAS_DRY * dlc::TEMP0); 
  static const double zg = 0.0 / (dlc::W0 * dlc::TIME0);                    // dimensionless z value at ground level
  static const double lpsrate = 0.0062 / dlc::TEMP0 * dlc::W0 * dlc::TIME0; // dimensionless moist adiabatic lapse rate
  static const double tempg = 273.15 / dlc::TEMP0;                          // dimensionless temperature at zg
  static const double pg = 100000.0 / dlc::P0;                              // dimensionless pressure at zg
  static const double gamma = DC::G / (DC::RGAS_DRY * 0.0062) - 1.0;        // constant in dry adiabatic expansion

  profile = 1.0 - lpsrate / tempg * (w * t - zg); // characteristic function for pressure profile as
  profile = pow(profile, gamma);                  //      a funciton of time (ie. height via z=w*t)

  pdot = -dp_dt_const * pg / tempg * profile;

  return pdot;
}

static double cp_moist(const double &qv, const double &qc)
/* effective specific heat capacity of moist parcel
  of air (dry + water vapour + liquid water) */
{
  return dlc::Cp_dry + dlc::Cp_v * qv + dlc::C_l * qc;
}

static double dtemp_dt_adia(const double &pdot, const N_Vector &y)
/* dtemp/dt differential equation describing
  temperature evolution solely due to pressure
  changes in parcel for adiabatic process (no heat loss).
  Parcel has water vapour mass mixing ratio (m_v/m_dry) = qv and
  liquid water mass mixing ratio (m_c/m_dry) = qc.
  note: True dTemp/dt = dtemp * TEMP0/TIME0  */
{
  double tempdot, rho_d, cp_m;
  double p, temp, qv, qc;

  p = NV_Ith_S(y, 0);
  temp = NV_Ith_S(y, 1);
  qv = NV_Ith_S(y, 2);
  qc = NV_Ith_S(y, 3);

  rho_d = dlc::Mr_ratio / (dlc::Mr_ratio + qv) * p / temp; // density of dry parcel (p_dry/temp)

  cp_m = cp_moist(qv, qc); // moist specific heat capacity

  tempdot = dlc::Rgas_dry / (rho_d * cp_m) * pdot;

  return tempdot;
}
