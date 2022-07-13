// Author: Clara Bayley
// File: runge_kutta_condensation.hpp
/* methods for solving condensation-diffusion 
growth equation for superdroplets using 
runge kuuta 4th order method */

#ifndef RUNGE_KUTTA_CONDENSATION_HPP
#define RUNGE_KUTTA_CONDENSATION_HPP


#include <iostream>
#include <math.h>

#include "../version2.0/claras_SDinit.hpp"
#include "../version2.0/claras_SDconstants.hpp"
#include "../version2.0/src/superdrop_solver/common2allsuperdrops.hpp"
#include "../version2.0/src/superdrop_solver/superdrop.hpp"

// using namespace std;
namespace dlc = dimless_constants;

double pv2qv(const double pv, const double p)
/* Calculate mass mixing ratio
  qv = m_v/m_dry = rho_v/rho_dry
  given vapour pressure pv = p_v/p_tot. */
{
  return dlc::Mr_ratio * pv / (p - pv);
}

double cp_moist(const double qv, const double qc)
/* (dimensionless) specific heat 
  capacity of moist parcel of air */
{
  return dlc::Cp_dry + dlc::Cp_v * (qv) + dlc::C_l * (qc);
}

double saturation_pressure(const double temp)
/* Calculate the equilibrium vapor pressure
  of water over liquid water ie. the
  saturation pressure (psat). Equation taken from
  python module typhon.physics.thermodynamics.e_eq_water_mk
  with conversion to real temp /K = T*Temp0 and from
  real psat to dimensionless psat = psat/P0. */
{
  double T, lnpsat;
  T = temp * dlc::TEMP0; // real T [K]

  if (T <= 0.0)
  {
    std::cout << "psat ERROR: T must be larger than 0K. T = " << std::endl;
  }

  lnpsat = (54.842763 // ln(psat) [Pa]
      - 6763.22 / T - 4.21 * log(T) + 0.000367 * T 
      + tanh(0.0415 * (T - 218.8)) 
      * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T));

  return exp(lnpsat) / dlc::P0; // dimensionless psat
}


void condensation_onto_superdroplets_RK4(const double delt, double &p,
              double &temp, double &qv, double qc, 
              double &deltemp, double &delqv, double &delqc, 
              Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers)
/* Change to superdroplet radii and temp, qv and 
  qc due to sum of radii changes via diffusion and 
  condensation of water vapour during timestep delt. 
  Using equations from "An Introduction To 
  Clouds...." (see note at top of file) */
{

  const double psat = saturation_pressure(temp);
  const double s_ratio = (p * qv) / (dlc::Mr_ratio + (qv)*psat); // supersaturation ratio s = pv/psat
  double epsdelr;

  /* superdroplet radii changes for timestep delt */
  for (int i = 0; i < nsupers; i++)
  {
    epsdelr = superdrops_arr[i].runge_kutta_condensation_growth(p, temp, psat, s_ratio, delt); 
  }
  
}




#endif //RUNGE_KUTTA_CONDENSATION_HPP