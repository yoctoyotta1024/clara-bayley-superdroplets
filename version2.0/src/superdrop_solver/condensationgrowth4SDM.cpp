// Author: Clara Bayley
// File: condensationgrowth4SDM.cpp
/* Functionality for modelling condensation-
   diffusional growth of superdroplets. Equations
   referenced as (eqn [X.YY]) are from "An Introduction
   To Clouds From The Microscale to Climate" by
   Lohmann, Luond and Mahrt, 1st edition. */

#include "condensationgrowth4SDM.hpp"

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

void condensation_onto_superdroplets(const double delt, double &p,
              double &temp, double &qv, double qc, 
              double &deltemp, double &delqv, double &delqc, 
              Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers)
/* Change to superdroplet radii and temp, qv and 
  qc due to sum of radii changes via diffusion and 
  condensation of water vapour during timestep delt. 
  Using equations from "An Introduction To 
  Clouds...." (see note at top of file) */
{

  double deltemp_k, delqv_k, delqc_k;
  double epsdelr, delm, psat, s_ratio;
  double tot_delrhov = 0.0;
  static const double eqnc = 4.0 * M_PI * dlc::Rho_l * pow(dlc::R0, 3.0);
  static const double dmdt_const = eqnc / init::DROPVOL;

  psat = saturation_pressure(temp);
  s_ratio = (p * qv) / (dlc::Mr_ratio + (qv)*psat); // supersaturation ratio s = pv/psat

  /* superdroplet radii changes for timestep delt */
  for (int i = 0; i < nsupers; i++)
  {
    epsdelr = superdrops_arr[i].condensation_growth(p, temp,
                                                    psat, s_ratio, delt);
    delm = dmdt_const * pow(superdrops_arr[i].r, 2.0) * epsdelr; // dm/dt * delta t (eqn [7.22])
    tot_delrhov += delm;                                         // drho_condensed_vapour/dt * delta t

  }
  
  delqc_k = tot_delrhov / dlc::Rho_dry; // change to temp, qv and qc as a result of
  delqv_k = -(delqc_k);                 //    condensation at kth (small) timestep
  deltemp_k = (dlc::Latent_v / cp_moist(qv, qc)) * (delqc_k);

  /* update to temp, qv and qc given timestep delt */
  temp += deltemp_k;
  qv += delqv_k;
  qc += delqc_k;

  /* cumulative changes from k=0 to k=k time step
  (calculated for coupling kinematics to CVODE ODE solver) */
  delqc += deltemp_k;
  delqv += delqv_k;
  deltemp += delqc_k;

}
