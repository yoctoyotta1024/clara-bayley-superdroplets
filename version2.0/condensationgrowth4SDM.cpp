// Author: Clara Bayley
// File: condensationgrowth4SDM.cpp
/* Functionality for modelling condensation- 
   diffusional growth of superdroplets */

#include "condensationgrowth4SDM.hpp"



double pv2qv(const double pv, const double p)
	/* Calculate mass mixing ratio
	qv = m_v/m_dry = rho_v/rho_dry
	given vapour pressure pv = p_v/p_tot. */
{  
  return dlc::Mr_ratio * pv/(p-pv);
}



double cp_moist(const double qv, const double qc)
	/* (dimensionless) specific heat capacity of 
		moist parcel of air */
{
  return dlc::Cp_dry + dlc::Cp_v*(qv) + dlc::C_l*(qc);
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
  T = temp*dlc::TEMP0;                               // real T [K]
  
  if(T <= 0)
  {
      std::cout << "psat ERROR: T must be larger than 0K. T = " <<std::endl;
  }

  lnpsat = (54.842763                               // ln(psat) [Pa]
        - 6763.22 / T
        - 4.21 * log(T)
        + 0.000367 * T
        + tanh(0.0415 * (T - 218.8))
        * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T));

  return exp(lnpsat)/dlc::P0;                      // dimensionless psat
}






void condensation_onto_superdroplets(const double delt, double &p, 
    double &temp, double &qv, double qc,
    double &deltemp, double &delqv, double &delqc,
    Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers)
	/* Change to temp, qv and qc due to change in radius of
	all superdroplets via diffusion and condensation 
	of water vapour during timestep delt. */
{

  double deltemp_k, delqv_k, delqc_k;
  double epsdelr, delm, psat, s_ratio;
	double tot_delrhov = 0;
	static const double eqnc = 4*M_PI*dlc::Rho_l*pow(dlc::R0, 3.0);
	static const double dmdt_const = eqnc/init::DROPVOL;                     //constant value in dm_dt diffusion growth equation

  psat = saturation_pressure(temp);
  s_ratio = (p*qv)/(dlc::Mr_ratio+(qv) * psat);                            // parcel supersaturation ratio s = pv/psat 

  /* superdroplet radii changes for timestep delt */
  for(int i=0; i<nsupers; i++)
  { 
		epsdelr = superdrops_arr[i].condensation_growth(p, temp,
					 														psat, s_ratio, delt);
		delm = dmdt_const*pow(superdrops_arr[i].r,2.0)*epsdelr;                // dm/dt * delta t
    tot_delrhov += delm;                                                   // drho_condensed_vapour/dt * delta t
  }
  delqc_k = tot_delrhov/dlc::Rho_dry;                                      // change to temp, qv and qc as a result of 
  delqv_k =  -(delqc_k);                                                     //       condensation at kth small timestep
  deltemp_k =  (dlc::Latent_v/cp_moist(qv, qc))*(delqc_k);

  /* temp, qv and qc change for timestep delt */
  temp += deltemp_k;
  qv += delqv_k;
  qc += delqc_k;

  /* cumulative changes from k=0 to k=k time step 
  (calculated for coupling to kinematics ODE solver) */
  delqc += deltemp_k; 
  delqv += delqv_k;
  deltemp += delqc_k;

}




void variables_b4tstep(const double y[4], double &p, double &temp, double &qv, 
      double &qc, double &deltemp, double &delqv,
      double &delqc)
	/* copies values of kinematic variables from CVODE
  ODE sovler (p, temp, qv, qc) to SDM before timesteps of
  SDM begin. Also resets deltas=0 for start of timestep */
{

  // double p = NV_Ith_S(y,0);                                             
  // double temp = NV_Ith_S(y,1); 
  // double qv = NV_Ith_S(y,2);
  // double qc = NV_Ith_S(y,3);  

  p = y[0];                                       
  temp = y[1];     
  qv = y[2];     
  qc = y[3];     
                                            
  deltemp = 0;
  delqv = 0;
  delqc = 0;
}