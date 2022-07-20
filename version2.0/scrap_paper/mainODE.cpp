// Author: Clara Bayley
// File: mainSDM.cpp
/* This file runs the the CVODE ode solver 
to solve the ODEs for the kinetics
(p, temp, qv and qc) over time*/



#include <iostream>
//#include <stdio.h>
#include <math.h>

#include "../claras_SDinit.hpp"
#include "../claras_SDconstants.hpp"
#include "cvode_odesolver.hpp"


using namespace std;


/* ------------------------------- Helper Functions ------------------------------- */

double vapourpressure_2_massmixratio(const double pv, const double p)
	/* Calculate mass mixing ratio
	qv = m_v/m_dry = rho_v/rho_dry
	given vapour pressure pv = p_v/p_tot. */
{  
  return dlc::Mr_ratio * pv/(p-pv);
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







/* ------------------------------- Main Program ------------------------------- */
int main(){

  int cvode_iterfail;                                        // flag to break integration if cvode fails
 
  const int nout = init::nout;                                 // number of output times
  const double t0 = init::TSPAN[0]/dlc::TIME0;                       // initial time (dimensionless)          
  const double tstep = init::TSPAN[1]/nout/dlc::TIME0;               // output time step (dimensionless)     

  const double cond_tstep = init::COND_TSTEP/dlc::TIME0;                   
  const double coll_tstep = init::COLL_TSTEP/dlc::TIME0;                    
  const double min_tstep = min(coll_tstep, cond_tstep);         // smallest timestep is one to increment time by
  const double ItersPerTstep = tstep/min_tstep;                 // no. of increments = ceil(TSTEP/min_tstep) >= 1
  double tout, delta_t, dt_cond, dt_coll;
  
  /* initial conditions & setup of model */  
  const double w = init::iW/dlc::W0;                          // dimensionless w velocity for f(t,y,ydot,w,...)
  const double p_init = init::P_INIT/dlc::P0;                 // initial (dimensionless) kinetics
  const double temp_init = init::TEMP_INIT/dlc::TEMP0;
  const double pv_i = saturation_pressure(temp_init)*init::relh_init/100;
  const double qv_init = vapourpressure_2_massmixratio(pv_i, p_init);                 // Initial qv from relative humidity
  const double qc_init = init::qc_init;


  /* SETUP CVODE ODE SOLVER */
  CvodeOdeSolver cvode;

  cvode.init_userdata(w, init::doThermo);
  const double y_init[4] = {p_init, temp_init, qv_init, qc_init};
  cvode.setup_ODE_solver(init::rtol, init::atols, y_init, t0);
  
  
  /* START CVODE ODE SOLVER */
  //tout = T0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver
  cvode.print_init_ODEdata(nout, t0, ceil(ItersPerTstep)*min_tstep, tstep);


  /* RUN CVODE ODE SOLVER FOR KINEMATICS
    (collecting data nout no. of times within init::TSPAN) */
  tout = t0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver
  dt_cond = 0;
  dt_coll = 0;
  for (int j=0; j<nout; j++)
  {

    /* a) divide each tstep into ItersPerTstep no. of 
      min_tsteps and model SD collisions and/or condensational growth */
    delta_t = 0;
    for(int k=0; k<ceil(ItersPerTstep); k++)                 // increment time for SDs simulation
    {
      delta_t += min_tstep;        
    }

    cvode_iterfail = cvode.advance_solution(tout);
    if(cvode_iterfail){ break; }
    
    /* Continute to next timestep */ 
    // if(init::doCond)
    // {
    //   yvec[1] += deltemp;
    //   yvec[2] += delqv;
    //   yvec[3] += delqc; 

    //   cvode_iterfail = cvode.reinitialise(tout, y);
    // }

    if(cvode_iterfail){ break; }
    tout += delta_t;

  }

  /* END CVODE ODE SOLVER */
  cvode.destroy_cvode();




 return(0);
}

