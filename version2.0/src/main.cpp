// Author: Clara Bayley
// File: mainSDM.cpp
/* This file runs the entire superdrop model (SDM)
coupled with a CVODE ode solver for the kinetics
(p, temp, qv and qc) over time*/


#include <iostream>
#include <math.h>

#include "../claras_SDinit.hpp"
#include "../claras_SDconstants.hpp"

#include "./superdrop_solver/common2allsuperdrops.hpp"
#include "./superdrop_solver/superdrop.hpp"
#include "./superdrop_solver/readwritefuncs4SDM.hpp"
#include "./superdrop_solver/collisions4SDM.hpp"
#include "./superdrop_solver/condensationgrowth4SDM.hpp"

#include "./kinetics_solver/cvode_odesolver.hpp"
#include "./kinetics_solver/differential_functions.hpp"


using namespace std;
namespace dlc = dimless_constants;




int main(){
  
  /* variables for time-stepping model */
  ofstream yfile, sdfile;                                       // files to write superdroplet and ODE output data to
  int cvode_iterfail;                                           // flag to break integration if cvode fails
  const int nout = init::nout;                                  // number of output times
  const double t0 = init::TSPAN[0]/dlc::TIME0;                  // initial time (dimensionless)          
  const double tstep = init::TSPAN[1]/nout/dlc::TIME0;          // output time step (dimensionless)     
 
  const double cond_tstep = init::COND_TSTEP/dlc::TIME0;                   
  const double coll_tstep = init::COLL_TSTEP/dlc::TIME0;                    
  const double min_tstep = min(coll_tstep, cond_tstep);         // smallest timestep is one to increment time by
  const double ItersPerTstep = tstep/min_tstep;                 // no. of increments = ceil(TSTEP/min_tstep) >= 1
  double tout, delta_t, dt_cond, dt_coll;

  /* variables required for superdroplet model */
  Superdrop superdrops_arr[init::NSUPERS];
  const int nsupers = init::NSUPERS;                           // allows flexible no. of superdroplets to be included in model (nsupers <= NSUPERS)
  double deltemp, delqv, delqc;                                // changes to kinematics due to superdroplets
  double p, temp, qv, qc;                                      // copies of kinematic variables from ode sovler passed to SDM 

  /* initial conditions for ODEs */ 
  const double w = init::iW/dlc::W0;                          // dimensionless w velocity for f(t,y,ydot,w,...)
  const double p_init = init::P_INIT/dlc::P0;                 // initial (dimensionless) kinetics
  const double temp_init = init::TEMP_INIT/dlc::TEMP0;
  const double pv_i = saturation_pressure(temp_init)*init::relh_init/100;
  const double qv_init = pv2qv(pv_i, p_init);                 //initial qv for solver
  const double qc_init = init::qc_init;



  /* Initialise Superdroplets using initdrops_csv file */
  for(int i=0; i<nsupers; i++)                  //for loop over all superdroplets
  {
    superdrops_arr[i] = Superdrop(init::iRho_l, 
                init::iRho_sol, init::iMr_sol, init::iIONIC); 
  }
  initialise_Superdrop_instances(init::initdrops_csv,
                                  superdrops_arr, nsupers);

  /* Get nhalf, scale_p and pvec (index list) given constant no. of nsupers */
  const int nhalf = floor(nsupers/2.0);
  vector<int> pvec(nsupers);
  for(int i=0; i<nsupers; i++)
  {
    pvec[i] = i;
  }
  
  /* setup CVODE ODE solver */
  CvodeOdeSolver cvode;
  cvode.init_userdata(w, init::doThermo);
  const double y_init[4] = {p_init, temp_init, qv_init, qc_init};
  cvode.setup_ODE_solver(init::rtol, init::atols, y_init, t0);
  cvode.print_init_ODEdata(nout, t0, ceil(ItersPerTstep)*min_tstep, tstep);
    
  /* assign variables from CVODE ODE solver to SDM kinematic variables */
  double &t = cvode.t;
  cvode.get_variables_b4tstep(p, temp, qv, qc, deltemp, delqv, delqc);

  /* write header to .csv files and open in preparation for writing data */ 
  write_outputheader(init::solution_csv);
  write_superdrop_outputheader(init::solutionSD_csv); 
  yfile.open(init::solution_csv, ios::app);
  sdfile.open (init::solutionSD_csv, ios::app);
  write_output(yfile, t, p, temp, qv, qc);
  write_superdrop_output(sdfile, superdrops_arr, nsupers);


  /* run superdroplet model (SDM) coupled to 
    CVODE ODE SOLVER for kinematics (collecting data 
    nout no. of times within init::TSPAN) */
  tout = t0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver
  dt_cond = 0;
  dt_coll = 0;
  for (int j=0; j<nout; j++)
  {

    /* SDM kinematic variables at timestep */
    // cout << " -- t of SDM: " << t << endl;
    cvode.get_variables_b4tstep(p, temp, qv, qc, deltemp, delqv, delqc);
    print_output(t, p, temp, qv, qc);
    
    /* SECTION 1: run SUPERDROPLET model */
    /* (a) divide each tstep into ItersPerTstep no. of 
      min_tsteps and model SD collisions and/or condensational growth */
    delta_t = 0;
    for(int k=0; k<ceil(ItersPerTstep); k++)                 // increment time for SDs simulation
    {
      delta_t += min_tstep;
      dt_cond += min_tstep;                                  // change in time since last condensation event
      dt_coll += min_tstep;                                  // change in time since last collision event

      /* (b) Superdroplet Condensation-Diffusion Growth */
      if(init::doCond){
        if(dt_cond >= cond_tstep)
        {
          //cout << "cond @ " << t+delta_t << endl;
          condensation_onto_superdroplets(cond_tstep, p, temp, qv, qc,
              deltemp, delqv, delqc, superdrops_arr, nsupers);
          dt_cond=0;
        }
      }

      /* (c) Superdroplet Collisions */
      if (init::doColl)
      {
        if(dt_coll >= coll_tstep)
        {
          //cout << "coll @ " << t+delta_t << endl;
          collide_droplets(nsupers, nhalf, pvec, superdrops_arr);
          dt_coll=0;
        }
      }
        
    }


    /* SECTION 2: run CVODE ODE solver */
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


    /* SECTION 3: write data and proceed to next time step */
    if(cvode_iterfail){ break; }
    tout += delta_t;

    /* Output solution and error after every large timestep */
    write_superdrop_output(sdfile, superdrops_arr, nsupers);
    write_output(yfile, t, p, temp, qv, qc);


  }




  /* end CVODE ODE solver and close data files */
  sdfile.close();
  yfile.close();
  cvode.destroy_cvode();



  return 0;
}