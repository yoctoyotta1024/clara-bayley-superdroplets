// Author: Clara Bayley
// File: mainSDM.cpp
/* This file runs the entire superdrop model (SDM)
coupled with a CVODE ode solver for the kinetics
(p, temp, qv and qc) over time*/


//#include <iostream>

#include "claras_SDinit.hpp"
#include "claras_SDconstants.hpp"

#include "common2allsuperdrops.hpp"
#include "superdrop.hpp"
#include "readwritefuncs4SDM.hpp"
#include "collisions4SDM.hpp"
#include "condensationgrowth4SDM.hpp"


//using namespace std;
namespace dlc = dimless_constants;




int main()
{
  ofstream yfile, sdfile;                                     // files to write superdroplet and ODE output data to

  Superdrop superdrops_arr[init::NSUPERS];
  const int nsupers = init::NSUPERS;                           // allows flexible no. of superdroplets to be included in model (nsupers <= NSUPERS)
  double deltemp, delqv, delqc;                                // changes to kinematics due to superdroplets
  double p, temp, qv, qc;                                      // copies of kinematic variables from ode sovler passed to SDM 

  const int nout = init::nout;                                 // number of output times
  const double t0 = init::TSPAN[0]/dlc::TIME0;                       // initial time (dimensionless)          
  const double tstep = init::TSPAN[1]/nout/dlc::TIME0;               // output time step (dimensionless)     

  const double cond_tstep = init::COND_TSTEP/dlc::TIME0;                   
  const double coll_tstep = init::COLL_TSTEP/dlc::TIME0;                    
  const double min_tstep = min(coll_tstep, cond_tstep);         // smallest timestep is one to increment time by
  const double ItersPerTstep = tstep/min_tstep;                 // no. of increments = ceil(TSTEP/min_tstep) >= 1
  double tout, delta_t, dt_cond, dt_coll;
  
  /* Get nhalf, scale_p and pvec (index list) given constant no. of nsupers */
  const int nhalf = floor(nsupers/2.0);
  vector<int> pvec(nsupers);
  for(int i=0; i<nsupers; i++)
  {
    pvec[i] = i;
  }
  
  /* Initialise Superdroplets using initdrops_csv file */
  for(int i=0; i<nsupers; i++)                  //for loop over all superdroplets
  {
    superdrops_arr[i] = Superdrop(init::iRho_l, 
                init::iRho_sol, init::iMr_sol, init::iIONIC); 
  }
  initialise_Superdrop_instances(init::initdrops_csv,
                                  superdrops_arr, nsupers);


  /* Initial qv and qc from relative humidity */
  const double p_init = init::P_INIT/dlc::P0;                 // initial (dimensionless) kinetics
  const double temp_init = init::TEMP_INIT/dlc::TEMP0;
  const double pv_i = saturation_pressure(temp_init)*init::relh_init/100;
  const double qv_init = pv2qv(pv_i, p_init);                 //initial qv for solver
  const double qc_init = init::qc_init;

  double t = 0;    //!!!!!!!!!!@@@@@@@@@@@###########     
  double yvec[4] = {p_init, temp_init, qv_init, qc_init};                             



  /* write header to .csv files and open in preparation for writing data */ 
  write_outputheader(init::solution_csv);
  write_superdrop_outputheader(init::solutionSD_csv); 
  yfile.open(init::solution_csv, ios::app);
  sdfile.open (init::solutionSD_csv, ios::app);
  write_output(yfile, t, yvec);
  write_superdrop_output(sdfile, superdrops_arr, nsupers);


  /* RUN SUPERDROPLET MODEL COUPLED TO CVODE ODE SOLVER FOR KINEMATICS
    (collecting data nout no. of times within init::TSPAN) */
  tout = t0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver
  dt_cond = 0;
  dt_coll = 0;
  
  for (int j=0; j<nout; j++)
  {
    /* variables and their changes due to SDM at timestep */
    // cout << " -- t of SDs: " << t << endl;
    print_output(t, yvec);
    
    /* RUN SUPERDROPLET MODEL */
    /* kinematic variables coupled to superdroplets */
    variables_b4tstep(yvec, p, temp, qv, qc, deltemp, delqv, delqc);

    /* a) divide each tstep into ItersPerTstep no. of 
      min_tsteps and model SD collisions and/or condensational growth */
    delta_t = 0;
    for(int k=0; k<ceil(ItersPerTstep); k++)                 // increment time for SDs simulation
    {
      delta_t += min_tstep;
      dt_cond += min_tstep;                                  // change in time since last condensation event
      dt_coll += min_tstep;                                  // change in time since last collision event

      /* (b)) Superdroplet Condensation-Diffusion Growth */
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

    /* RUN CVODE ODE SOLVER */
    t = tout ;     //!!!!!!!!!!@@@@@@@@@@@###########     
    

    /* Continute to next timestep */
    if(init::doCond)
    {
      yvec[1] += deltemp;
      yvec[2] += delqv;
      yvec[3] += delqc; 
      // /* Reinitialize the solver after discontinuous change in y */
      // retval = CVodeReInit(cvode_mem, tout, y);
      // if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);
    }

    //if (retval == CV_SUCCESS)
    //{
    tout += delta_t;
    //}



    /* Output solution and error after every large timestep */
    write_superdrop_output(sdfile, superdrops_arr, nsupers);
    write_output(yfile, t, yvec); 
    
  }


  sdfile.close();




  return 0;
}




