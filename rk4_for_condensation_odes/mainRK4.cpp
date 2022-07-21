// Author: Clara Bayley
// File: mainRK4.cpp
/* This file solves condensation-diffusion 
growth equation for superdroplets using 
runge kuuta 4th order method */

#include <iostream>
#include <math.h>

/* coupled model parameter & input files */
#include "../version2/claras_SDinit.hpp"
#include "../version2/claras_SDconstants.hpp"

#include "runge_kutta_condensation.hpp"

/* Superdroplet Model (SDM) files */
#include "../version2/src/superdrop_solver/common2allsuperdrops.hpp"
#include "../version2/src/superdrop_solver/superdrop.hpp"
#include "../version2/src/superdrop_solver/readwritefuncs4SDM.hpp"
//#include "../version2/src/superdrop_solver/collisions4SDM.hpp"
//#include "../version2/src/superdrop_solver/condensationgrowth4SDM.hpp"

// using namespace std;
// namespace dlc = dimless_constants;



int main()
{

  /* ---------- same as main() in ../version2/src/ -------------- */
  /* files to write output data to */
  const std::string initdrops_csv2 = "../../version2/dimlessSDinit.csv"; //.csv filename for initialisation of superdrops

  ofstream yfile, sdfile; // files to write superdroplet and ODE output data to
  string setupfiles[] = {"../../version2/claras_SDinit.hpp",
                         "../../version2/claras_SDconstants.hpp"}; // files to copy to setup_txt file

  /* variables for time-stepping model */
  const int nout = init::nout;                             // number of output times
  const double t0 = init::TSPAN[0] / dlc::TIME0;           // initial time (dimensionless)
  const double tstep = init::TSPAN[1] / nout / dlc::TIME0; // output time step (dimensionless)

  const double cond_tstep = init::COND_TSTEP / dlc::TIME0;
  const double coll_tstep = init::COLL_TSTEP / dlc::TIME0;
  const double min_tstep = min(coll_tstep, cond_tstep); // smallest timestep is one to increment time by
  const double ItersPerTstep = tstep / min_tstep;       // no. of increments = ceil(TSTEP/min_tstep) >= 1
  double tout, delta_t, dt_cond, dt_coll;

  /* variables required for superdroplet model */
  Superdrop superdrops_arr[init::NSUPERS];
  const int nsupers = init::NSUPERS; // allows flexible no. of superdroplets to be included in model (nsupers <= NSUPERS)
  double deltemp, delqv, delqc;      // changes to kinematics due to superdroplets
  double p, temp, qv, qc;            // copies of kinematic variables from ode sovler passed to SDM

  /* initial conditions for ODEs */
  const double p_init = init::P_INIT / dlc::P0; // initial (dimensionless) kinetics
  const double temp_init = init::TEMP_INIT / dlc::TEMP0;
  const double pv_i = saturation_pressure(temp_init) * init::relh_init / 100;
  const double qv_init = pv2qv(pv_i, p_init);    // initial qv for solver
  const double qc_init = init::qc_init;

  /* Initialise Superdroplets using initdrops_csv file */
  for (int i = 0; i < nsupers; i++)
  {
    superdrops_arr[i] = Superdrop(init::iRho_l,
                                  init::iRho_sol, init::iMr_sol, init::iIONIC);
  }
  initialise_Superdrop_instances(initdrops_csv2,
                                 superdrops_arr, nsupers);

  // /* Get nhalf and pvec (list of indicies) given constant no. of nsupers */
  // const int nhalf = floor(nsupers / 2.0);
  // vector<int> pvec(nsupers);
  // for (int i = 0; i < nsupers; i++)
  // {
  //   pvec[i] = i;
  // }

  /* ------------------------------------------------------------ */

  double t = 0;
  p = p_init;
  temp = temp_init; 
  qv = qv_init;
  qc = qc_init;
  deltemp = 0;
  delqv = 0;
  delqc = 0;

  /* write setup to .txt file */
  WriteSetup2Txt(setupfiles, init::setup_txt);

  /* write header and initial data to .csv files */
  write_outputheader(init::solution_csv);
  write_superdrop_outputheader(init::solutionSD_csv);
  yfile.open(init::solution_csv, ios::app);
  sdfile.open(init::solutionSD_csv, ios::app);
  write_output(yfile, t, p, temp, qv, qc);
  write_superdrop_output(sdfile, superdrops_arr, nsupers);


  /* *** run superdroplet model (SDM) coupled to
    CVODE ODE SOLVER for kinematics. Writing output
    data nout no. of times within init::TSPAN) *** */
  tout = t0 + ceil(ItersPerTstep) * min_tstep; // first output time of ODE solver
  dt_cond = 0;
  dt_coll = 0;
  for (int j = 0; j < nout; j++)
  {
    print_output(t, p, temp, qv, qc); // print variables at timestep t

    /* SECTION 1: run SUPERDROPLET model */
    /* (a) divide each tstep into ItersPerTstep no. of min_tsteps and
      model SD collisions and/or condensational growth for each min_tstep */
    delta_t = 0;
    for (int k = 0; k < ceil(ItersPerTstep); k++) // increment time of SDM
    {
      delta_t += min_tstep;
      dt_cond += min_tstep; // change in time since last condensation event
      dt_coll += min_tstep; // change in time since last collision event

      /* (b) model superdroplet condensation-diffusion growth */
      if (init::doCond)
      {
        if (dt_cond >= cond_tstep)
        {
          condensation_onto_superdroplets_RK4(cond_tstep, p, temp, qv, qc,
                            deltemp, delqv, delqc, superdrops_arr, nsupers);
          dt_cond = 0; //reset time since last condensation event
        }
      }

    }
  write_superdrop_output(sdfile, superdrops_arr, nsupers);
  write_output(yfile, t, p, temp, qv, qc);
  t = tout;
  tout += delta_t;

  }





  return 0;
}