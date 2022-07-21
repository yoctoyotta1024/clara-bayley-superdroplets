// Author: Clara Bayley
// File: mainSDM.cpp
/* This file runs the entire superdrop model (SDM)
coupled with a CVODE ode solver for the kinetics
(p, temp, qv and qc) over time */

#include <iostream>
#include <math.h>

/* coupled model parameter & input files */
#include "../claras_SDinit.hpp"
#include "../claras_SDconstants.hpp"

/* Superdroplet Model (SDM) files */
#include "./superdrop_solver/common2allsuperdrops.hpp"
#include "./superdrop_solver/superdrop.hpp"
#include "./superdrop_solver/readwritefuncs4SDM.hpp"
#include "./superdrop_solver/collisions4SDM.hpp"
#include "./superdrop_solver/condensationgrowth4SDM.hpp"

/* CVODE ODE kinetics solver files */
#include "./kinetics_solver/cvode_odesolver.hpp"
#include "./kinetics_solver/differential_functions.hpp"

using namespace std;
namespace dlc = dimless_constants;

int main()
{

  /* files to write output data to */
  ofstream kinetics_datafile, superdrops_datafile; // files to write superdroplet and ODE output data to
  string setupfiles2copy[] = {"../claras_SDinit.hpp",
                         "../claras_SDconstants.hpp"}; // files with contents to copy to setup_txt output file

  /* variables for time-stepping model */
  int is_cvodefail;                                        // flag to break integration if CVODE ODE solver fails
  const int nout = init::nout;                             // number of output times
  const double t_init = init::TSPAN[0] / dlc::TIME0;       // initial time (dimensionless)
  const double tstep = init::TSPAN[1] / nout / dlc::TIME0; // output time step (dimensionless)

  const double cond_tstep = init::COND_TSTEP / dlc::TIME0;
  const double coll_tstep = init::COLL_TSTEP / dlc::TIME0;
  const double min_tstep = min(coll_tstep, cond_tstep); // smallest timestep is one to increment time by
  const double iters_per_tstep = tstep / min_tstep;     // no. of increments = ceil(TSTEP/min_tstep) >= 1
  double tout;                                          // time at which data is written to output datafiles
  double deltat;                                        // change in time during superdrop model
  double deltat_cond;                                   // time difference since last condensation event
  double deltat_coll;                                   // time difference since last collision event

  /* variables required for superdroplet model (SDM) */
  Superdrop superdrops_array[init::NSUPERS]; // array of instances of Superdroplet class (superdroplet objects)
  const int nsupers = init::NSUPERS;       // no. of superdroplets to be included in model (nsupers <= init::NSUPERS)
  double deltemp, delqv, delqc;            // changes to kinematic variables due to SDM
  double p, temp, qv, qc;                  // copies of kinematic variables from ode sovler passed to SDM

  /* initial kinematic conditions for ODEs */
  const double w = init::iW / dlc::W0;          // dimensionless w velocity passed to kinetics ODEs eg. dp_dt(t,y,ydot,w,...)
  const double p_init = init::P_INIT / dlc::P0;
  const double temp_init = init::TEMP_INIT / dlc::TEMP0;
  const double vapourp_init = saturation_pressure(temp_init) * init::relh_init / 100; //initial vapour pressure for calculating qv_init
  const double qv_init = vapourpressure_2_massmixratio(vapourp_init, p_init);
  const double qc_init = init::qc_init;

  /* Initialise Superdroplets using initdrops_csv file */
  for (int i = 0; i < nsupers; i++)
  {
    superdrops_array[i] = Superdrop(init::iRho_l,
                                  init::iRho_sol, init::iMr_sol, init::iIONIC);
  }
  initialise_Superdrop_instances(init::initdrops_csv,
                                 superdrops_array, nsupers);

  /* Get nhalf and pvec (list of indicies) given constant no. of nsupers */
  const int nhalf = floor(nsupers / 2.0);
  vector<int> pvec(nsupers);
  for (int i = 0; i < nsupers; i++)
  {
    pvec[i] = i;
  }

  /* setup CVODE ODE solver */
  CvodeOdeSolver cvode;
  cvode.init_userdata(w, init::doThermo);
  const double y_init[4] = {p_init, temp_init, qv_init, qc_init};
  cvode.setup_ODE_solver(init::rtol, init::atols, y_init, t_init);
  cvode.print_init_ODEdata(nout, t_init, ceil(iters_per_tstep) * min_tstep, tstep);

  /* copy kinematic variables from CVODE ODE solver to SDM */
  double &t = cvode.t;
  cvode.get_variables_b4tstep(p, temp, qv, qc, deltemp, delqv, delqc);

  /* write setup to .txt file */
  WriteSetup2Txt(setupfiles2copy, init::setup_txt);

  /* write header and initial data to .csv files */
  write_outputheader(init::solution_csv);
  write_superdrop_outputheader(init::solutionSD_csv);
  kinetics_datafile.open(init::solution_csv, ios::app);
  superdrops_datafile.open(init::solutionSD_csv, ios::app);
  write_output(kinetics_datafile, t, p, temp, qv, qc);
  write_superdrop_output(superdrops_datafile, superdrops_array, nsupers);

  /* *** run superdroplet model (SDM) coupled to
    CVODE ODE SOLVER for kinematics. Writing output
    data nout no. of times within init::TSPAN) *** */
  tout = t_init + ceil(iters_per_tstep) * min_tstep; // first output time of ODE solver
  deltat_cond = 0;
  deltat_coll = 0;
  for (int j = 0; j < nout; j++)
  {
    print_output(t, p, temp, qv, qc); // print variables at timestep t

    /* SECTION 1: run SUPERDROPLET model */
    /* (a) divide each tstep into iters_per_tstep no. of min_tsteps and
      model SD collisions and/or condensational growth for each min_tstep */
    deltat = 0;
    for (int k = 0; k < ceil(iters_per_tstep); k++) // increment time of SDM
    {
      deltat += min_tstep;
      deltat_cond += min_tstep; // change in time since last condensation event
      deltat_coll += min_tstep; // change in time since last collision event

      /* (b) model superdroplet condensation-diffusion growth */
      if (init::doCond)
      {
        if (deltat_cond >= cond_tstep)
        {
          condensation_onto_superdroplets(cond_tstep, p, temp, qv, qc,
                                          deltemp, delqv, delqc, superdrops_array, nsupers);
          deltat_cond = 0; // reset time since last condensation event
        }
      }

      /* (c) model superdroplet collisions */
      if (init::doColl)
      {
        if (deltat_coll >= coll_tstep)
        {
          collide_droplets(nsupers, nhalf, pvec, superdrops_array);
          deltat_coll = 0; // reset time since last collision event
        }
      }
    }

    /* SECTION 2: run timestep of CVODE ODE solver */
    is_cvodefail = cvode.advance_solution(tout);
    if (is_cvodefail)
    {
      break;
    }

    /* reinitialise solver with changes to temp
      qv and qc due to SDM condensation */
    if (init::doCond)
    {
      is_cvodefail = cvode.reinitialise(tout, deltemp, delqv, delqc);
    }

    /* SECTION 3: write output data and proceed to next time step */
    cvode.get_variables_b4tstep(p, temp, qv, qc, deltemp, delqv, delqc);
    write_superdrop_output(superdrops_datafile, superdrops_array, nsupers);
    write_output(kinetics_datafile, t, p, temp, qv, qc);

    if (is_cvodefail)
    {
      break;
    }
    tout += deltat;
  }

  /* end CVODE ODE solver and close data files */
  superdrops_datafile.close();
  kinetics_datafile.close();
  cvode.destroy_cvode();

  return 0;
}