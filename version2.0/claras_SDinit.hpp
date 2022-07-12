// Author: Clara Bayley
// File: claras_SDinit.hpp
/* File containing initial conditions
for CVODE ODE solver and SDM */

#ifndef CLARAS_SDINIT_HPP
#define CLARAS_SDINIT_HPP

#include <iostream>

#include "claras_SDconstants.hpp"
namespace dlc = dimless_constants;

namespace init
{

/* namespace containing paramaters for
initialisation of model */

/* filenames */
const std::string setup_txt = "sol_setup.txt";
const std::string initdrops_csv = "../dimlessSDinit.csv"; //.csv filename for initialisation of superdrops
const std::string solutionSD_csv = "solSD.csv";           //.csv filename for initialisation of superdrops
const std::string solution_csv = "sol.csv";               //.csv filename for initialisation of superdrops

/* Model Settings */
const bool doThermo = true; // enable condensational growth of superdroplets
const bool doCond = false;  // enable condensational growth of superdroplets
const bool doColl = true;   // enable collisions of superdroplets

/* initial parcel conditions */
const double TEMP_INIT = 273.15; // initial parcel temperature [T]
const double P_INIT = 100000;    // initial pressure [Pa]
const double relh_init = 60;     // initial relative humidity (%)
const double qc_init = 0;        // initial liquid water content []

/* droplet init params. Note: Superdroplet inital eps, r0
and m_sol data is made using python "createSDinitcsv.py" */
const int NSUPERS = 8192; // max. no. distinct superdrop objects in array
//const int NSUPERS = 2;
const double DROPVOL = 1e6; // volume of parcel occupied by superdroplets [m^3]
const double iRho_l = dlc::Rho_l;
const double iRho_sol = dlc::Rho_sol;
const double iMr_sol = dlc::Mr_sol;
const int iIONIC = dlc::IONIC;

/* Model timestep parameters */
const double COND_TSTEP = 1;       // time between SD condensation events  = ceil(coll/cond)*min(coll,cond) [s]
const double COLL_TSTEP = 1;       // time between SD collision events = ceil(coll/cond)*min(coll,cond) [s]
const int nout = 40;               // No. time points to write output data at
const double TSPAN[2] = {0, 4000}; // time span of integration [s]

/* CVODE ODE solver paramters */
const double iW = 0.5;                            // vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
const double rtol = 1e-6;                         // relative tolerance (tol) for integration
const double atols[4] = {1e-6, 1e-6, 1e-6, 1e-6}; // absolute tolerances for kinetics [P, T, qv, qc]

}

#endif // CLARAS_INIT_HPP