// Author: Clara Bayley
// File: differential_functions.hpp
/* Header file for ODEs which are solved by
 CVODE ode solver to model evolution of the
 kinetics (p, temp, qv and qc) over time */

#ifndef DIFFERENTIAL_FUNCTIONS_HPP
#define DIFFERENTIAL_FUNCTIONS_HPP

// #include <iostream>
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */

#include "../../claras_SDconstants.hpp"

// using namespace std;
namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/* user data structure for passing
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct
{
  realtype w;
  bool doThermo;
} * UserData;

int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
/* Simple function f(t,y, ydot) called by ODE solver to
  solve differential equations over time. */

static double dp_dt(const double &t, const double &w);
/* dp/dt differential equation (dimensionless)
  describing pressure evolution over time.
  note: true dP/dt = dp/dt * P0/TIME0 */

static double cp_moist(const double &qv, const double &qc);
/* effective specific heat capacity of moist parcel
  of air (dry + water vapour + liquid water) */

static double dtemp_dt_adia(const double &pdot, const N_Vector &y);
/* dtemp/dt differential equation describing
  temperature evolution solely due to pressure
  changes in parcel for adiabatic process (no heat loss).
  Parcel has water vapour mass mixing ratio (m_v/m_dry) = qv and
  liquid water mass mixing ratio (m_c/m_dry) = qc.
  note: True dTemp/dt = dtemp * TEMP0/TIME0  */

#endif // DIFFERENTIAL_FUNCTIONS_HPP