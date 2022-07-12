// Author: Clara Bayley
// File: cvode_odesolver.hpp
/* Header file for CVODE ode solver
which models evolution of the kinetics
(p, temp, qv and qc) over time*/

#ifndef CVODE_ODESOLVER_HPP
#define CVODE_ODESOLVER_HPP

#include <iostream>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

#include "differential_functions.hpp"

using namespace std;

class CvodeOdeSolver
/* CVODE ODE Solver class */
{
private:
public:
  /* SUNDIALS solver stuff */
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval; // reusable return flag

  /* ODE problem stuff */
  static const int NEQ = 4;
  UserData data;
  realtype t;
  N_Vector y;
  N_Vector re_y;
  N_Vector ATOLS;
  realtype RTOL;

  /* constructor functions */
  CvodeOdeSolver();

  int setup_ODE_solver(const double i_rtol, const double i_atols[NEQ],
                       const double y_init[NEQ], const double t0);

  void init_userdata(const double w, const bool doThermo);

  void get_variables_b4tstep(double &p, double &temp, double &qv,
                             double &qc, double &deltemp, double &delqv, double &delqc);

  int advance_solution(const double tout);

  int reinitialise(const double tout, const double deltemp,
                   const double delqv, const double delqc);

  int check_retval(void *returnvalue, const char *funcname, int opt);

  int print_init_ODEdata(const int nout, const double t0,
                         const double t1, const double tstep);

  void destroy_cvode();
};

#endif // CVODE_ODESOLVER_HPP