// Author: Clara Bayley
// File: cvode_odesolver.cpp
/* This file contains the CVODE ode solver 
for the evolution of the kinetics 
(p, temp, qv and qc) over time*/


#include "cvode_odesolver.hpp"


void cvode_odesolver{

  // Create the SUNDIALS solver stuff
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;                                           // reusable return flag

  // Create problem stuff
  realtype t, tout;
  N_Vector y; //, e;
  N_Vector abstol;

  // initialise vectors, matrix and solver
  UserData data;
  data = (UserData) malloc(sizeof *data);
  y = NULL;
  abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;


}