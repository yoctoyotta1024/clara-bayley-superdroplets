// Author: Clara Bayley
// File: cvode_odesolver.hpp
/* Header file for CVODE ode solver 
which models evolution of the kinetics 
(p, temp, qv and qc) over time*/


#ifndef CVODE_ODESOLVER_HPP
#define CVODE_ODESOLVER_HPP 


#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */








#endif //CVODE_ODESOLVER_HPP  