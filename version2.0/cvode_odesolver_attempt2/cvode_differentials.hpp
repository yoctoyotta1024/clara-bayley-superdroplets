// Author: Clara Bayley
// File: cvode_differentials.hpp
/* Header file for ODEs which are solved by
 CVODE ode solver to model evolution of the
 kinetics (p, temp, qv and qc) over time */


#ifndef CVODE_DIFFERENTIALS_HPP
#define CVODE_DIFFERENTIALS_HPP



// #include <iostream>

#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */

// using namespace std;



/* user data structure for passing 
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct {
  realtype w;
  bool doThermo;
} *UserData;


/* Simple function f(t,y, ydot) called by ODE solver to 
  solve differential equations over time. */
int odes_func(realtype t, N_Vector y, N_Vector ydot, void* user_data);








#endif //CVODE_DIFFERENTIALS_HPP