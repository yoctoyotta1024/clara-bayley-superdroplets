#ifndef JACOBIAN 
#define JACOBIAN 

#include "differentials.hpp"

/* Some Definitions */
/* User-defined matrix accessor macro: IJth 
  IJth(A,i,j) references the (i,j)th element (row, column) of the 
  dense matrix A, where i and j are in the range [1..NEQ]. 
*/
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1)    // (i,j)-th matrix component i,j=1..NEQ
#define ZERO  RCONST(0.0)

/* -------------------------------
 * Functions called by the solver for Jacobian of ODEs
 *------------------------------- */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
// ---------------------------------------------------------------------- //



/*
 * Jacobian routine. Compute J(t,y) = df/dy.*
 */

static int Jac(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int NEQ = NV_LENGTH_S(y);

  realtype p, temp, qv, qc, dp, dtemp;
  p = Ith(y,1); 
  temp = Ith(y,2);
  qv = 0.0;  //Ith(y,3);
  qc = 0.0;  //Ith(y,4);

  dp =  Ith(ydot,1);
  dtemp =  Ith(ydot,2);
  
  for (int j = 1; j <= NEQ; j++){
    IJth(J,1,j) = RCONST(0.0);                 //dp/dt if f(t) only so all Jacobian elements = 0
  }
  IJth(J,2,1) = -dtemp/p;
  IJth(J,2,2) = dtemp/temp;

  
  return(0);
}









#endif //JACOBIAN 