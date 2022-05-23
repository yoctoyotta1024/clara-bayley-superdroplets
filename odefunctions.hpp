#include <cmath>


/* -------------------------------
 * Functions called by the solver
 *------------------------------- */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

//static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
// ---------------------------------------------------------------------- //


/* Some Definitions */
#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */
#define ZERO  RCONST(0.0)


// /*
//  * g routine. Compute functions g_i(t,y) for i = 0,1.
//  */

// static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
// {
//   realtype y1, y3;

//   y1 = Ith(y,1); y3 = Ith(y,3);
//   gout[0] = y1 - RCONST(0.0001);
//   gout[1] = y3 - RCONST(0.01);

//   return(0);
// }



/*
 * Simple function f(t,y) that is the differential equations.
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype v, y1, y2, yd1, yd2;

  v = sin(pow(t, 0.2));
  y1 = Ith(y,1); 
  y2 = Ith(y,2);

  yd1 = Ith(ydot,1) = y2;
  yd2 = Ith(ydot,2) = v*y2*(1-pow(y1,2)) - y1;

  return(0);
}




/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype v, y1, y2;
  v = sin(pow(t, 0.2));
  y1 = Ith(y,1);
  y2 = Ith(y,2);
  
  IJth(J,1,1) = RCONST(0.0);
  IJth(J,1,2) = RCONST(1.0);

  IJth(J,2,1) = RCONST(-2.0)*v*y1*y2 - RCONST(1.0);
  IJth(J,2,2) = v*(RCONST(1.0)-pow(y1,2));
  

  return(0);
}
