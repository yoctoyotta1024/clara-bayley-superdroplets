#include <iostream>
#include <stdio.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

#include "init.hpp"
#include "constants.hpp"
#include "differentials.hpp"
// #include "jacobian.hpp"               // file with Jacobian function (is using optional step 12.)
#include "cvodehelpers.hpp"



namespace dlc = dimless_constants;
using namespace dlc;


/* User-defined vector and matrix accessor macros: Ith, IJth 

  Ith(v,i) references the ith component of the vector v, where i is in
  the range [1..NEQ] (and NEQ is No. of equations ie. variables). 

  IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
  i and j are in the range [1..NEQ]. 
*/
#define Ith(v,i)    NV_Ith_S(v,i-1)         // i-th vector component i=1..NEQ
//#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) // (i,j)-th matrix component i,j=1..NEQ







/* Problem Constants */
#define NEQ   2                                  // number of equations
#define Y1    RCONST(p_init/dlc::P0)             // initial y components
#define Y2    RCONST(temp_init/dlc::TEMP0)

#define RTOL  RCONST(rtol)                       // scalar relative tolerance
#define ATOL1 RCONST(atols[0])                   // vector absolute tolerance components
#define ATOL2 RCONST(atols[0])

#define T0    RCONST(tspan[0]/dlc::TIME0)        // initial time (dimensionless)          
#define TSTEP RCONST(tspan[1]/nout/dlc::TIME0)   // output time step (dimensionless)     
#define NOUT  nout                               // number of output times
#define ZERO  RCONST(0.0)




/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(){


  // TODO de-dimesionalise: w = w/W0
  // TODO de-dimesionalise: z/(W0*T0))


  /* Project name and savefile names*/
  std::string PROJNAME, SAVENAME, STATSNAME, savey, saverr;
  PROJNAME = "2-species dynamics problem";
  SAVENAME = "sundials2";
  STATSNAME = SAVENAME+"_stats.csv";


  // Create the SUNDIALS solver stuff//
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;           /* reusable return flag */
  int retvalr;
  int rootsfound[2];

  // Create problem stuff
  realtype t, tout;
  N_Vector y, e;
  N_Vector abstol;

  // Output files to write to
  FILE* SFID;                // integration stats output file 
  FILE *YFID = NULL;         // solution output file 
  FILE *EFID = NULL;         // error output file   

  // initialise vectors, matrix and solver
  y = NULL;
  abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;


  /* 0. Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /*  1. Initialize parallel or multi-threaded environment, if appropriate. */
  // ---------------------------------------------------------------------- //
    
  /* 2. Initial conditions */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Create serial vector to store the solution error */
  e = N_VClone(y);
  if (check_retval((void *)y, "N_VClone", 0)) return(1);
  N_VConst(ZERO, e);                           /* Set initial error vector */

  /* Initialize y vector */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  // Ith(y,3) = Y3;

  /* 3. Set the vector absolute tolerance */
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

  Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;

  /* 4. Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* 5. Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* 6. Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, RTOL, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  /* 7. Optional Inputs:
  Call CVodeRootInit to specify the root function g with 2 components */
  //retval = CVodeRootInit(cvode_mem, 2, g);
  //if (check_retval(&retval, "CVodeRootInit", 1)) return(1);

  /* 8. Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* 9. Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);
  
  // 10. Set linear solver optional inputs.
  // ---------------------------------------------------------------------------

  /* 11. Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  // /* 12. Set the user-supplied Jacobian function Jac (optional) */
  // retval = CVodeSetJacFn(cvode_mem, Jac);
  // if (check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* 13. Print problem setup and write initial conditions*/
  retval = PrintINITData(PROJNAME, NEQ, y, 
                  RTOL, abstol, NOUT, T0, TSTEP, TSTEP);
  if(check_retval(&retval, "PrintINITData", 1)) return(1);

  /* Output initial conditions */
  savey = SAVENAME+"_sol.csv";
  saverr = SAVENAME+"_err.csv";
  YFID = fopen(savey.c_str(),"w");
  EFID = fopen(saverr.c_str(),"w");
  WriteOutput(1, YFID, EFID);
  WriteOutput(T0, y, e, 1, YFID, EFID);
  
  /* Open Integration Statistics File in preparation for writing */ 
  SFID = fopen(STATSNAME.c_str(), "w");

  /* 14. RUN SOLVER
  For NOUT iterations, call CVode, test for error then
  print and write results and iterate (if sucess) */
  
  tout = T0 + TSTEP;   // first output time = t0 + TSTEP
  for (int iout = 0; iout < NOUT; iout++){

    PrintOutput(t, y);
    
    /* 14(a) Advance solution in time */
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;

    if (retval == CV_ROOT_RETURN) {
      retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
      if (check_retval(&retvalr, "CVodeGetRootInfo", 1)) return(1);
      PrintRootInfo(rootsfound[0],rootsfound[1]);
    }

    /* 14(b) Output solution and error */
    //retval = ComputeError(t, y, e, &ec, udata);
    //if (check_retval(&retval, "ComputeError", 1)) break;
    WriteOutput(t, y, e, 1, YFID, EFID);
    if (check_retval(&retval, "WriteOutput", 1)) break;

    /* 14(c) Continute to next timestep */
    if (retval == CV_SUCCESS) {
      tout += TSTEP;
    }

    retval = CVodePrintAllStats(cvode_mem, SFID, SUN_OUTPUTFORMAT_CSV);
 
  }

  fclose(YFID);
  fclose(EFID);
  fclose(SFID);


  /* 15. Print final statistics to the screen */
  printf("\nLast Iteration Statistics:\n");
  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);


  /* 16. Free memory */
  N_VDestroy(y);                            /* Free y vector */
  N_VDestroy(e);                            /* Free error vector */
  N_VDestroy(abstol);                       /* Free abstol vector */
  CVodeFree(&cvode_mem);                    /* Free CVODE memory */
  SUNLinSolFree(LS);                        /* Free the linear solver memory */
  SUNMatDestroy(A);                         /* Free the matrix memory */
  SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */

  return(0);
}

