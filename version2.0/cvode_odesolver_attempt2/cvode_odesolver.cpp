// Author: Clara Bayley
// File: cvode_odesolver.cpp
/* This file contains the CVODE ode solver 
functionality for the evolution of the kinetics 
(p, temp, qv and qc) over time */


#include "cvode_odesolver.hpp"




CvodeOdeSolver::CvodeOdeSolver()
{
  /* initialise vectors, matrix and solver */
  data = (UserData) malloc(sizeof *data);
  y = NULL;
  ATOLS = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;
}





void CvodeOdeSolver::init_userdata(const double w, const bool doThermo)
{
  data->w = w;
  data->doThermo = doThermo;
};




int CvodeOdeSolver::setup_ODE_solver(const double i_rtol, const double i_atols[NEQ],
                      const double y_init[NEQ], const double t0)
{
 
  /* 0. Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)){ retval = 1;}

  /*  1. Initialize parallel or multi-threaded environment, if appropriate. */
  // ---------------------------------------------------------------------- //
  
  /* 2. Set the scalar relative and vector absolute tolerances */
  RTOL = i_rtol;
  ATOLS = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)ATOLS, "N_VNew_Serial", 0)) return(1);

  for(int i=0; i<NEQ; i++)
  {
    NV_Ith_S(ATOLS,i) = i_atols[i];
  }

  /* 3. initialise y vector with initial conditions */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
  for(int i=0; i<NEQ; i++)
  {
    NV_Ith_S(y,i) = y_init[i];
  }
 
  /* 4. Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula (CV_BDF) */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* 5. Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, odes_func, t0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);
  
  /* 6. Set linear solver optional inputs. 
   * Provide user data which can be accessed in user provided routines */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval((void *)&retval, "CVodeSetUserData", 1)) return(1);
  
  /* 7. Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, RTOL, ATOLS);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  /* 8. Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* 9. Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);
  
  /* 10. Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);
 
  return 0;
};




int CvodeOdeSolver::advance_solution(const double tout)
/* Advance ODE solution in time */
{

  retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1)) return 1;
  //cout << "t of ODE SOLVER -> " << tout << endl;

  return 0;
}



// int CvodeOdeSolver::reinitialise(const double tout, const double re_y[NEQ])
// /* Reinitialize the solver after discontinuous change in y */
// {
//   retval = CVodeReInit(cvode_mem, tout, re_y);
//   if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);

//   return 0;
// }




/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */
int CvodeOdeSolver::check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n"
            << funcname << " " << endl;
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed with retval = %d\n\n"
             << funcname << " " << *retval << endl;
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    cout << stderr << "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n"
            << funcname << endl;
    return(1); }

  return(0);
}



int CvodeOdeSolver::print_init_ODEdata(const int nout, const double t0,
        const double t1, const double tstep)
{

  if (y == NULL)
  { 
    retval = -1;
    if(check_retval(&retval, "print_init_ODEdata", 1)) return(1);
    return retval;
  }

  cout << "---------------------------------------------" << endl;
  cout << "No. Equations (NEQ) = " << NEQ << endl;
  cout << "y0      = " << NV_Ith_S(y,0) << endl;
  cout << "y1      = " << NV_Ith_S(y, 1) << endl;
  cout << "y2      = " << NV_Ith_S(y, 2) << endl;
  cout << "y3      = " << NV_Ith_S(y, 3) << endl;
  cout << "---------------------------------------------" << endl;
  cout << "RTOL        = " << RTOL << endl;
  cout << "ATOLS[0]    = " << NV_Ith_S(ATOLS, 0) << endl;
  cout << "ATOLS[1]    = " << NV_Ith_S(ATOLS, 1) << endl;
  cout << "ATOLS[2]    = " << NV_Ith_S(ATOLS, 2) << endl;
  cout << "ATOLS[3]    = " << NV_Ith_S(ATOLS, 3) << endl;
  cout << "---------------------------------------------" << endl;
  cout << "inital t    = " << t0 << endl;
  cout << "first t     = " << t1 << endl;
  cout << "t step      = " << tstep << endl;
  cout << "no. outputs = " << nout << endl;
  cout << "---------------------------------------------\n" << endl;


  retval = 0;
  if(check_retval(&retval, "print_init_ODEdata", 1)) return(1);
  return retval;
};




void CvodeOdeSolver::destroy_cvode()
{

  /* print final statistics to the terminal screen */
  printf("\nLast Iteration Statistics:\n");
  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  /* free memory */
  N_VDestroy(y);                            /* Free y vector */
  N_VDestroy(ATOLS);                       /* Free abstol vector */
  free(data);                               /* free user_data pointer struc */
  CVodeFree(&cvode_mem);                    /* Free CVODE memory */
  SUNLinSolFree(LS);                        /* Free the linear solver memory */
  SUNMatDestroy(A);                         /* Free the matrix memory */
  SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */

}

