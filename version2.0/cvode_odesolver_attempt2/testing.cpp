#include <iostream>
//#include <stdio.h>
#include <math.h>

#include "../claras_SDinit.hpp"
#include "../claras_SDconstants.hpp"
#include "cvode_odesolver.hpp"


namespace dlc = dimless_constants;



// // --------------------------------- cvode header ------------------------- //
// //#include <iostream>

// #include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
// #include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
// #include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
// #include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

// //#include "differentials.hpp"
// #include "cvode_differentials.hpp"

// using namespace std;


// class CvodeOdeSolver
// /* CVODE ODE Solver class */
// {
//   private:

//   public:
//     /* SUNDIALS solver stuff */
//     SUNContext sunctx;
//     SUNMatrix A;
//     SUNLinearSolver LS;
//     void *cvode_mem;
//     int retval;                                           // reusable return flag
    
//     /* ODE problem stuff */
//     static const int NEQ = 4;
//     UserData data;
//     realtype t;
//     N_Vector y;
//     N_Vector ATOLS; 
//     realtype RTOL;

//     /* constructor functions */
//     CvodeOdeSolver();
        
//     int setup_ODE_solver(const double i_rtol, const double i_atols[NEQ],
//                       const double y_init[NEQ], const double t0);

//     void init_userdata(const double w, const bool doThermo);

//     int advance_solution(const double tout);

//     int reinitialise(const double tout, const double re_y[NEQ]);

//     int check_retval(void *returnvalue, const char *funcname, int opt);
    
//     int print_init_ODEdata(const int nout, const double t0,
//         const double t1, const double tstep);

//     void destroy_cvode(); 
// };



// // -------------------------------------------------------------------- //


// // --------------------------------- cvode cpp ------------------------- //




// CvodeOdeSolver::CvodeOdeSolver()
// {
//   /* initialise vectors, matrix and solver */
//   data = (UserData) malloc(sizeof *data);
//   y = NULL;
//   ATOLS = NULL;
//   A = NULL;
//   LS = NULL;
//   cvode_mem = NULL;
// }





// void CvodeOdeSolver::init_userdata(const double w, const bool doThermo)
// {
//   data->w = w;
//   data->doThermo = doThermo;
// };




// int CvodeOdeSolver::setup_ODE_solver(const double i_rtol, const double i_atols[NEQ],
//                       const double y_init[NEQ], const double t0)
// {
 
//   /* 0. Create the SUNDIALS context */
//   retval = SUNContext_Create(NULL, &sunctx);
//   if (check_retval(&retval, "SUNContext_Create", 1)){ retval = 1;}

//   /*  1. Initialize parallel or multi-threaded environment, if appropriate. */
//   // ---------------------------------------------------------------------- //
  
//   /* 2. Set the scalar relative and vector absolute tolerances */
//   RTOL = i_rtol;
//   ATOLS = N_VNew_Serial(NEQ, sunctx);
//   if (check_retval((void *)ATOLS, "N_VNew_Serial", 0)) return(1);

//   for(int i=0; i<NEQ; i++)
//   {
//     NV_Ith_S(ATOLS,i) = i_atols[i];
//   }

//   /* 3. initialise y vector with initial conditions */
//   y = N_VNew_Serial(NEQ, sunctx);
//   if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
//   for(int i=0; i<NEQ; i++)
//   {
//     NV_Ith_S(y,i) = y_init[i];
//   }
 
//   /* 4. Call CVodeCreate to create the solver memory and specify the
//    * Backward Differentiation Formula (CV_BDF) */
//   cvode_mem = CVodeCreate(CV_BDF, sunctx);
//   if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

//   /* 5. Call CVodeInit to initialize the integrator memory and specify the
//    * user's right hand side function in y'=f(t,y), the initial time T0, and
//    * the initial dependent variable vector y. */
//   retval = CVodeInit(cvode_mem, odes_func, t0, y);
//   if (check_retval(&retval, "CVodeInit", 1)) return(1);
  
//   /* 6. Set linear solver optional inputs. 
//    * Provide user data which can be accessed in user provided routines */
//   retval = CVodeSetUserData(cvode_mem, data);
//   if (check_retval((void *)&retval, "CVodeSetUserData", 1)) return(1);
  
//   /* 7. Call CVodeSVtolerances to specify the scalar relative tolerance
//    * and vector absolute tolerances */
//   retval = CVodeSVtolerances(cvode_mem, RTOL, ATOLS);
//   if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

//   /* 8. Create dense SUNMatrix for use in linear solves */
//   A = SUNDenseMatrix(NEQ, NEQ, sunctx);
//   if (check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

//   /* 9. Create dense SUNLinearSolver object for use by CVode */
//   LS = SUNLinSol_Dense(y, A, sunctx);
//   if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);
  
//   /* 10. Attach the matrix and linear solver to CVODE */
//   retval = CVodeSetLinearSolver(cvode_mem, LS, A);
//   if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);
 
//   return 0;
// };




// int CvodeOdeSolver::advance_solution(const double tout)
// /* Advance ODE solution in time */
// {

//   //cout << " -- t of ODE SOLVER: " << t << endl;
//   retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
//   if (check_retval(&retval, "CVode", 1)) return 1;
//   //cout << "t of ODE SOLVER -> " << tout << endl;

//   return 0;
// }



// // int CvodeOdeSolver::reinitialise(const double tout, const double re_y[NEQ])
// // /* Reinitialize the solver after discontinuous change in y */
// // {
// //   retval = CVodeReInit(cvode_mem, tout, re_y);
// //   if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);

// //   return 0;
// // }




// /*
//  * Check function return value...
//  *   opt == 0 means SUNDIALS function allocates memory so check if
//  *            returned NULL pointer
//  *   opt == 1 means SUNDIALS function returns an integer value so check if
//  *            retval < 0
//  *   opt == 2 means function allocates memory so check if returned
//  *            NULL pointer
//  */
// int CvodeOdeSolver::check_retval(void *returnvalue, const char *funcname, int opt)
// {
//   int *retval;

//   /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
//   if (opt == 0 && returnvalue == NULL) {
//     cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n"
//             << funcname << " " << endl;
//     return(1); }

//   /* Check if retval < 0 */
//   else if (opt == 1) {
//     retval = (int *) returnvalue;
//     if (*retval < 0) {
//       cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed with retval = %d\n\n"
//              << funcname << " " << *retval << endl;
//       return(1); }}

//   /* Check if function returned NULL pointer - no memory allocated */
//   else if (opt == 2 && returnvalue == NULL) {
//     cout << stderr << "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n"
//             << funcname << endl;
//     return(1); }

//   return(0);
// }



// int CvodeOdeSolver::print_init_ODEdata(const int nout, const double t0,
//         const double t1, const double tstep)
// {

//   if (y == NULL)
//   { 
//     retval = -1;
//     if(check_retval(&retval, "print_init_ODEdata", 1)) return(1);
//     return retval;
//   }

//   cout << "---------------------------------------------" << endl;
//   cout << "No. Equations (NEQ) = " << NEQ << endl;
//   cout << "y0      = " << NV_Ith_S(y,0) << endl;
//   cout << "y1      = " << NV_Ith_S(y, 1) << endl;
//   cout << "y2      = " << NV_Ith_S(y, 2) << endl;
//   cout << "y3      = " << NV_Ith_S(y, 3) << endl;
//   cout << "---------------------------------------------" << endl;
//   cout << "RTOL        = " << RTOL << endl;
//   cout << "ATOLS[0]    = " << NV_Ith_S(ATOLS, 0) << endl;
//   cout << "ATOLS[1]    = " << NV_Ith_S(ATOLS, 1) << endl;
//   cout << "ATOLS[2]    = " << NV_Ith_S(ATOLS, 2) << endl;
//   cout << "ATOLS[3]    = " << NV_Ith_S(ATOLS, 3) << endl;
//   cout << "---------------------------------------------" << endl;
//   cout << "inital t    = " << t0 << endl;
//   cout << "first t     = " << t1 << endl;
//   cout << "t step      = " << tstep << endl;
//   cout << "no. outputs = " << nout << endl;
//   cout << "---------------------------------------------\n" << endl;


//   retval = 0;
//   if(check_retval(&retval, "print_init_ODEdata", 1)) return(1);
//   return retval;
// };




// void CvodeOdeSolver::destroy_cvode()
// {

//   /* print final statistics to the terminal screen */
//   printf("\nLast Iteration Statistics:\n");
//   retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

//   /* free memory */
//   N_VDestroy(y);                            /* Free y vector */
//   N_VDestroy(ATOLS);                       /* Free abstol vector */
//   free(data);                               /* free user_data pointer struc */
//   CVodeFree(&cvode_mem);                    /* Free CVODE memory */
//   SUNLinSolFree(LS);                        /* Free the linear solver memory */
//   SUNMatDestroy(A);                         /* Free the matrix memory */
//   SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */

// }




// // ---------------------------------------------------------------- //








// --------------------------------- main ------------------------- //


double pv2qv(const double pv, const double p)
	/* Calculate mass mixing ratio
	qv = m_v/m_dry = rho_v/rho_dry
	given vapour pressure pv = p_v/p_tot. */
{  
  return dlc::Mr_ratio * pv/(p-pv);
}



double saturation_pressure(const double temp)
	/* Calculate the equilibrium vapor pressure 
	of water over liquid water ie. the
	saturation pressure (psat). Equation taken from
	python module typhon.physics.thermodynamics.e_eq_water_mk
	with conversion to real temp /K = T*Temp0 and from 
	real psat to dimensionless psat = psat/P0. */
{
	double T, lnpsat;
  T = temp*dlc::TEMP0;                               // real T [K]
  
  if(T <= 0)
  {
      std::cout << "psat ERROR: T must be larger than 0K. T = " <<std::endl;
  }

  lnpsat = (54.842763                               // ln(psat) [Pa]
        - 6763.22 / T
        - 4.21 * log(T)
        + 0.000367 * T
        + tanh(0.0415 * (T - 218.8))
        * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T));

  return exp(lnpsat)/dlc::P0;                      // dimensionless psat
}



/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main(){

  int cvode_iterfail;                                        // flag to break integration if cvode fails
 
  const int nout = init::nout;                                 // number of output times
  const double t0 = init::TSPAN[0]/dlc::TIME0;                       // initial time (dimensionless)          
  const double tstep = init::TSPAN[1]/nout/dlc::TIME0;               // output time step (dimensionless)     

  const double cond_tstep = init::COND_TSTEP/dlc::TIME0;                   
  const double coll_tstep = init::COLL_TSTEP/dlc::TIME0;                    
  const double min_tstep = min(coll_tstep, cond_tstep);         // smallest timestep is one to increment time by
  const double ItersPerTstep = tstep/min_tstep;                 // no. of increments = ceil(TSTEP/min_tstep) >= 1
  double tout, delta_t, dt_cond, dt_coll;
  
  /* initial conditions & setup of model */  
  const double w = init::iW/dlc::W0;                          // dimensionless w velocity for f(t,y,ydot,w,...)
  const double p_init = init::P_INIT/dlc::P0;                 // initial (dimensionless) kinetics
  const double temp_init = init::TEMP_INIT/dlc::TEMP0;
  const double pv_i = saturation_pressure(temp_init)*init::relh_init/100;
  const double qv_init = pv2qv(pv_i, p_init);                 // Initial qv from relative humidity
  const double qc_init = init::qc_init;


  /* SETUP CVODE ODE SOLVER */
  CvodeOdeSolver cvode;

  cvode.init_userdata(w, init::doThermo);
  const double y_init[4] = {p_init, temp_init, qv_init, qc_init};
  cvode.setup_ODE_solver(init::rtol, init::atols, y_init, t0);
  
  
  /* START CVODE ODE SOLVER */
  //tout = T0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver
  cvode.print_init_ODEdata(nout, t0, ceil(ItersPerTstep)*min_tstep, tstep);


  /* RUN CVODE ODE SOLVER FOR KINEMATICS
    (collecting data nout no. of times within init::TSPAN) */
  tout = t0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver
  dt_cond = 0;
  dt_coll = 0;
  for (int j=0; j<nout; j++)
  {

    /* a) divide each tstep into ItersPerTstep no. of 
      min_tsteps and model SD collisions and/or condensational growth */
    delta_t = 0;
    for(int k=0; k<ceil(ItersPerTstep); k++)                 // increment time for SDs simulation
    {
      delta_t += min_tstep;        
    }

    cvode_iterfail = cvode.advance_solution(tout);
    if(cvode_iterfail){ break; }
    
    /* Continute to next timestep */ 
    // if(init::doCond)
    // {
    //   yvec[1] += deltemp;
    //   yvec[2] += delqv;
    //   yvec[3] += delqc; 

    //   cvode_iterfail = cvode.reinitialise(tout, y);
    // }

    if(cvode_iterfail){ break; }
    tout += delta_t;


  }

  /* END CVODE ODE SOLVER */
  cvode.destroy_cvode();




 return(0);
}

