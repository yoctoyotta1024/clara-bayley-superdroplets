#include <iostream>
#include <stdio.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

#include "init.hpp"
#include "constants.hpp"
#include "superdroplets.hpp"
#include "differentials.hpp"
#include "collisions.hpp"
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

#define SDloop(i,nsupers) for(int i=0; i<nsupers; i++)  //for loop over all superdroplets





/* ODE Solving Constants */
#define NEQ   4                                  // number of equations
#define Y1    RCONST(p_init/dlc::P0)             // initial y components
#define Y2    RCONST(temp_init/dlc::TEMP0)
#define Y4    RCONST(qc_init)
const realtype  w = iW/dlc::W0;                  // dimensionless w velocity for f(t,y,ydot,w,...)
//const realtype Vol = iVOL/dlc::VOL0;             // dimensionless volume of parcel

#define ZERO  RCONST(0.0)
#define RTOL  RCONST(rtol)                       // scalar relative tolerance
#define ATOLS RCONST(atols)                     // vector absolute tolerance components for y1 to y4

#define T0    RCONST(tspan[0]/dlc::TIME0)        // initial time (dimensionless)          
#define TSTEP RCONST(tspan[1]/nout/dlc::TIME0)   // output time step (dimensionless)     
#define NOUT  nout                               // number of output times

#define COND_TSTEP RCONST(cond_tstep/dlc::TIME0)                    // No. droplet collisions events
#define COLL_TSTEP RCONST(coll_tstep/dlc::TIME0)                    // No. droplet collisions events


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main(){

  /* Project name and savefile names*/
  std::string PROJNAME, SAVENAME, STATSNAME;
  PROJNAME = "Simple Rising Parcel (with Drop Condensation tbc)";          // project title sometimes used
  SAVENAME = "./bin/sundials2";                                                  // solution written to SAVENAME_sol.csv file
  STATSNAME = SAVENAME+"_stats.csv";                                       // integration statistics written to this .csv

  string INITDROPSCSV;
  INITDROPSCSV = "dimlessSDinit.csv";                                       // file to read for SD eps, r and m_sol initialisation 

  // ------------------ ODE Solver stuff beyond this line ------------------- //
  
  // Create the SUNDIALS solver stuff//
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;          // reusable return flags
  //int retvalr;           
  //int rootsfound[2];

  // Create problem stuff
  realtype t, tout;
  N_Vector y; //, e;
  N_Vector abstol;
  realtype ItersPerTstep, min_tstep; 
  realtype delta_t, dt_cond, dt_coll;

  // Output files to write to
  FILE* STSFID;                // integration stats output file 
  FILE *IFID = NULL;           // setup (init.hpp and consts.hpp) written to this file 
  FILE *YFID = NULL;           // solution output file 
  FILE *SDFID = NULL;          // solution output file 
  //FILE *EFID = NULL;         // error output file   

  // initialise vectors, matrix and solver
  UserData data;
  data = (UserData) malloc(sizeof *data);
  y = NULL;
  abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Initial qv and qc from relative humidity */
  realtype pv_init;
  pv_init = saturation_pressure(Y2)*relh_init/100;
  realtype Y3 = pv2qv(pv_init, Y1);      //initial qv for solver

  /* Initialise Superdroplets using INITDROPSCSV .csv file */
  Superdrop drops_arr[nsupers];
  SDloop(i, nsupers)
  {
    drops_arr[i] = Superdrop(iRho_l, iRho_sol, iMr_sol, iIONIC); 
  }
  Superdrop* ptr; 
  ptr = &drops_arr[0];
  initialise_Superdrop_instances(INITDROPSCSV, ptr, nsupers);
  
  /* set values of pointers given in user_data to f() ODE function */
  InitUserData(data, w, doCond, nsupers, ptr);

  /* Get nhalf, scale_p and pvec (index list) given nsupers */
  int nhalf = floor(nsupers/2);
  int scale_p = nsupers*(nsupers-1)/(2*nhalf);
  vector<int> pvec(nsupers);
  SDloop(i, nsupers){
    pvec[i] = i;
  }

  /* 0. Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /*  1. Initialize parallel or multi-threaded environment, if appropriate. */
  // ---------------------------------------------------------------------- //

  /* 2. Initial conditions */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  // /* Create serial vector to store the solution error */
  // e = N_VClone(y);
  // if (check_retval((void *)y, "N_VClone", 0)) return(1);
  // N_VConst(ZERO, e);                           /* Set initial error vector */

  /* Initialize y vector */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;
  Ith(y,4) = Y4;

  /* 3. Set the vector absolute tolerance */
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return(1);

  for(int i=1; i<5; i++)
  {
    Ith(abstol,i) = ATOLS[i-1];
  }

  /* 4. Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* 5. Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);
  
  /* 10. Set linear solver optional inputs. 
      Provide user data which can be accessed in user provided routines */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval((void *)&retval, "CVodeSetUserData", 1)) return(1);
  
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

  /* Write Setup (init.hpp and constants.hpp) to IFID file */ 
  string setupfiles[] = {"init.hpp", "constants.hpp"}; 
  IFID = fopen((SAVENAME+"_setup.txt").c_str(),"w");
  WriteSetup2Txt(&setupfiles[0], 2, IFID);
  fclose(IFID);

  /* Output initial conditions */
  YFID = fopen((SAVENAME+"_sol.csv").c_str(),"w");
  SDFID = fopen((SAVENAME+"_SDsol.csv").c_str(),"w");
  //EFID = fopen((SAVENAME+"_err.csv".c_str(),"w"); 
  HeaderWriteOutput(YFID, SDFID);   
  WriteOutput(T0, y, nsupers, ptr, YFID, SDFID, NULL); 

  /* Open Integration Statistics File in preparation for writing */ 
  STSFID = fopen(STATSNAME.c_str(), "w");


  /* 14. RUN SOLVER
  For NOUT iterations, call CVode, test for error then
  print and write results and iterate (if sucess) */
  min_tstep = min(COLL_TSTEP, COND_TSTEP);         // smallest timestep is one to increment time by
  ItersPerTstep = TSTEP/min_tstep;                 // no. of increments = ceil(TSTEP/min_tstep) >= 1
  tout = T0+ceil(ItersPerTstep)*min_tstep;         // first output time of ODE solver

  realtype dtemp, dqv, dqc; 
  for (int iout = 0; iout < NOUT; iout++){

    PrintOutput(t, y);

    /* Run SD Model */
    //cout << " -- t of SDs: " << t << endl;
    delta_t = 0;
    dtemp = 0;
    dqv = 0;
    dqc = 0; 
    for(int istep=0; istep<ceil(ItersPerTstep); istep++)                 // increment time for SDs simulation
    {
      delta_t += min_tstep;
      dt_cond += min_tstep;  // change in time since last condensation event
      dt_coll += min_tstep;  // change in time since last collision event

      /* 14(b) Simulate Superdroplet Diffusion Growth */
      if(doCond){
        if(dt_cond >= COND_TSTEP)
        {
          //cout << "cond @ " << t+delta_t << endl;
          condensation_droplet_growth(COND_TSTEP, &Ith(y,1), &Ith(y,2), 
                      &Ith(y,3), &Ith(y,4), &dtemp, &dqv, &dqc, ptr, nsupers);
          dt_cond=0;
        }
      }
    
      /* 14(c) Simulate Superdroplet Collisions */
      if (doColl){
        if(dt_coll >= COLL_TSTEP)
        {
          //cout << "coll @ " << t+delta_t << endl;
          collide_droplets(nsupers, nhalf, scale_p, ptr, pvec);
          dt_coll=0;
        }
      }
        
    }
    //cout << "t of SDs -> " << t+delta_t << endl;

    /* 14(a) Advance ODE solution in time */
    //cout << " -- t of ODE SOLVER: " << t << endl;
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;
    //cout << "t of ODE SOLVER -> " << tout << endl;

    // if (retval == CV_ROOT_RETURN) {
    //   retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
    //   if (check_retval(&retvalr, "CVodeGetRootInfo", 1)) return(1);
    //   PrintRootInfo(rootsfound[0],rootsfound[1]);
    // }
     
    /* 14(d) Continute to next timestep */
    if(doCond){
      Ith(y,2) += dtemp; 
      Ith(y,3) += dqv;
      Ith(y,4) += dqc;
      /* Reinitialize the solver after discontinuous change in y */
      retval = CVodeReInit(cvode_mem, tout, y);
      if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);
    }

    if (retval == CV_SUCCESS) {
      tout += delta_t;
    }
    
    /* 14(e) Output solution and error after every large timestep */
    //retval = ComputeError(t, y, e, &ec, udata);
    //if (check_retval(&retval, "ComputeError", 1)) break;
    WriteOutput(t, y, nsupers, ptr, YFID, SDFID, NULL);
    if (check_retval(&retval, "WriteOutput", 1)) break;

    retval = CVodePrintAllStats(cvode_mem, STSFID, SUN_OUTPUTFORMAT_CSV);
 
  }

  fclose(YFID);
  fclose(SDFID);
  //fclose(EFID);
  fclose(STSFID);


  /* 15. Print final statistics to the screen */
  printf("\nLast Iteration Statistics:\n");
  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);


  /* 16. Free memory */
  N_VDestroy(y);                            /* Free y vector */
  // N_VDestroy(e);                         /* Free error vector */
  N_VDestroy(abstol);                       /* Free abstol vector */
  free(data);                               /* free user_data pointer struc */
  CVodeFree(&cvode_mem);                    /* Free CVODE memory */
  SUNLinSolFree(LS);                        /* Free the linear solver memory */
  SUNMatDestroy(A);                         /* Free the matrix memory */
  SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */

  return(0);
}

