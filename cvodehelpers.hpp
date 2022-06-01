#ifndef CVODEHELPERS
#define CVODEHELPERS

#include <iostream>
using namespace std;

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */
/* Private functions to output results */
static int PrintINITData(string PROJNAME, int NEQ, N_Vector y, realtype rtol,
          N_Vector abstol, int NOUT, realtype T0, realtype T1, realtype TSTEP);
static void PrintOutput(realtype t, N_Vector y);
static void PrintRootInfo(int root_f1, int root_f2);
static int HeaderWriteOutput(FILE* YFID, int index0, int indexend, bool fPrintTimeHeader);
static int WriteOutput(realtype t, N_Vector y, int nsupers, 
              FILE* YFID, FILE* SDFID, FILE* EFID);

/* Private function to check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);
// ---------------------------------------------------------------------- //



/* Some Definitions */
#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */



static int PrintINITData(string PROJNAME, int NEQ, N_Vector y, realtype rtol,
          N_Vector abstol, int NOUT, realtype T0, realtype T1, realtype TSTEP)
{

  if (y == NULL) return(-1);

  printf("\n%s\n\n", PROJNAME.c_str());
  printf("---------------------------------------------\n");
  printf("No. Equations (NEQ) = %d\n", NEQ);
  printf("y1      = %.4f\n", Ith(y,1));
  printf("y2      = %.4f\n", Ith(y, 2));
  printf("---------------------------------------------\n");
  printf("rtol       = %.3g\n", rtol);
  printf("abstol1    = %.3g\n", Ith(abstol, 1));
  printf("abstol2    = %.3g\n", Ith(abstol, 2));
  printf("---------------------------------------------\n");
  printf("inital t    = %.3g\n", T0);
  printf("first t     = %.3g\n", T1);
  printf("t step      = %.3g\n", TSTEP);
  printf("no. outputs = %d\n", NOUT);
  printf("---------------------------------------------\n\n");

  return(0);
}





static void PrintOutput(realtype t, N_Vector y)
{
realtype y1, y2, y3, y4;
y1 = Ith(y,1);
y2 = Ith(y,2);
y3 = Ith(y,3);
y4 = Ith(y,4);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e    y =[%14.6e  %14.6e  %14.6e %14.6e, ...]\n",t, y1,y2,y3,y4);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}



static void PrintRootInfo(int root_f1, int root_f2)
{
  printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);

  return;
}


static int HeaderWriteOutput(FILE* YFID, int index0, int indexend, bool fPrintTimeHeader)
{
  string whtspc = ",    ";
  string HEADERSTR;
  
  /* check file pointers */
  if (YFID == NULL) return(1);

  /* write header to file */
  HEADERSTR = "/* columns are:  ";
  if(fPrintTimeHeader){
  HEADERSTR += "t"+whtspc;
  }
  for(int i=index0; i<=indexend; i++){
    HEADERSTR += "y"+to_string(i)+whtspc;
  }
  HEADERSTR += " */\n";

  /* output header string to disk */
  fprintf(YFID, "%s", HEADERSTR.c_str());


  return(0);
}




/* Output the solution to disk (or terminal) */
static int WriteOutput(realtype t, N_Vector y, int nsupers, 
                      FILE* YFID, FILE* SDFID, FILE* EFID)
{
  realtype *ydata = N_VGetArrayPointer(y);
  //realtype *edata = N_VGetArrayPointer(e);

  /* check file pointers */
  if (YFID == NULL) return(1);

  /* output solution to disk */
  fprintf(YFID, "%24.14e,%24.14e,%24.14e,%24.14e,%24.14e\n",
          t, ydata[0], ydata[1], ydata[2], ydata[3]);
  
  /* output superdroplet solution to disk if SDFID!=NULL */
  if(SDFID)
  {
    for(int i=0; i<nsupers-1; i++)
    {
      fprintf(SDFID, "%24.14e,", ydata[i+4]);
    }
      fprintf(SDFID, "%24.14e\n", ydata[nsupers+3]);
  }
  
  // /* output error to disk */
  // if(EFID)
  // {
  // fprintf(EFID, "%24.14e,%24.14e,%24.14e,%24.14e,%24.14e\n",
  //         t, edata[0], edata[1], edata[2], edata[3]);
  // }


  return(0);
}








/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}









#endif // CVODEHELPERS