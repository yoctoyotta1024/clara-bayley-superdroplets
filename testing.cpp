// run with g++ -I /usr/local/sundials-6/include  testing.cpp -std=c++11 
#include <iostream>
#include <vector>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include "constants.hpp"
#include "init.hpp"
#include "superdroplets.hpp"
#include "collisions.hpp"

namespace dlc = dimless_constants;
using namespace dlc;
using namespace std;



#define T0    RCONST(tspan[0]/dlc::TIME0)        // initial time (dimensionless)          
#define TSTEP RCONST(tspan[1]/nout/dlc::TIME0)   // output time step (dimensionless)     
#define NOUT  nout                               // number of output times
#define COLL_TSTEP RCONST(coll_tstep/dlc::TIME0)                    // No. droplet collisions events

#define SDloop(i,nsupers) for(int i=0; i<nsupers; i++)  //for loop over all superdroplets


static int WriteSetup2Txt(string INITHPP, string CONSTSHPP, FILE* IFID);


int main(){

  string SAVENAME = "sundials2";                                                  // solution written to SAVENAME_sol.csv file

  FILE *IFID = NULL;           // solution output file 
  IFID = fopen((SAVENAME+"_setup.txt").c_str(),"w");

  string INITHPP, CONSTSHPP;
  INITHPP = "init.hpp";                                       // file to read for SD eps, r and m_sol initialisation 
  CONSTSHPP = "constants.hpp";                                       // file to read for SD eps, r and m_sol initialisation 
  
  WriteSetup2Txt(INITHPP, CONSTSHPP, IFID); 



  return 0;
}








  static int WriteSetup2Txt(string INITHPP, string CONSTSHPP, FILE* IFID)
{

  ifstream initfile;
  ofstream writefile;
  string header, line;
  string breakheader = "// ----------------------------- //\n";

  /* check file to write to pointers */
  if (IFID == NULL) return(1);

  /* read .hpp init file */ 
  string readfiles[] = {INITHPP, CONSTSHPP};

  for(int i=0; i<2; i++){
    
    initfile.open(readfiles[i]);
    cout << i << " writing " << readfiles[i]<< " to XXX_setup.txt"<<endl;

    fprintf(IFID, "%s", breakheader.c_str()); 
    header = "// --------- "+readfiles[i]+" --------- //\n";
    fprintf(IFID, "%s", header.c_str());
    fprintf(IFID, "%s", breakheader.c_str()); 

    while(getline(initfile, line)){   // read file line by line
    
      /* output lines to .txt file on disk */
      fprintf(IFID, "%s", (line+"\n").c_str());
    }

    fprintf(IFID, "%s", (breakheader+"\n\n\n").c_str()); 

    initfile.close();
  }





  return(0);
}