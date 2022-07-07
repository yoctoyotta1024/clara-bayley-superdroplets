// run with g++ -I /usr/local/sundials-6/include  testing.cpp -std=c++11 
#include <iostream>
#include <vector>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include "constants.hpp"
#include "init.hpp"
#include "superdroplets.hpp"
#include "collisions.hpp"

namespace dlc = dimless_constants;




#define T0    RCONST(tspan[0]/dlc::TIME0)        // initial time (dimensionless)          
#define TSTEP RCONST(tspan[1]/nout/dlc::TIME0)   // output time step (dimensionless)     
#define NOUT  nout                               // number of output times
#define COLL_TSTEP RCONST(coll_tstep/dlc::TIME0)                    // No. droplet collisions events

#define SDloop(i,nsupers) for(int i=0; i<nsupers; i++)  //for loop over all superdroplets




#define T0    RCONST(tspan[0]/dlc::TIME0)        // initial time (dimensionless)          
#define TSTEP RCONST(tspan[1]/nout/dlc::TIME0)   // output time step (dimensionless)     
#define NOUT  nout                               // number of output times

#define COLL_TSTEP RCONST(coll_tstep/dlc::TIME0)                    // No. droplet collisions events


int main(){

  realtype t, tout, CollsPerTstep;
  
  double tspan[2]    = {0, 4000};                   // time span of integration [s]
  int nout           = 20;                         // No. time points to evaluate (save data at)
  double coll_tstep  = 100;                         // maximum time between each droplet collisions event [s]

  CollsPerTstep = TSTEP/(COLL_TSTEP);
  tout = T0+COLL_TSTEP;   // first output time = t0 + TSTEP/CollsPerTstep

  cout << T0 << endl;
  double delt = 0;
  

  while(delt < COLL_TSTEP){
      
    cout << "do integration" << endl;

    delt += COLL_TSTEP;
      
    cout << delt << endl;
  }

    tout += delt;
  
  }


  return 0;
}





