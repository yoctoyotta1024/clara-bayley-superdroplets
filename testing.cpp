// run with g++ -I /usr/local/sundials-6/include  testing.cpp -std=c++11 
#include <iostream>
#include <cmath>
#include <vector>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include "constants.hpp"
#include "init.hpp"

namespace dlc = dimless_constants;
using namespace dlc;
using namespace std;



#define T0    RCONST(tspan[0]/dlc::TIME0)        // initial time (dimensionless)          
#define TSTEP RCONST(tspan[1]/nout/dlc::TIME0)   // output time step (dimensionless)     
#define NOUT  nout                               // number of output times


#define COLL_TSTEP RCONST(coll_tstep/dlc::TIME0)                    // No. droplet collisions events
int main(){

  realtype tout, CollsPerTstep;
  
  cout << TSTEP << endl;
  cout << COLL_TSTEP<< endl;
  CollsPerTstep = TSTEP/(COLL_TSTEP);

  cout << CollsPerTstep << endl;


  tout = T0;   // first output time = t0 + TSTEP
  cout << "writing output data" << endl;
  cout << "----- "<< T0*dlc::TIME0<<" -------" << endl;
  for (int iout = 0; iout < NOUT; iout++){
    
    for(int j=0; j<ceil(CollsPerTstep); j++){
      cout << "coliision event " << j << endl;
      cout << "Advance solution in time" << endl;
    
    tout += TSTEP/CollsPerTstep;
    cout << "tout "<< tout << endl;
    }

    cout << "writing output data, iout "<< iout << endl;
    cout << "----- "<< tout*dlc::TIME0<<" -------" << endl;
    //tout += TSTEP;

  }


}





