// run with g++ -I /usr/local/sundials-6/include  testing.cpp -std=c++11 
#include <iostream>
#include <cmath>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include "constants.hpp"


namespace dlc = dimless_constants;
using namespace dlc;
using namespace std;


/* ---- constants for use in dp_dt function ----- */
realtype dp_dt_const = dlc::W0*dlc::TIME0*DC::G/(DC::RGAS_DRY*dlc::TEMP0); 
realtype lpsrate = 0.0062/dlc::TEMP0*dlc::W0*dlc::TIME0;
realtype tempg = 273.15/dlc::TEMP0;
realtype pg = 100000/dlc::P0;
realtype zg = 0/(dlc::W0*dlc::TIME0);
realtype gamma = DC::G/(DC::RGAS_DRY*0.0062)-1;
/* ----------------------------------------------- */







int main(){
    
  realtype t, w, dp, profile;

  t = 0.001001001001001001;
  w = 0.5/dlc::W0;

  profile = 1 - lpsrate/tempg*(w*t-zg);
  profile = pow(profile, gamma);

  dp = -dp_dt_const*pg/tempg*profile;
  // dp = -dp_dt_const*pg/tempg*profile;

  cout << t << " " << w << endl;
  cout << profile << endl;
  cout << gamma << endl;
  cout << dp << endl;

}