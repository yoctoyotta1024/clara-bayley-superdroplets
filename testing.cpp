// run with g++ -I /usr/local/sundials-6/include  testing.cpp -std=c++11 
#include <iostream>
#include <cmath>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include "constants.hpp"
#include "init.hpp"


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




/* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants */
typedef struct {
  realtype w, drop1, drop2, drop3;
} *UserData;


static void InitUserData(UserData data);
realtype f(realtype z, void *user_data);



static void InitUserData(UserData data, realtype w)
{
  data->w = w;
  data->drop1 = 1;
  data->drop2 = 2;
  data->drop3 = 3;
}


int main(){

  realtype w = iW/dlc::W0;
  cout << w << endl;
  UserData data;
  data = (UserData) malloc(sizeof *data);
  InitUserData(data, w);

  realtype y;
  y = f(3, data);
    
  cout << "y: " << y<< endl;

  free(data);
}




realtype f(realtype z, void *user_data){

  UserData data;
  realtype w;
  data = (UserData) user_data;
  
  w = data -> w;

  realtype y = w+2;

  return y;
}