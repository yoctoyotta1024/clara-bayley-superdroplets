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



int diffusion_factors(realtype* fkl, realtype* fdl, realtype temp, realtype p, realtype psat)
{
  realtype Thermk, Diffuse_v;
  realtype Temp = temp*dlc::TEMP0;
  realtype P = p*dlc::P0;
  realtype Psat = psat*dlc::P0;

  /* dimensionless factors for condensation-diffusional growth equation */
  Thermk = 7.11756e-5*pow(Temp, 2.0) + Temp*4.38127686e-3;      // [eq.7.24 lohmann intro 2 clouds]
  Diffuse_v = (4.012182971e-5 / P * pow(Temp,1.94))/DC::RGAS_V;

  *fkl = (DC::LATENT_RGAS_V/Temp-1)*DC::LATENT_V/(Thermk*dlc::F0); 
  *fdl = Temp/(Diffuse_v*Psat)/dlc::F0;
    
  
  return 0;
}



int main(){

  realtype fkl, fdl;
  diffusion_factors(&fkl, &fdl, 1,1,611.2126978267946/P0);
  cout << fkl << "pointed"<<endl;
  cout << fdl <<endl;
}




