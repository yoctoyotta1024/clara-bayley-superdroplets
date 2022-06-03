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




int main(){

  int nhalf = floor(nsupers/2);
  int scale_p = nsupers*(nsupers-1)/(2*nhalf);

  int retval2;
  string INITDROPSCSV;
  INITDROPSCSV = "dimlessSDinit.csv";                                       // file to read for SD eps, r and m_sol initialisation 

  /* Initialise Superdroplets using INITDROPSCSV .csv file */
  Superdrop drops_arr[nsupers];
  SDloop(i, nsupers) 
  {
    drops_arr[i] = Superdrop(iRho_l, iRho_sol, iMr_sol, iIONIC); 
  }
  Superdrop* ptr; 
  ptr = &drops_arr[0];
  initialise_Superdrop_instances(INITDROPSCSV, ptr, nsupers);
  

  vector<int> numvec(nsupers);
  SDloop(i, nsupers){
    numvec[i] = i;
    cout << i << ": " << ptr + i << endl;
  }
  cout << "initial adresses" << endl;

  SDloop(i, nsupers){
    numvec[i] = i;
    cout << ptr + i  <<" initial eps: " << (ptr + i) -> eps;
    cout << ", initial r: " << (ptr + i) -> r;
    cout << ", initial m_sol: " << (ptr + i) -> m_sol << endl;
  } 
  cout << "---" << endl;
  retval2 = collide_droplets(nsupers, nhalf, scale_p, ptr, numvec);
  cout << "---" << endl;

  SDloop(i, nsupers){
    numvec[i] = i;
    cout << ptr + i  << " eps: " << (ptr + i) -> eps;
    cout << ", r: " << (ptr + i) -> r;
    cout << ", m_sol: " << (ptr + i) -> m_sol << endl;
  } 

  return 0;
}








