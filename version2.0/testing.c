// Author: Clara Bayley
// File: main.cpp
/* This file runs the entire superdrop (SD) model */


#include <iostream>
#include <fstream>

#include "claras_SDinit.hpp"
#include "superdrop.hpp"
#include "common2allsuperdrops.hpp"
#include "readwritefuncs4SDM.hpp"

#include "claras_SDconstants.hpp"


//using namespace std;
namespace dlc = dimless_constants;



void simple_reference(Superdrop (&arr)[init::NSUPERS], int nsupers)
{

  for(int i=0; i<floor(nsupers/2.0); i++)  //for loop over all superdroplets
  { 
    Superdrop& drop1 = arr[2*i];
    Superdrop& drop2 = arr[2*i+1];
    drop1.eps = 2*i;
    drop2.eps = 2*i+1;
  }
}






int main()
{

  Superdrop superdrops_arr[init::NSUPERS];
  const int nsupers = init::NSUPERS;                  // allows flexible no. of superdroplets to be included in model (nsupers <= NSUPERS)

  /* Initialise Superdroplets using initdrops_csv file */
  for(int i=0; i<nsupers; i++)                  //for loop over all superdroplets
  {
    superdrops_arr[i] = Superdrop(init::iRho_l, 
                init::iRho_sol, init::iMr_sol, init::iIONIC); 
  }
  initialise_Superdrop_instances(init::initdrops_csv,
                                  superdrops_arr, nsupers);


  // How to create a reference is this: Superdrop (&ref)[init::NSUPERS] = superdrops_arr;

  for(int i=0; i<nsupers; i++)  //for loop over all superdroplets
  { 
    cout << superdrops_arr[i].eps<< endl;   
  }

  simple_reference(superdrops_arr, nsupers);

  for(int i=0; i<nsupers; i++)  //for loop over all superdroplets
  { 
    cout << superdrops_arr[i].eps<< endl;   
  }








  return 0;
}

