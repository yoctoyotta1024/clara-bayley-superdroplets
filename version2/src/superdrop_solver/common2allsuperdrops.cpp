// Author: Clara Bayley
// File: common2allsuperdrops.cpp
/* Functionality file for overarching class
"common to all superdrops" that is parent
of superdrop class */

#include "common2allsuperdrops.hpp"

Common2AllSuperdrops::Common2AllSuperdrops(double rho_l,
           double rho_sol, double mr_sol, double ionic)
/* constructor function called everytime
    instance of class is created */
{
    setPrivates(rho_l, rho_sol, mr_sol, ionic);
};

void Common2AllSuperdrops::setPrivates(double aRho_l,
         double aRho_sol, double aMr_sol, double aIonic)
/* Combined Setter functions to set private attributes */
{
    rho_l = aRho_l;
    rho_sol = aRho_sol;
    mr_sol = aMr_sol;
    ionic = aIonic;
};