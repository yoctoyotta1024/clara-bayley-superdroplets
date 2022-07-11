// Author: Clara Bayley
// File: collisions4SDM.hpp
/* Header file for functions involved
in superdroplet model */



#ifndef COLLISIONS4SDM_HPP
#define COLLISIONS4SDM_HPP


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include "../../claras_SDinit.hpp"
#include "common2allsuperdrops.hpp"
#include "superdrop.hpp"

using namespace std;



void collide_droplets(const int nsupers, const int nhalf, vector<int> pvec,
       Superdrop (&superdrops_arr)[init::NSUPERS]);





#endif //COLLISIONS4SDM_HPP