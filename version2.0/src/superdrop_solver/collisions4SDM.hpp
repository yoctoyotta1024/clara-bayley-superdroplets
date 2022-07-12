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
/* superdroplet collision scheme according to shima et al. 2009.
  Monte-carlo collisions of random pairs of SDs from list of 
  all superdroplets in given volume */

static Superdrop& get_dropletref(const string whichdrop,
                Superdrop &dropA, Superdrop &dropB);
/* compare dropA.eps with dropB.eps and return 
  either drop1 or drop2 such that drop1.eps is always >
  drop2.eps */

#endif // COLLISIONS4SDM_HPP