// Author: Clara Bayley
// File: condensationgrowth4SDM.cpp
/* Functionality for modelling condensation-
   diffusional growth of superdroplets */

#ifndef CONDENSATIONGROWTH4SDM_HPP
#define CONDENSATIONGROWTH4SDM_HPP

#include <iostream>
#include <math.h>

#include "../../claras_SDinit.hpp"
#include "../../claras_SDconstants.hpp"
#include "common2allsuperdrops.hpp"
#include "superdrop.hpp"

// using namespace std;
namespace dlc = dimless_constants;

double pv2qv(const double pv, const double p);
/* Calculate mass mixing ratio
qv = m_v/m_dry = rho_v/rho_dry
given vapour pressure pv = p_v/p_tot. */

double cp_moist(const double qv, const double qc);
/* (dimensionless) specific heat capacity of
moist parcel of air */

double saturation_pressure(const double temp);
/* Calculate the equilibrium vapor pressure
of water over liquid water ie. the
saturation pressure (psat). Equation taken from
python module typhon.physics.thermodynamics.e_eq_water_mk
with conversion to real temp /K = T*Temp0 and from
real psat to dimensionless psat = psat/P0. */

void condensation_onto_superdroplets(const double delt, double &p,
                                     double &temp, double &qv, double qc,
                                     double &deltemp, double &delqv, double &delqc,
                                     Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers);
/* Change to temp, qv and qc due to change in radius of
all superdroplets via diffusion and condensation
of water vapour. */

#endif // CONDENSATIONGROWTH4SDM_HPP