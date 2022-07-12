// Author: Clara Bayley
// File: superdrops.hpp
/* Header file for superdrops class */

#ifndef SUPERDROP_HPP
#define SUPERDROP_HPP

#include <math.h>

#include "../../claras_SDconstants.hpp"
#include "common2allsuperdrops.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

class Superdrop : public Common2AllSuperdrops
/* Superdroplet Class. Is child class so each
  Superdrop instance inherits properties
  from Common2AllSuperdrops */
{
private:
  double eps0;   // initial droplet multiplicity
  double r0;     // initial droplet radius
  double m_sol0; // initial droplet solute mass

public:
  double eps;   // multiplicity of droplet
  double r;     // radius of droplet
  double m_sol; // mass of solute dissovled
  double b;     // kohler b factor

  Superdrop(){};
  /* completely empty constructor function */

  Superdrop(double aRho_l, double aRho_sol,
            double aMr_sol, double aIonic) : Common2AllSuperdrops(aRho_l, aRho_sol,
                                                                  aMr_sol, aIonic) {}
  /* constructor function with parent initialised */

  Superdrop(double aEps, double aR, double aM_sol, double aRho_l, double aRho_sol,
            double aMr_sol, double aIonic);
  /* constructor function with parent and instance initialised */

  void setSuperdropInitials(double aEps, double aR, double aM_sol);
  /* Combined Setter functions to set private attributes */

  double getR0() { return r0; }

  // double getEps0() { return eps0;}

  // double getM_sol0() { return m_sol0;}

  // double rhoeff();
  //   /* calculates effective density of droplet
  //   so mass_droplet = m = 4/3*pi*r^3 * rhoeff */

  double dry_r();
  /* calculate radius of dry droplet, ie.
  radius if it was entirely made of solute */

  double mass();
  /* calculate total mass of droplet
        mass = (water + dry areosol)  */

  double vol();
  /* volume of droplet */

  // double m_w();
  //   /* mass of only water in droplet */

  double akohler_factor(double temp);
  /* calculate a in raoult factor (exp^(a/r)) to
  account for effect of dissolved solute
  on radial growth of droplet. Using eq.6.24
  and eq.6.22 of lohmann, luond
  and mahrt intro 2 clouds textbook */

  double bkohler_factor();
  /* calculate b in kelvin factor (1-b/r^3)
  to account for curvature on radial growth
  of droplet. Using eq.6.24 and eq.6.22 of lohmann,
  luond and mahrt intro 2 clouds textbook */

  void diffusion_factors(double &fkl, double &fdl,
                         const double &p, const double &temp, const double &psat);
  /* Calculate dimensionless Fkl and Fdl
  * heat and vapour diffusion factors in
  equation for radial growth of droplets
  according to eq.7.27 lohmann, luond and
  mahrt intro 2 clouds textbook
  */

  double condensation_growth(const double &p, const double &temp,
                             const double &psat, const double &s_ratio, const double &delt);
  /* radial growth/shrink of each superdroplet due to
  condensation and diffusion of water vapour.
  [eq.7.27 lohmann intro 2 clouds] */
};

#endif // SUPERDROP_HPP