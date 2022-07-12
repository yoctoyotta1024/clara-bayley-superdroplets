// Author: Clara Bayley
// File: common2allsuperdrops.hpp
/* Header file for overarching class
"common to all superdrops" that is parent
of superdrop class */

#ifndef COMMON2ALLSUPERDROPS_HPP
#define COMMON2ALLSUPERDROPS_HPP

class Common2AllSuperdrops
/* Parent class for all superdroplets. Each
      Superdrop instance inherits these properties'' */

{
private:
    double rho_l;   // density of liquid in droplets [Kg/m^3]
    double rho_sol; // density of solute in droplets [Kg/m^3]
    double mr_sol;  // Mr of solute [g/mol]
    double ionic;   // degree ionic dissociation (van't Hoff factor)

public:
    Common2AllSuperdrops(){};
    /* constructor function called everytime
        empty instance of class is created */

    Common2AllSuperdrops(double rho_l, double rho_sol,
                         double mr_sol, double ionic);
    /* constructor function called everytime
        instance of class is created */

    void setPrivates(double aRho_l, double aRho_sol,
                     double aMr_sol, double aIonic);
    /* Combined Setter functions to set
    private attributes */

    double getRho_l() { return rho_l; }

    double getRho_sol() { return rho_sol; }

    double getMr_sol() { return mr_sol; }

    double getIonic() { return ionic; }
};

#endif // COMMON2ALLSUPERDROPS_HPP