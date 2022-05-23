#ifndef SUPERDROPLETCLASSES
#define SUPERDROPLETCLASSES

// #include <iostream>
#include <math.h>


// using namespace std;
// using namespace DropletConstants;



/* constants that belong in parts of equations
Eg. kohler_factor a = KOH_A/temp */
const double KOH_A    = 3.3e-7;             // used in kohler factor a calculation
const double KOH_B    = 43e-6;        // used in kohler factor b calculation





class Common2AllSuperdrops
/* Parent class for all superdroplets. Each 
      Superdrop instance inherits these properties'' */

{
    private:
        double rho_l;         //density of liquid in droplets [Kg/m^3]
        double rho_sol;       // density of solute in droplets [Kg/m^3]
        double mr_sol;        // Mr of solute [g/mol]
        double ionic;         // degree ionic dissociation (van't Hoff factor)

    public:
        

    Common2AllSuperdrops(double rho_l, double rho_sol, double mr_sol, double ionic)
    // constructor function called everytime instance of class is created
    {
        setPrivates(rho_l, rho_sol, mr_sol, ionic);

    };

    void setPrivates(double aRho_l, double aRho_sol, double aMr_sol, double aIonic)
    // Combined Setter functions to set private attributes
    {
        rho_l = aRho_l;
        rho_sol = aRho_sol;
        mr_sol = aMr_sol;
        ionic = aIonic;
    }

    double getRho_l(){
        return rho_l;
    }

    double getRho_sol(){
        return rho_sol;
    }
    
    double getMr_sol(){
        return mr_sol;
    }
    
    double getIonic(){
        return ionic;
    }

};




class Superdrop : public Common2AllSuperdrops
/* Superdroplet Class. Is child class so each 
      Superdrop instance inherits properties
      from Common2AllSuperdrops */
{
    private:
        double dry_r;          // dry radius of droplet
        double addsol;         // contribution to density due to solute = addsol/r^3
        double r0;             // initial droplet radius
        double m0;             // initial droplet mass
        double m_w0;           // initial dropelt water content


    public:
        int eps;               // multiplicity of droplet
        double r;              // radius of droplet [m]
        double m_sol;          // mass of solute dissovled [Kg]
        double b;              // kohler b factor = 43e-6 m^3/mol * mass_solute/ Mr_solute 
        

    Superdrop(int aEps, double aR, double aM_sol, double aRho_l, double aRho_sol, 
                double aMr_sol, double aIonic) : Common2AllSuperdrops(aRho_l, aRho_sol, 
                            aMr_sol, aIonic)
    {



       eps = aEps;       
       r = aR;
       m_sol = aM_sol; 

       b = KOH_B*m_sol*getIonic()/getMr_sol(); 
       
       setSuperdropPrivates(aR);

    }



   void setSuperdropPrivates(double aR)
    // Combined Setter functions to set private attributes
    {
       dry_r = pow(3*m_sol/(4*M_PI*getRho_sol()), 1.0/3);
       addsol = (getRho_sol() - getRho_l())*pow(dry_r, 3.0);
       r0 = aR;
       m_w0 = 4.0/3 * M_PI * getRho_l() * pow((r0 - dry_r), 3.0);
       m0 = m_sol + m_w0;
    }
    
    double getAddsol(){
        return addsol;
    }

    double getDry_r(){  
        return dry_r;
    }
    
    double getR0(){
        return r0;
    }
   
    double getM_w0(){
        return m_w0;
    }

    double getM0(){
        return m0;
    }



    double rhoeff(){
    /* calculates effective density of droplet
    so mass_droplet = m = 4/3*pi*r^3 * rhoeff */
    
        double rhoeff;
        rhoeff = getRho_l() + addsol/pow(r,3.0);

        return rhoeff;
    }
    

    double m(){
    // total mass of droplet (water + dry areosol) 
        return 4.0/3 * M_PI * rhoeff() * pow(r,3.0);
    }


    double m_w(){
    // mass of only water in droplet
        double v_w;
        v_w = 4.0/3 * M_PI *  (pow(r,3.0) - pow(dry_r, 3.0));
        return getRho_l() * v_w;
    }
    
    double vol(){
    // volume of droplet
        return 4.0/3 * M_PI * pow(r,3.0);
    } 
    

    void kohler_factor(double temp, double* ad_a, double* ad_b)
    /* calculate b in kelvin factor (1-b/r^3)
        and a in raoult factor (exp^(a/r)) to
        account for curvature and effect of solute
        on radial growth of droplet respectively.
        Using eq.6.24 and eq.6.22 of lohmann, luond
        and mahrt intro 2 clouds textbook */
    {
        *ad_a = KOH_A/temp;                      // 2*surface_tens (=0.0756N/m)/(rho_l*rgas_v*temp) [eq.6.24]
        *ad_b = b;
    }



};












#endif //SUPERDROPLETCLASSES