#ifndef SUPERDROPLETCLASSES
#define SUPERDROPLETCLASSES

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;



/* constants that belong in parts of equations
Eg. kohler_factor a = KOH_A/temp */
//#include "constants.hpp"
//namespace dlc = dimless_constants;
const double KOH_A    = 3.3e-7/(dlc::TEMP0*dlc::R0);             // used in dimensionless kohler factor a calculation
const double KOH_B    = 43e-6*dlc::RHO0/dlc::MR0;                // used in dimensionless kohler factor b calculation




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
        
    Common2AllSuperdrops() {}
    
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
        double eps0;           // initial droplet mulitplicity
        double r0;             // initial droplet radius
        double m_sol0;         // initial droplet solute mass

    public:
        double eps;            // multiplicity of droplet
        double r;              // radius of droplet
        double m_sol;          // mass of solute dissovled
        double b;              // kohler b factor
        

    Superdrop() {}

    Superdrop(double aRho_l, double aRho_sol, 
                double aMr_sol, double aIonic) : Common2AllSuperdrops(aRho_l, aRho_sol, 
                            aMr_sol, aIonic){}

 
    Superdrop(double aEps, double aR, double aM_sol, double aRho_l, double aRho_sol, 
                double aMr_sol, double aIonic) : Common2AllSuperdrops(aRho_l, aRho_sol, 
                            aMr_sol, aIonic)
    {

       eps = aEps;       
       r = aR;
       m_sol = aM_sol; 

       setSuperdropPrivates(aEps, aR, aM_sol);

    }



   void setSuperdropPrivates(double aEps, double aR, double aM_sol)
    // Combined Setter functions to set private attributes
    {
       dry_r = pow(3*m_sol/(4*M_PI*getRho_sol()), 1.0/3);
       eps0 = aEps;
       r0 = aR;
       m_sol0 = aM_sol;
    }
    

    double getDry_r(){  
        return dry_r;
    }
    
    double getR0(){
        return r0;
    }

    double getEps0(){
        return eps0;
    }

    double getM_sol0(){
        return m_sol0;
    }
   

    double rhoeff(){
    /* calculates effective density of droplet
    so mass_droplet = m = 4/3*pi*r^3 * rhoeff
    taking into account mass of liquid and mass of
    solute assuming solute occupies volume it
    would given (dry) density */
    
        double rhoeff;
        rhoeff = 3*m_sol/(4.0*M_PI*pow(r, 3.0))*(1-getRho_l()/getRho_sol());
        rhoeff += getRho_l();

        return rhoeff;
    }
    
    double vol(){
    // volume of droplet
        return 4.0/3 * M_PI * pow(r,3.0);
    } 

    double m(){
    /* total mass of droplet (water + dry areosol) =
        m =  4/3*pi*rho_l**3 + m_sol(1-rho_l/rho_sol) 
        ie. m = 4/3*pi*rhoeff*R**3 */
        
        double m;
        m = m_sol*(1-getRho_l()/getRho_sol());
        m += 4/3.0*M_PI*pow(r, 3.0)*getRho_l();                                  

        return m;
    }


    double m_w(){
    /* mass of only water in droplet */
        double v_w;
        v_w = -m_sol/getRho_sol();
        v_w += 4/3.0*M_PI*pow(r, 3.0);

        return getRho_l()*v_w;
    }
    

    

    double akohler_factor(double temp)
    /* calculate a in raoult factor (exp^(a/r)) to
        account for effect of dissolved solute
        on radial growth of droplet. Using eq.6.24 
        and eq.6.22 of lohmann, luond
        and mahrt intro 2 clouds textbook */
    {
        return KOH_A/temp;                      // 2*surface_tens (=0.0756N/m)/(rho_l*rgas_v*temp) [eq.6.24]
    }

    double bkohler_factor()
    /* calculate b in kelvin factor (1-b/r^3)
        to account for curvature on radial growth 
        of droplet. Using eq.6.24 and eq.6.22 of lohmann, 
        luond and mahrt intro 2 clouds textbook */
    {
        return KOH_B*m_sol*getIonic()/getMr_sol();
    }



};










Superdrop* initialise_Superdrop_instances(string FNAME, Superdrop* ptr, int nsupers)
{

  /* read CSV file for eps, r and m_sol then create 
         instances of Superdroplet class */ 
  ifstream file;
  double eps, r, m_sol;
  double arr[3] = {0,0,0};  
  eps=arr[0]; r=arr[1]; m_sol=arr[2];
  string line, substr;
  vector<string> result;
  bool headerend= false;
  int hend, k;
  int l = 0; 
  int n = 0;

  file.open(FNAME);
  cout << "Reading in Superdroplet Initialisation Data... \n" << endl;

  while(!file.eof()){   // while reading the file returns non-empty lines
    
    file>>line;
      
    if(line == "*/"){
      // find end of header by first time that line = */
      headerend = true;
      hend = l;
    }

    if(line == ""){
      // stop reading data if a line is empty
      break;
    }
        
    if(headerend && l > hend){
      // all lines after end of header line are data
      cout<< "SD"<< l-hend << ": "<<line<< endl;    
    
      stringstream ss(line);
      
      k = 0;    
      while(ss.good())
      {
        getline(ss, substr, ',' );
        result.push_back(substr);
        arr[k] = atof(substr.c_str());
        k++;
      }
      if(n < nsupers){
        (ptr+n) -> eps = arr[0];
        (ptr+n) -> r = arr[1];
        (ptr+n) -> m_sol = arr[2];
        (ptr+n) -> setSuperdropPrivates(arr[0], arr[1], arr[2]);
        cout <<"      eps = "<< (ptr+n) -> eps;
        cout << ", r = " << (ptr+n) -> r;
        cout << ", m_sol = " << (ptr+n) -> m_sol << endl;   
      n++;
      }
      else{
        cout << "Superdrop not made" << endl;
      }
    }
    
    line = "";
    l++;
  }

  cout << "\n ------------------------------- " << endl;
  cout << " -- Datafile reading complete -- " << endl;
  cout << "  No. Superdroplets Created = " << nsupers << " out of "<< l-hend-1 <<endl;
  if(n < nsupers){
    cout << "\n!!!! !!!! !!!! !!!! ERROR WARNING !!!! !!!! !!!! !!!!"<< endl;
    cout << "  !! Some bad superdrops have been created !!" << endl;
    cout << "     No. superdrops in array (nsupers) greater" << endl;
    cout << "     than No. lines of data in .csv file." << endl;
    cout << "!!!! !!!! !!!! !!!! ERROR WARNING !!!! !!!! !!!! !!!!\n"<< endl;
  }
  else if (l-hend-1 > nsupers){
    cout << "\n!!!! !!!! !!!! !!!! WARNING !!!! !!!! !!!! !!!!"<< endl;
    cout << "Fewer superdrops have been created than is possible" << endl;
    cout << "    No. superdrops in array (nsupers) less than" << endl;
    cout << "      No. lines of data in .csv file." << endl;
    cout << "!!!! !!!! !!!! !!!! WARNING !!!! !!!! !!!! !!!!\n"<< endl;
  }
  cout << " ------------------------------- \n" << endl;
  
  file.close();

  return ptr;
}







#endif //SUPERDROPLETCLASSES