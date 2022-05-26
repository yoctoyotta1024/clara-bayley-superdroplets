#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "superdropletclasses.hpp"
#include "constants.hpp"

namespace dlc = dimless_constants;
using namespace std;


/* droplet init params. First create superdroplet eps, r0
and m_sol data using python "create_superdrop_init.py" */
int nsupers       = 20;                // no. distinct superdrop objects in array
double Rho_l      = dlc::Rho_l;
double Rho_sol    = dlc::Rho_sol;
double Mr_sol     = dlc::Mr_sol;
int IONIC         = dlc::IONIC;                 


int create_drop_arr(int nsupers, double Rho_l, double Rho_sol,
                    double Mr_sol, int IONIC); 



int main(){

  
  create_drop_arr(nsupers, Rho_l, Rho_sol, Mr_sol, IONIC);
  
  return 0;

}






int create_drop_arr(int nsupers, double Rho_l, double Rho_sol,
                    double Mr_sol, int IONIC)
{

  /* read CSV file for eps, r and m_sol then create 
         instances of Superdroplet class */

  ifstream file;
  string FNAME;
  FNAME = "dimlessinit_superdroplets.csv";

  Superdrop drops_arr[nsupers];

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

      eps = arr[0]; r=arr[1]; m_sol = arr[2];
      if(n < nsupers){
        drops_arr[n] = Superdrop(eps, r, m_sol, Rho_l, Rho_sol, Mr_sol, IONIC); 
        cout <<"      eps = "<< drops_arr[n].eps;
        cout << ", r = " << drops_arr[n].r;
        cout << ", m_sol = " << drops_arr[n].m_sol << endl;   
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

  return 0;
}



