// Author: Clara Bayley
// File: superdropmodel_functions.cpp
/* Functions involved in SD model */



#include "readwritefuncs4SDM.hpp"



void print_output(const double t, const double p, 
      const double temp, const double qv, const double qc)
/* print t and kinematic data (p, temp, qv, qc) to terminal */
{
  int prec = 4; 
  cout << "t=" << fixed << setprecision(prec) << t << ", ";     
  cout << "y=[" << scientific << setprecision(prec) << p << ", ";     
  cout << scientific << setprecision(prec) << temp  << ", ";     
  cout << scientific << setprecision(prec) << qv  << ", ";     
  cout << scientific << setprecision(prec) << qc  << "]\n";     
}




void initialise_Superdrop_instances(const string FNAME, 
        Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers)
{
  /* read CSV file for eps, r and m_sol then create 
         instances of Superdroplet class */ 
 
  double eps, r, m_sol;
  double arr[3] = {0,0,0};  
  eps=arr[0]; r=arr[1]; m_sol=arr[2];
  
  ifstream file;
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
        superdrops_arr[n].eps = arr[0];
        superdrops_arr[n].r = arr[1];
        superdrops_arr[n].m_sol = arr[2];
        superdrops_arr[n].setSuperdropInitials(arr[0], arr[1], arr[2]);
        cout <<"      eps = "<< superdrops_arr[n].eps;
        cout << ", r = " << superdrops_arr[n].r;
        cout << ", m_sol = " << superdrops_arr[n].m_sol << endl;   
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


}



void write_outputheader(const string solution_csv)
/* Create new .csv file with header for writing data to */
{
  ofstream wfile;
  wfile.open(solution_csv, ios::trunc);

  wfile <<  "/* columns are dimensionless:  t,    "
            "pressure (y[0]),    temperature (y[1]),"
            "    qv (y[2]),   qc (y[3]),     */\n";

  wfile.close();

}




void write_output(ofstream &wfile, const double t, const double p, 
      const double temp, const double qv, const double qc)
/* Output t and kinematic data (p, temp, qv, qc) to wfile on disk */
{
  int prec = 16;     // set precision

  wfile << scientific << setprecision(prec) << t << ",";      
  wfile << scientific << setprecision(prec) << p << ",";      
  wfile << scientific << setprecision(prec) << temp << ",";      
  wfile << scientific << setprecision(prec) << qv << ",";      
  wfile << scientific << setprecision(prec) << qc << "\n";         
}





void write_superdrop_outputheader(const string solutionSD_csv)
/* Create new .csv file with header for writing superdroplet data to */
{
  ofstream wfile;
  wfile.open(solutionSD_csv, ios::trunc);

  wfile << "/* columns are dimensionless Superdroplet:" 
          "SD eps_1, eps_2, eps_3, ... eps_nsupers," 
          "r_1, r_2, r_3, ... r_nsupers," 
          "m_sol_1, m_sol_2, m_sol_3, ... m_sol_nsupers   */\n";

  wfile.close();

}




void write_superdrop_output(ofstream &wfile, 
            Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers)
/* Create new .csv file with header for writing superdroplet data to */
{
  
  int p = 16;     // set precision

  for(int i=0; i<nsupers; i++)  
  {
    //write SD eps output in columns[0:nsuper]
    wfile << scientific << setprecision(p) << superdrops_arr[i].eps << ",";      
  }
  
  for(int i=0; i<nsupers; i++)  
  {
    // write SD r output in columns[nsuper:2*nsuper]
    wfile << scientific << setprecision(p) << superdrops_arr[i].r << ",";
  }
  
  for(int i=0; i<nsupers-1; i++)  
  {
    // SD m_sol output in columns[2*nsuper:]
    wfile << scientific << setprecision(p) << superdrops_arr[i].m_sol << ",";
  }
  // writes last SD m_sol output in column[2*nsuper:] with '\n' instead of ','
  wfile << scientific << setprecision(p) << superdrops_arr[nsupers-1].m_sol << "\n";  
  
}



