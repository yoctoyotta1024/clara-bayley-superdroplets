// Author: Clara Bayley
// File: cvode_helpers.cpp
/* Helper function for using CVODE
(sundials) ODE solver module */



#include "cvode_helpers.hpp"



void print_output(const double t, const double y[4])
/* print t, y to terminal */
{

  int p = 4; 
  cout << "t=" << fixed << setprecision(p) << t << ", ";     
  cout << "y=[" << scientific << setprecision(p) << y[0] << ", ";     
  cout << scientific << setprecision(p) << y[1]  << ", ";     
  cout << scientific << setprecision(p) << y[2]  << ", ";     
  cout << scientific << setprecision(p) << y[3]  << "]\n";     
}



void write_outputheader(const string solution_csv)
/* Create new .csv file with header for writing data to */
{
  ofstream wfile;
  wfile.open(solution_csv, ios::trunc);

  wfile <<  "/* columns are dimensionless:  t,    P (y[0]),    T (y[1]),"
              "    qv (y[2]),   qc (y[3]),     */\n";

  wfile.close();

}



void write_output(ofstream &wfile, const double t, const double y[4])
/* Write output t, y to wfile on disk */
{
  int p = 16;     // set precision

  wfile << scientific << setprecision(p) << t << ",";      
  wfile << scientific << setprecision(p) << y[0] << ",";      
  wfile << scientific << setprecision(p) << y[1] << ",";      
  wfile << scientific << setprecision(p) << y[2] << ",";      
  wfile << scientific << setprecision(p) << y[3] << "\n";         
}

