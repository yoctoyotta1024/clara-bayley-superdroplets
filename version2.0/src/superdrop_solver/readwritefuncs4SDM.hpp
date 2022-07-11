// Author: Clara Bayley
// File: readwritefuncs4SDMM.hpp
/* Header file for functions involved
in superdroplet model */



#ifndef READWRITEFUNCS4SDM_HPP
#define READWRITEFUNCS4SDM_HPP



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


#include "../../claras_SDinit.hpp"
#include "common2allsuperdrops.hpp"
#include "superdrop.hpp"

using namespace std;




void print_output(const double t, const double y[4]);
/* print t, y to terminal */


void initialise_Superdrop_instances(string FNAME, 
        Superdrop (&superdrops_arr)[init::NSUPERS], int nsupers);
/* read .csv file for initial eps, r and m_sol and then create 
	instances of Superdroplet class */ 


void write_outputheader(const string solution_csv);
/* Create new .csv file with header for writing data to */


void write_output(ofstream &wfile, const double t, const double y[4]);
/* Write output t, y to wfile on disk */


void write_superdrop_outputheader(const string solutionSD_csv);
/* Create new .csv file with header for writing superdroplet data to */


void write_superdrop_output(ofstream &wfile, 
          Superdrop (&superdrops_arr)[init::NSUPERS], const int nsupers);
/* Create new .csv file with header for writing superdroplet data to */








#endif //READWRITEFUNCS4SDM_HPP