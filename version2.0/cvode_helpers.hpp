// Author: Clara Bayley
// File: cvode_helpers.hpp
/* Header file for helper function for
CVODE (sundials) ODE solver module*/


#ifndef CVODE_HELPERS_HPP
#define CVODE_HELPERS_HPP


#include <iostream>
#include <fstream>
#include <sstream>


#include "claras_SDinit.hpp"




using namespace std;


void print_output(const double t, const double y[4]);
/* print t, y to terminal */

void write_outputheader(const string solution_csv);
/* Create new .csv file with header for writing data to */


void write_output(ofstream &wfile, const double t, const double y[4]);
/* Write output t, y to wfile on disk */



#endif //CVODE_HELPERS_HPP 