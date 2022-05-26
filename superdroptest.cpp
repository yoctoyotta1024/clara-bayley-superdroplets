#include <iostream>
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */

#include "superdropletclasses.hpp"
#include "constants.hpp"


namespace cnst = dimmed_constants;
using namespace std;

int main(){


    cout << "making drops" << endl;


    double eps = 2;
    double r = 6;
    double m_sol = 9*8800*M_PI;
    double RHO_SOL = 2200;
    double MR_SOL = 0.058443;
    double IONIC = 2;

    Superdrop drop1(eps, r, m_sol, cnst::RHO_L, RHO_SOL, MR_SOL, IONIC);

    cout << drop1.m_sol << endl;
    cout << drop1.getM_w0() << endl;
    cout << drop1.getM0() << endl;
    cout << drop1.getDry_r() << endl;
    

    cout << "-- initialise a droplet ---" << endl;
    cout << drop1.getR0() << endl;
    cout << drop1.r << endl;
    cout << drop1.rhoeff() << endl;
    cout << drop1.m() << endl;
    drop1.r = 5;
    cout << "-- change r of a droplet ---" << endl;
    cout << drop1.getR0() << endl;
    cout << drop1.r << endl;
    cout << drop1.rhoeff() << endl;
    cout << drop1.m() << endl;
    cout << "-----" << endl;


    double a = 3;
    double b = 4;
    double temp = 300;
    cout << a << " a&b " << b << endl;
    drop1.kohler_factor(temp, &a, &b);
    cout << a << " a&b " << b << endl;
    cout << "--- calc kohler factors using pointers to assign values --" << endl;


    double eps2 = 2;
    double eps3 = 2;
    double r2 = 0.654;
    double r3 = 0.33;
    Superdrop drop2(eps2, r2, m_sol, cnst::RHO_L, RHO_SOL, MR_SOL, IONIC);
    Superdrop drop3(eps3, r3, m_sol, cnst::RHO_L, RHO_SOL, MR_SOL, IONIC);
    cout << drop2.r << " ,"<< drop3.r << endl;
    cout << drop2.eps << " ," << drop3.eps << endl;
    cout << "-- initialse many drops ---" << endl;


    Superdrop obj[2] = {Superdrop(eps2, r2, m_sol, cnst::RHO_L, RHO_SOL, MR_SOL, IONIC),
                    Superdrop(eps3, r3, m_sol, cnst::RHO_L, RHO_SOL, MR_SOL, IONIC) };

    Superdrop drops[2] = {drop2, drop3};
    
    for(int i=0; i<2; i++)
    {
        cout << drops[i].r << endl;
        cout << obj[i].r << endl;
        cout << drops[i].eps << endl;
        cout << obj[i].eps<< endl;
    }
    cout << "-- acces attributes of many drops once they are stored in array ---" << endl;


    FILE *YFID = NULL;        /* solution output file */
    savey = SAVENAME+"_sol.csv";
    YFID = fopen(savey.c_str(),"w");
    EFID = fopen(saverr.c_str(),"w");

    fprintf(YFID, "%24.14e,%24.14e,%24.14e\n", eps2, r2, m_sol);

    fclose(YFID)

    return 0;
}