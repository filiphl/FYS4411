#include <iostream>
#include "solver.h"
#include <iomanip>
#include <time.h>

using namespace std;

int main()
{
    double t0 = clock();
    int N = 500;
    int D = 1;
    double E;

    Solver s(N,D);

    double E_sum = s.localenergy();
    double E_sum2 = E_sum*E_sum;
    double anal = s.Analytical();
    cout << "Numerical: "<< E_sum <<endl;
    cout << "Analytical: "<< anal<<endl;

    double localE;
    double expE = 0.0;


    // Monte Carlo loop
    int nCycles = 1000;
    for (int i=0; i<nCycles; i++){
        if (i%100 ==0) {
            cout << "  step # "<< i << "\r";
            fflush(stdout);
        }
        s.metropolis_step();
        E = s.localenergy();
        localE = s.localenergy();
        expE += localE;
        E_sum += E;
        E_sum2 += E*E;
    }
    cout << endl;


    nCycles++;
    expE = expE/nCycles;

    E = E_sum/nCycles;
    double E2 = E_sum2/nCycles;
    cout << "energy=" << E<< endl;
    cout << "variance=" << E2 - E*E  << endl; // <E^2> - <E>^2
    cout << "accepted=" << s.m_accepted/nCycles << endl;

    double t1 = clock();
    cout << "Time used: "<< (t1-t0)/1e6 << " s"<<endl;

}

