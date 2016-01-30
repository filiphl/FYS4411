#include <iostream>
#include "solver.h"

using namespace std;

int main()
{
    int N = 2;
    int D = 3;
    mat r = zeros<mat>(N,D);
    r(0,0) = 1;
    r(1,1) = 3;
    cout << r;
    Solver s;
    s.m_nDimensions = D;
    s.m_nParticles = N;
    double E = s.localenergy(r);
    double anal = s.Analytical(r);
    cout << "Numerical: "<< E <<endl;
    cout << "Analytical: "<< anal<<endl;

    int nCycles = 15;
    for (int i=0; i<nCycles; i++){
        r = s.metropolis_step(r);

        cout<<"Total saved   " << s.localenergy(r) <<endl;
    }


}

