#include <iostream>
#include "solver.h"

using namespace std;

int main()
{
    int N = 1;
    int D = 3;
    mat r = zeros<mat>(N,D);
    r(0,0) = 1;
    cout << r;
    Solver s;
    s.m_nDimensions = 3;
    s.m_nParticles = 1;
    double E = s.localenergy(r);
    cout << E <<endl;
}

