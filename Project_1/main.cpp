#include <iostream>
#include "solver.h"
#include <iomanip>

using namespace std;

int main()
{
    int N = 1;
    int D = 1;
    mat r = zeros<mat>(N,D);
    r(0,0) = 1;
    //r(1,1) = 0;
    cout << r;
    Solver s;
    s.m_nDimensions = D;
    s.m_nParticles = N;
    double E_sum = s.localenergy(r);
    double E_sum2 = E_sum*E_sum;
    double anal = s.Analytical(r);
    cout << "Numerical: "<< E_sum <<endl;
    cout << "Analytical: "<< anal<<endl;

    double E;
    ofstream myfile;    // Write2file
    myfile.open("energies.txt");
    int nCycles = 2000;

    double Prob;
    double localE;
    double expE = 0.0;


    for (int i=0; i<nCycles; i++){
        r = s.metropolis_step(r);

        E = s.localenergy(r);
        //cout<<"Total saved   " << E <<endl;

        Prob = s.wavefunction(r)*s.wavefunction(r);
        localE = s.localenergy(r);
        expE += localE;
        cout << localE<< endl;
        E_sum += E;
        E_sum2 += E*E;
        //myfile << E <<endl;
    }
    myfile.close();
    //system("python readdate.py");
    //cout << "average energy: " << E_sum/nCycles<<endl;

    nCycles++;
    expE = expE/nCycles;

    E = E_sum/nCycles;
    double E2 = E_sum2/nCycles;
    cout << "energy=" << E<< endl;
    cout << "variance=" << E2 - E*E  << endl; // <E^2> - <E>^2
    cout << "accepted=" << s.m_accepted/nCycles << endl;



}

