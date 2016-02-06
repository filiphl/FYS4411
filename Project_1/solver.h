#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>
#include <stdlib.h>
#include "random.h"


using namespace arma;


class Solver
{
private:
    double m_m = 1;
    double m_w = 1;
    double m_alpha = 0.5;
    int NumberOfAtoms = 0;
    int NumberOfDimensions = 0;
    double m_D = 0.5;
    double m_dt = 0.05;

public:
    Solver();
    Solver(int N, int D);
    mat r;
    mat qForceOld;
    mat qForceNew;
    int m_nParticles = 0;
    int m_nDimensions = 0;
    double m_accepted = 0;
    double m_dx = 1;
    double localenergy();
    double wavefunction();
    double placeParticles(double a);
    double Analytical();
    void metropolis_step();
    void qForce(mat &);
};

#endif // SOLVER_H
