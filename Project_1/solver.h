#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>

using namespace arma;


class Solver
{
private:
    double m_m = 1;
    double m_w = 1;
    double m_alpha = 1;

public:
    int m_nParticles = 0;
    int m_nDimensions = 0;
    Solver();
    double localenergy(mat);
    double wavefunction(mat);
    void addparticle();
    double Analytical(mat);

};

#endif // SOLVER_H
