#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>

using namespace arma;


class Solver
{
private:
    double m_m;
    double m_w;
    int m_nParticles;
    int m_nDimensions;
    double m_alpha;

public:
    Solver();
    double localenergy(mat r);
    double wavefunction(mat r);

};

#endif // SOLVER_H
