#ifndef SOLVER_H
#define SOLVER_H
#include <armadillo>

using namespace arma;


class Solver
{
private:
    double m_m;
    int m_nParticles;
    int m_nDimensions;
    double m_alpha;

public:
    Solver();
    double localenergy(mat r);
    double wavefucntion(mat r);

};

#endif // SOLVER_H
