#ifndef SYSTEM_H
#define SYSTEM_H
#include <armadillo>
using namespace arma;

class System
{
private:
    int m_nParticles;
    int m_nDimensions;
    double m_alpha;

public:
    System();
    void hamiltonian();
    double waveFucntion(mat);

};

#endif // SYSTEM_H
