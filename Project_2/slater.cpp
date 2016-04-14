#include "slater.h"


double Slater::getR() const
{
    return m_R;
}

void Slater::setR(double R)
{
    m_R = R;
}

Slater::Slater(System *system)
{
    m_system = system;

    // Initialize determinants
    int height = m_system->getNumberOfParticles()/2;
    m_up.resize(height);
    m_down.resize(height);
    m_upInverse.resize(height);
    m_downInverse.resize(height);
    for (int i = 0; i<height; i++){
        m_up[i].resize(height);
        m_down[i].resize(height);
        m_upInverse[i].resize(height);
        m_downInverse[i].resize(height);
        for (int j=0; j<height; j++){
            m_up          [i][j] = 0;
            m_down        [i][j] = 0;
            m_upInverse   [i][j] = 0;
            m_downInverse [i][j] = 0;
        }
    }
}

void Slater::update(int row, double value)
{
    cout << "The algorithms for updating the determinants are not implemented yet."<<endl;
    exit(0);
}



