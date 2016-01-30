#include "solver.h"

Solver::Solver()
{
    Solver::localenergy(mat r){

    }

    Solver::wavefunction(mat r){
        double psi = 1;
        for (int p=0; p<m_nParticles; p++){
            double r_single = 0;
            double g = 0;
            for (int d=0; d<m_nDimensions; d++){
                r_single += r(p,d) * r(p,d);
            }
            g = exp(-m_alpha * r_single );
            double psi = psi*g;
        total += r_single;
        }
        return exp(-m_alpha * total )
    }
}

