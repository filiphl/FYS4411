#include "system.h"

using namespace arma;


System::System()
{

    System::waveFucntion(mat pos){
        double r = 0;
        for(int p=0; p<m_nParticles; p++){
            double singleParticle = 0;
            for (int d=0; d<m_nDimensions; d++){
                singleParticle += pos(p,d)*pos(p,d);
            }
            r += sqrt(singleParticle);
        }
        return exp(-m_alpha * r);
    }

}

