#include "solver.h"
using namespace arma;

Solver::Solver(){
    m_dx = 1;
}

Solver::Solver(int N, int D){
    m_nParticles = N;
    m_nDimensions = D;
    r = zeros<mat>(N,D);
    qForceOld = zeros<mat>(N, D);
    qForceNew = zeros<mat>(N, D);
    placeParticles(0.5);
}



double Solver::Analytical(){
    double energy = 0;
    double h2 = 1;
    for (int p=0; p<m_nParticles; p++){
        double r_single = 0;
        for (int d=0; d<m_nDimensions; d++){
            r_single += r(p,d)*r(p,d);
        }
        energy += h2*m_alpha*(-2*m_alpha*r_single + m_nDimensions) + 0.5*m_m*m_w*m_w*r_single;
    }
    return energy;
}


double Solver::localenergy(){
    double kinetic = 0;
    double potential = 0;
    double hbar = 1;
    double dstep = 0.0001;
    mat rplus = r;
    mat rminus = r;
    double psiplus;
    double psiminus;
    double psi = wavefunction();


    for (int p = 0; p<m_nParticles; p++){
        double r_single = 0;
        for (int d=0; d<m_nDimensions; d++){
            r(p,d) += dstep;
            psiplus = wavefunction();
            r(p,d) -= 2*dstep;
            psiminus = wavefunction();
            r(p,d) += dstep;

            /*rplus(p,d) += dstep;
            rminus(p,d) -= dstep;

            psiplus = wavefunction(rplus);
            psiminus = wavefunction(rminus);
            */
            kinetic += (psiplus -2*psi + psiminus);

            r_single += r(p,d)*r(p,d);

            rplus(p,d) = r(p,d);
            rminus(p,d) = r(p,d);
        }
        potential += r_single;
    }


    kinetic *= -0.5*(hbar*hbar)/(m_m*dstep*dstep*psi);
    potential *= 0.5*m_m*m_w*m_w;


    return kinetic +potential;
}

double Solver::wavefunction(){
    double psi = 1;
    for (int p=0; p<m_nParticles; p++){
        double r_single = 0;
        double g = 0;
        for (int d=0; d<m_nDimensions; d++){
            r_single += r(p,d)*r(p,d);
        }
        //g = exp(-m_alpha * r_single );
        //psi *= g;
        g = -m_alpha * r_single;
        psi += g;
    }
    return exp(psi);
}

double Solver::placeParticles(double a)
{
    for (int p=0; p<m_nParticles; p++){
        for (int d=0; d<m_nDimensions; d++){
            if (Random::nextDouble()>0.5){
                r(p,d) += Random::nextDouble()*a;
            }
            else{
                r(p,d) -= Random::nextDouble()*a;
            }

        }
    }
}


void Solver::metropolis_step(){

    int p = Random::nextInt(m_nParticles);
    double oldwavefunction = wavefunction();

    double* dx = new double[3];
    qForce(qForceOld);

    for (int d=0; d<m_nDimensions; d++){
        dx[d] = m_dx*Random::nextGaussian(0,0.2);
        r(p,d) += dx[d];
    }

    double newwavefunction = wavefunction();
    qForce(qForceNew);

    double GreensFunction = 0.0;
    for (int p=0; p < m_nParticles; p++) {
        for (int d=0; d< m_nDimensions; d++){
            GreensFunction += 0.5*(qForceOld(p,d)+qForceNew(p,d))*
                    (m_D*m_dt*0.5*(-qForceOld(p,d)+qForceNew(p,d))-r(p,d)+dx[d]);
        }
    }
    GreensFunction = exp(GreensFunction);

    double prob = GreensFunction*(newwavefunction*newwavefunction)/(oldwavefunction*oldwavefunction);
    double mynt = Random::nextDouble(); // Uniform [0,1]

    if (mynt < prob){
        m_accepted++;
    }
    else{
        for (int d=0; d<m_nDimensions; d++){
            //newr(p,d) = oldr(p,d);
            r(p,d) -= dx[d];
        }
    }

}

void Solver::qForce(mat &qForceX)
{
    for (int p=0; p<m_nParticles; p++){
        for (int d=0; d<m_nDimensions; d++){
            r(p,d) -= m_dx;
            double waveFunctionOld = wavefunction();
            r(p,d) += 2*m_dx;
            double waveFunctionNew = wavefunction();
            r(p,d) -= m_dx;
            qForceX(p,d) = (waveFunctionNew - waveFunctionOld)/m_dx;    //factor 2 cancels with factor 2 in quantum force
        }
    }
    qForceX /= wavefunction();
}

