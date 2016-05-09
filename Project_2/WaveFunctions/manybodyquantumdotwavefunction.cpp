#include "manybodyquantumdotwavefunction.h"

ManyBodyQuantumDotWaveFunction::ManyBodyQuantumDotWaveFunction(System *system, double omega, double a, double beta) :
    WaveFunction(system)
{
    m_omega = omega;
    m_beta = beta;
    m_npHalf =  m_system->getNumberOfParticles()/2;

    m_parameters.push_back(beta);
    m_parameters.push_back(omega);
    m_parameters.push_back(a);
    m_numberOfParameters = 3;

    m_a          = zeros<mat>(m_npHalf*2, m_npHalf*2);
    m_slaterUp   = zeros<mat>(m_npHalf, m_npHalf);
    m_slaterDown = zeros<mat>(m_npHalf, m_npHalf);
    m_distances  = zeros<mat>(m_npHalf*2, m_npHalf*2);

    m_quantumNumbers = zeros<mat>(10,2);
    m_quantumNumbers(0 ,0) = 0;  m_quantumNumbers(0, 1) = 0;
    m_quantumNumbers(1 ,0) = 1;  m_quantumNumbers(1, 1) = 0;
    m_quantumNumbers(2 ,0) = 0;  m_quantumNumbers(2, 1) = 1;
    m_quantumNumbers(3 ,0) = 2;  m_quantumNumbers(3, 1) = 0;
    m_quantumNumbers(4 ,0) = 1;  m_quantumNumbers(4, 1) = 1;
    m_quantumNumbers(5 ,0) = 0;  m_quantumNumbers(5, 1) = 2;
    m_quantumNumbers(6 ,0) = 3;  m_quantumNumbers(6, 1) = 0;
    m_quantumNumbers(7 ,0) = 2;  m_quantumNumbers(7, 1) = 1;
    m_quantumNumbers(8 ,0) = 1;  m_quantumNumbers(8, 1) = 2;
    m_quantumNumbers(9 ,0) = 0;  m_quantumNumbers(9, 1) = 3;

    for (int i=0; i<system->getNumberOfParticles(); i++){
        for (int j=0; j<system->getNumberOfParticles(); j++){
            if ((i<m_npHalf)&&(j<m_npHalf))        { m_a(i,j)=1./3; }
            else if ((i>=m_npHalf)&&(j>=m_npHalf)) { m_a(i,j)=1./3; }
            else                                   { m_a(i,j)=1.0;  }
        }
    }
    setupSlater();
}


double ManyBodyQuantumDotWaveFunction::computeLaplacian(std::vector<class Particle*> particles) {
    /*We are only interested in a single vector component, since we only move one particle and one dimension
     * at a time. Therefore, only a scalar is returned.
    double value = 0;
    for (int i=0;i<m_npHalf;i++){
        for (int j=0;j<m_npHalf;j++){
            value += ;
        }
    }*/
    double slaterLap = 0;
    double corrLap = 0;
    double crossTerm = 0;
    int k = m_newlyMoved;
    //if (k<m_npHalf){ TO BE EFFECTIFIED
    for (int i=0; i<m_npHalf; i++){

        for (int l=0; l<m_npHalf; l++){
            slaterLap += ddSingleParticleWF(i,l)*m_slaterUpInverse(l,i);
            //            cout << "ddSingleParticleWF: "<<ddSingleParticleWF(i,l)<<endl;
        }

        corrLap += correlationLap(particles, i);

        for (int d=0; d<2; d++){
            crossTerm += 2 * correlationGrad(particles, i, d) * slaterGrad(particles, i, d) * m_RSD;
        }

    }

    for (int i=m_npHalf; i<m_npHalf*2; i++){

        for (int l=m_npHalf; l<m_npHalf*2; l++){
            slaterLap += ddSingleParticleWF(i,l)*m_slaterDownInverse(l-m_npHalf,i-m_npHalf);
        }

        corrLap += correlationLap(particles, i);

        for (int d=0; d<2; d++){
            crossTerm += 2 * correlationGrad(particles, i, d) * slaterGrad(particles, i, d);
            //cout << slaterGrad(particles, i, d) << endl;
        }
    }
    /*}
    else{
        for (int i=0; i<m_npHalf*2; i++){

            for (int l=0; l<m_npHalf; l++){
                slaterLap += ddSingleParticleWF(i,l)*m_slaterUpInverse(j,i);
            }

            corrLap += correlationLap(particles, i);

            for (int d=0; d<2; d++){
                crossTerm += 2*correlationGrad(particles, i, d)*slaterGrad(particles, i, d);
            }
        }
    }*/
    //cout << crossTerm << endl;
    return corrLap + slaterLap + crossTerm;
}

double ManyBodyQuantumDotWaveFunction::computeGradient(std::vector<Particle *> particles, int particle, int dimension){

    //    cout << "SlaterGrad: "<< slaterGrad(particles, particle, dimension)<<endl;
    //    cout << "corgrad: "<< correlationGrad(particles, particle, dimension)<<endl;
    return slaterGrad(particles, particle, dimension) + correlationGrad(particles, particle, dimension);
}


double ManyBodyQuantumDotWaveFunction::evaluate(std::vector<class Particle*> particles)
{
    cout << "Entered evaluate in ManyBodyQuantumDotWaveFunction. This is strange?"<<endl;
    cout << m_distances(10000,0)<<endl;
}




double ManyBodyQuantumDotWaveFunction::SingleParticleWF(int nx, int ny, double x, double y)
{
    return  hermite(nx, x) *
            hermite(ny, y) *
            exp(-m_omega*(x*x + y*y)*0.5);
}

double ManyBodyQuantumDotWaveFunction::ddSingleParticleWF(int i, int j)
{
    int nx, ny;
    double x = m_system->getParticles()[i]->getNewPosition()[0];
    double y = m_system->getParticles()[i]->getNewPosition()[1];
    double ri2 = x*x + y*y;

    if (j<m_npHalf){    // Spin up
        nx = m_quantumNumbers(j,0);
        ny = m_quantumNumbers(j,1);
    }
    else {              // Spin down
        nx = m_quantumNumbers(j-m_npHalf,0);
        ny = m_quantumNumbers(j-m_npHalf,1);
    }

    double   Hx = hermite(nx, m_system->getParticles()[i]->getNewPosition()[0]);
    double   Hy = hermite(ny, m_system->getParticles()[i]->getNewPosition()[1]);
    double  dHx = hermiteDerivative(nx, m_system->getParticles()[i]->getNewPosition()[0]);
    double  dHy = hermiteDerivative(ny, m_system->getParticles()[i]->getNewPosition()[1]);
    double ddHx = hermiteDoubleDerivative(nx, m_system->getParticles()[i]->getNewPosition()[0]);
    double ddHy = hermiteDoubleDerivative(ny, m_system->getParticles()[i]->getNewPosition()[1]);


    return (Hy*ddHx + Hx*ddHy -
            2*m_omega*(x*Hy*dHx + y*Hx*dHy) -
            m_omega*Hx*Hy*(m_system->getNumberOfDimensions()-m_omega*ri2))*
            exp(-m_omega*ri2*0.5);
}

double ManyBodyQuantumDotWaveFunction::slaterGrad(std::vector<Particle *> particles, int k, int j)
{
    double element = 0;
    double slater = 0;


    double expRK = exp(-0.5*m_omega*(particles[k]->getNewPosition()[0]*particles[k]->getNewPosition()[0]
            + particles[k]->getNewPosition()[1]*particles[k]->getNewPosition()[1]));

    if (k < m_npHalf){
        for (int i=0; i<m_npHalf; i++){
            int nx = m_quantumNumbers(i,0);
            int ny = m_quantumNumbers(i,1);

            if (j == 0){
                element = hermite(ny, particles[i]->getNewPosition()[1])
                        *hermiteDerivative(nx, particles[i]->getNewPosition()[0]);
            }
            else{
                element = hermite(nx, particles[i]->getNewPosition()[0])
                        *hermiteDerivative(ny, particles[i]->getNewPosition()[1]);
            }
            element *= expRK;
            slater += (element - m_omega*particles[k]->getNewPosition()[j]
                       *SingleParticleWF(m_quantumNumbers(i,0),
                                         m_quantumNumbers(i,1),
                                         particles[i]->getNewPosition()[0],
                       particles[i]->getNewPosition()[1]))
                    *m_slaterUpInverse(i,k);

        }
        /* element *= expRK;
        slater += element - m_omega*particles[k]->getNewPosition()[j]
                *SingleParticleWF(m_quantumNumbers(k,0),m_quantumNumbers(k,1),particles[k]->getNewPosition()[0],particles[k]->getNewPosition()[1]);*/
    }
    else{
        for (int i=m_npHalf; i<m_system->getNumberOfParticles(); i++){
            if (j == 0){
                element = hermite(m_quantumNumbers(i-m_npHalf,1), particles[i]->getNewPosition()[1])
                        *hermiteDerivative(m_quantumNumbers(i-m_npHalf,0), particles[i]->getNewPosition()[0]);
            }
            else{
                element = hermite(m_quantumNumbers(i-m_npHalf,0), particles[i]->getNewPosition()[0])
                        *hermiteDerivative(m_quantumNumbers(i-m_npHalf,1), particles[i]->getNewPosition()[1]);
            }
            element *= expRK;
            slater += (element - m_omega*particles[k]->getNewPosition()[j]
                       *SingleParticleWF(m_quantumNumbers(i-m_npHalf,0),
                                         m_quantumNumbers(i-m_npHalf,1),
                                         particles[i]->getNewPosition()[0],
                       particles[i]->getNewPosition()[1]))
                    *m_slaterDownInverse(i-m_npHalf, k-m_npHalf);

        }
        /* element *= expRK;
        slater += element - m_omega*particles[k]->getNewPosition()[j]
                *SingleParticleWF(m_quantumNumbers(k-m_npHalf,0),m_quantumNumbers(k-m_npHalf,1),particles[k]->getNewPosition()[0],particles[k]->getNewPosition()[1]);*/
    }
    cout << "element: "<< setw(8) << element << "    slater: "<<setw(8)<<setprecision(5)<<slater<<endl;
    //cout <<slater<< endl;
    return slater/m_R;
}


double ManyBodyQuantumDotWaveFunction::correlationGrad(std::vector<Particle *> particles, int k, int d){
    //There might be an error in sign here, PLS CHECK OMGLOL
    double correlation = 0;
    for (int i=0; i<k; i++){
        double betaFrac = m_a(i,k)/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        correlation += (particles[k]->getNewPosition()[d] -
                        particles[i]->getNewPosition()[d])*
                (betaFrac/m_distances(i,k));
    }
    for (int j=k+1; j<m_npHalf*2; j++){
        double betaFrac = m_a(k,j)/((1+m_beta*m_distances(k,j))*(1+m_beta*m_distances(k,j)));
        correlation -= (particles[j]->getNewPosition()[d] -
                        particles[k]->getNewPosition()[d])*
                (betaFrac/m_distances(k,j));
    }
    //    cout << "correlation: "<<correlation<<endl;
    return correlation;
}

double ManyBodyQuantumDotWaveFunction::correlationLap(std::vector<Particle *> particles, int k)
{
    double correlation = 0;

    /*// This is needed because corrGrad calculates grad(Psi)/Psi, so (grad(Psi)/Psi)^2 must be multiplied by Psi
    double exponent = 0;
    for (int i=0; i<m_npHalf*2; i++){
        for (int j=i+1; j<m_npHalf*2; j++){
            exponent += m_a(i,j)*m_distances(i,j)/(1+m_beta*m_distances(i,j));
        }
    }
    double factor = exp(exponent);*/

    for (int d; d<2; d++){
        correlation += correlationGrad(particles,k,d)*correlationGrad(particles,k,d);
    }

    for (int i=0; i<k; i++){
        double betaFrac2 = m_a(i,k)/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        double betaFrac3 = -2*m_beta*m_a(i,k)/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        correlation += (1/m_distances(i,k))*betaFrac2 + betaFrac3;
    }

    for (int j=k+1; j<m_npHalf*2; j++){
        double betaFrac2 = m_a(k,j)/((1+m_beta*m_distances(k,j))*(1+m_beta*m_distances(k,j)));
        double betaFrac3 = -2*m_beta*m_a(k,j)/((1+m_beta*m_distances(k,j))*(1+m_beta*m_distances(k,j))*(1+m_beta*m_distances(k,j)));
        correlation += (1/m_distances(k,j))*betaFrac2 + betaFrac3;
    }

    return correlation;
}



double ManyBodyQuantumDotWaveFunction::hermite(int energyLevel, double position)
{
    double x = position;
    double sW = sqrt(m_omega);
    if (energyLevel == 0){
        return 1;
    }
    else if (energyLevel == 1) {
        return 2*x*sW;
    }
    else if (energyLevel == 2) {
        return 4*x*sW*x*sW - 2;
    }
    else if (energyLevel == 3){
        return 8*x*sW*x*sW*x*sW - 12*x*sW;
    }
    else if (energyLevel == 4){
        return 16*x*sW*x*sW*x*sW*x*sW - 48*x*sW*x*sW + 12;
    }
    else {
        cout << "Too high energy level!"<<endl;
        exit(0);
    }
}

double ManyBodyQuantumDotWaveFunction::hermiteDerivative(int energyLevel, double position)
{
    double x = position;
    double sW = sqrt(m_omega);
    if (energyLevel == 0){
        return 0;
    }
    else if (energyLevel == 1) {
        return 2*sW;
    }
    else if (energyLevel == 2) {
        return 8*x*sW*sW;
    }
    else if (energyLevel == 3){
        return 24*x*x*sW*sW*sW - 12*sW;
    }
    else if (energyLevel == 4){
        return 64*x*x*x*sW*sW*sW*sW - 96*x*sW*sW;
    }
    else {
        cout << "Too high energy level!"<<endl;
        exit(0);
    }
}

double ManyBodyQuantumDotWaveFunction::hermiteDoubleDerivative(int energyLevel, double position)
{
    double x = position;
    double sW = sqrt(m_omega);
    if (energyLevel == 0){
        return 0;
    }
    else if (energyLevel == 1) {
        return 0;
    }
    else if (energyLevel == 2) {
        return 8*sW*sW;
    }
    else if (energyLevel == 3){
        return 48*x*sW*sW*sW;
    }
    else if (energyLevel == 4){
        return 192*x*x*sW*sW*sW*sW - 96*sW*sW;
    }
    else {
        cout << "Too high energy level!"<<endl;
        exit(0);
    }
}

double ManyBodyQuantumDotWaveFunction::f(int i, int j)
{
    return m_a(i,j)*m_distances(i,j)/(1+m_beta*m_distances(i,j));
}






double ManyBodyQuantumDotWaveFunction::computeRatio(std::vector<class Particle*> particles, int i, int j, double change)
{
    //    cout << "chang: "<< change<<endl;
    particles[i]->adjustNewPosition(change, j);
    double RSD = 0;
    double RC = 0;

    double x =  particles[i]->getNewPosition()[0];
    double y =  particles[i]->getNewPosition()[1];

    if (i<m_npHalf){
        for (int k=0; k<m_npHalf; k++){
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            RSD += m_slaterUpInverse(k,i)*SingleParticleWF(nx, ny, x, y);
        }
    }
    else{
        for (int k=0; k<m_npHalf; k++){
            int nx = m_quantumNumbers(k,0);
            int ny = m_quantumNumbers(k,1);
            RSD += m_slaterDownInverse(k,i-m_npHalf)*SingleParticleWF(nx, ny, x, y);
        }

    }

    //double exponentN = 0;
    //double exponentO = 0;
    /*for(int k=0; k<m_npHalf*2; k++){
        for (int p=k+1; p<m_npHalf*2; p++){
            double r_kpN = 0;
            double r_kpO = 0;
            for (int d=0; d<m_system->getNumberOfDimensions(); d++){
                r_kpN += (particles[k]->getNewPosition()[d]-particles[p]->getNewPosition()[d])
                        *(particles[k]->getNewPosition()[d]-particles[p]->getNewPosition()[d]);
                r_kpO += (particles[k]->getOldPosition()[d]-particles[p]->getOldPosition()[d])
                        *(particles[k]->getOldPosition()[d]-particles[p]->getOldPosition()[d]);
            }
            r_kpN = sqrt(r_kpN);
            r_kpO = sqrt(r_kpO);
            exponentN += m_a(k,p)*r_kpN/(1+m_beta*r_kpN);
            exponentO += m_a(k,p)*r_kpO/(1+m_beta*r_kpO);
        }
    }*/

    double fn = 0;
    double fo = 0;
    for (int k=0; k<i; k++){
        double r_kiN = 0;
        double r_kiO = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            r_kiN += (particles[k]->getNewPosition()[d]-particles[i]->getNewPosition()[d])
                    *(particles[k]->getNewPosition()[d]-particles[i]->getNewPosition()[d]);
            r_kiO += (particles[k]->getOldPosition()[d]-particles[i]->getOldPosition()[d])
                    *(particles[k]->getOldPosition()[d]-particles[i]->getOldPosition()[d]);
        }
        r_kiN = sqrt(r_kiN);
        r_kiO = sqrt(r_kiO);
        fn += m_a(k,i)*r_kiN/(1+m_beta*r_kiN);
        fo += m_a(k,i)*r_kiO/(1+m_beta*r_kiO);
    }

    for (int k=i+1; k<m_npHalf*2; k++){
        double r_ikN = 0;
        double r_ikO = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            r_ikN += (particles[i]->getNewPosition()[d]-particles[k]->getNewPosition()[d])
                    *(particles[i]->getNewPosition()[d]-particles[k]->getNewPosition()[d]);
            r_ikO += (particles[i]->getOldPosition()[d]-particles[k]->getOldPosition()[d])
                    *(particles[i]->getOldPosition()[d]-particles[k]->getOldPosition()[d]);
        }
        r_ikN = sqrt(r_ikN);
        r_ikO = sqrt(r_ikO);
        fn += m_a(i,k)*r_ikN/(1+m_beta*r_ikN);
        fo += m_a(i,k)*r_ikO/(1+m_beta*r_ikO);
    }

    RC = exp(fn - fo);
    m_RSD = RSD;
    m_R = RSD*RC;
    // cout << RC<<endl;
    psiBeta = 0;
    for (int i=0; i<m_npHalf*2; i++){
        for (int j=i+1; j<m_npHalf*2; j++){
            psiBeta -= m_a(i,j)*m_distances(i,j)*m_distances(i,j)/((1+m_beta*m_distances(i,j))*(1+m_beta*m_distances(i,j)));
        }
    }
    return m_R;
}





void ManyBodyQuantumDotWaveFunction::updateSlater(int i)
{
    // Up
    if (i<m_npHalf){
        mat upOldInverse = m_slaterUpInverse;
        for (int k=0; k<m_npHalf; k++){
            for (int j=0; j<m_npHalf; j++){
                if (j!=i){
                    double sum = 0;
                    for (int l=0; l<m_npHalf; l++){
                        int nx = m_quantumNumbers(l,0);
                        int ny = m_quantumNumbers(l,1);
                        sum += SingleParticleWF(nx,
                                                ny,
                                                m_system->getParticles()[i]->getNewPosition()[0],
                                m_system->getParticles()[i]->getNewPosition()[1])
                                *upOldInverse(l,j);

                    }
                    m_slaterUpInverse(k,j) = upOldInverse(k,j) - (upOldInverse(k,i)/m_R)*sum;
                }
                else{
                    m_slaterUpInverse(k,j) = upOldInverse(k,j)/m_R;
                }
            }
        }
    }

    // Down
    else{
        mat downOldInverse = m_slaterDownInverse;
        for (int k=0; k<m_npHalf; k++){
            for (int j=0; j<m_npHalf; j++){
                if (j!=i-m_npHalf){
                    double sum = 0;
                    for (int l=0; l<m_npHalf; l++){
                        int nx = m_quantumNumbers(l,0);
                        int ny = m_quantumNumbers(l,1);
                        sum += SingleParticleWF(nx,
                                                ny,
                                                m_system->getParticles()[i]->getNewPosition()[0],
                                m_system->getParticles()[i]->getNewPosition()[1])
                                *downOldInverse(l,j);

                    }
                    m_slaterDownInverse(k,j) = downOldInverse(k,j) - (downOldInverse(k,i-m_npHalf)/m_R)*sum;
                }
                else{
                    m_slaterDownInverse(k,j) = downOldInverse(k,j)/m_R;
                }
            }
        }
    }

    for (int k=0; k<i;k++){
        double rki = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            rki += (m_system->getParticles()[k]->getNewPosition()[d]-m_system->getParticles()[i]->getNewPosition()[d])
                    *(m_system->getParticles()[k]->getNewPosition()[d]-m_system->getParticles()[i]->getNewPosition()[d]);
        }
        m_distances(k,i) = sqrt(rki);
    }
    for (int k=i+1; k<m_npHalf*2;k++){
        double rik = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            rik += (m_system->getParticles()[k]->getNewPosition()[d]-m_system->getParticles()[i]->getNewPosition()[d])
                    *(m_system->getParticles()[k]->getNewPosition()[d]-m_system->getParticles()[i]->getNewPosition()[d]);
        }
        m_distances(i,k) = sqrt(rik);
    }

}

void ManyBodyQuantumDotWaveFunction::printParameters()
{
    cout << " Beta  :      "<< m_beta<<endl;
    cout << " Omega :      "<< m_omega<<endl;
    cout << " a     :      "<< m_parameters[2]<<endl;
}




void ManyBodyQuantumDotWaveFunction::setupSlater()
{

    for (int i=0; i<m_npHalf; i++){
        for (int j=0; j<m_npHalf; j++){
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            //cout << "nx=" << nx << "  ny=" << ny <<endl;
            m_slaterUp(i,j)   = SingleParticleWF(nx,ny,m_system->getParticles()[i]->getOldPosition()[0], m_system->getParticles()[i]->getOldPosition()[1]);
            m_slaterDown(i,j) = SingleParticleWF(nx,ny,m_system->getParticles()[i+m_npHalf]->getOldPosition()[0], m_system->getParticles()[i+m_npHalf]->getOldPosition()[1]);
        }
    }

    for (int i=0; i<m_npHalf*2; i++){
        for (int j=i+1; j<m_npHalf*2; j++){
            double rij = 0;
            for (int d=0; d<m_system->getNumberOfDimensions(); d++){
                rij += (m_system->getParticles()[i]->getNewPosition()[d]-m_system->getParticles()[j]->getNewPosition()[d])
                        *(m_system->getParticles()[i]->getNewPosition()[d]-m_system->getParticles()[j]->getNewPosition()[d]);
            }
            m_distances(i,j) = sqrt(rij);
        }
    }

    m_slaterUpInverse   = m_slaterUp.i();
    m_slaterDownInverse = m_slaterDown.i();
    cout << "Initial slater up:\n"         << m_slaterUp << endl;
    cout << "Initial slater up inverse:\n" << m_slaterUpInverse << endl;
    cout << "Initial slater down:\n"         << m_slaterDown << endl;
    cout << "Initial slater down inverse:\n" << m_slaterDownInverse << endl;

}




