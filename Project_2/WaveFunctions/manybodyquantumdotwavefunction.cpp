#include "manybodyquantumdotwavefunction.h"

ManyBodyQuantumDotWaveFunction::ManyBodyQuantumDotWaveFunction(System *system, double alpha, double omega, int a, double beta) :
    WaveFunction(system)
{
    m_np     = m_system->getNumberOfParticles();
    m_npHalf = m_np/2;

    if (m_np!=2 && m_np!=6 && m_np!=12 && m_np!=20){
        cout << "The Many body quantum dot is only implemented for 2, 6, 12 and 20 particles."<<endl;
        cout << "Number of particles given: "<< m_np<<endl;
        exit(0);
    }



    if (system->getNumberOfDimensions() != 2){
        cout << "The Many body quantum dot is only implemented for 2 dimensions." << endl;
        cout << "Number of dimensions given: "<< system->getNumberOfDimensions() << endl;
        exit(0);
    }

    name     = "Many body quantum dot";
    setAlpha(alpha);
    setOmega(omega);
    setBeta(beta);
    m_oa  = omega*alpha;
    m_npHalf =  m_system->getNumberOfParticles()/2;
    m_system->Jastrow = a;


    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(omega);
    m_parameters.push_back(a);
    m_numberOfParameters = 4;

    m_a          = zeros<mat>(m_np, m_np);
    m_slaterUp   = zeros<mat>(m_npHalf, m_npHalf);
    m_slaterDown = zeros<mat>(m_npHalf, m_npHalf);
    m_distances  = zeros<mat>(m_np, m_np);

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

    if (a != 0){//this test is nice for when we don't want interactions
        for (int i=0; i<system->getNumberOfParticles(); i++){
            for (int j=0; j<system->getNumberOfParticles(); j++){
                if ((i<m_npHalf)&&(j<m_npHalf))        { m_a(i,j)=1./3; }
                else if ((i>=m_npHalf)&&(j>=m_npHalf)) { m_a(i,j)=1./3; }
                else                                   { m_a(i,j)=1.0;  }
            }
        }
    }
    else{
        for (int i=0; i<system->getNumberOfParticles(); i++){
            for (int j=0; j<system->getNumberOfParticles(); j++){ m_a(i,j) = 0; }
        }
    }
    setupSlater();
}


double ManyBodyQuantumDotWaveFunction::evaluate(std::vector<class Particle*> particles)
{
    /*
    cout << "Entered evaluate in ManyBodyQuantumDotWaveFunction. This is strange?"<<endl;
    cout << m_distances(10000,0)<<endl;
    */
    mat slaterUp = zeros<mat>(m_npHalf, m_npHalf);
    mat slaterDown = zeros<mat>(m_npHalf, m_npHalf);

    for (int i=0; i<m_npHalf; i++){
        for (int j=0; j<m_npHalf; j++){
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double x =  particles[i]->getNewPosition()[0];
            double y =  particles[i]->getNewPosition()[1];
            slaterUp(i,j) = SingleParticleWF(nx,ny,x,y);
        }
    }

    for (int i=0; i<m_npHalf; i++){
        for (int j=0; j<m_npHalf; j++){
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            double x =  particles[i+m_npHalf]->getNewPosition()[0];
            double y =  particles[i+m_npHalf]->getNewPosition()[1];
            slaterDown(i,j) = SingleParticleWF(nx,ny,x,y);
        }
    }


    double exponent = 0;
    for (int i=0; i<m_np; i++){
        for (int j=i+1; j<m_np; j++){
/*
            double testr12 = 0;
            for (int k=0; k<2; k++){
                testr12 += (particles[i]->getNewPosition()[k]-particles[j]->getNewPosition()[k])*
                           (particles[i]->getNewPosition()[k]-particles[j]->getNewPosition()[k]);
            }
            testr12 = sqrt(testr12);
            cout << setprecision(10)<<testr12 << "    " << setprecision(10)<< m_distances(i,j)<< "    Diff: "<<testr12-m_distances(i,j) <<endl;
*/
            exponent += m_a(i,j)*m_distances(i,j)/(1+m_beta*m_distances(i,j));
        }
    }
    double correlation = exp(exponent);

    double detUp =   det(slaterUp);
    double detDown = det(slaterDown);
    //cout << correlation<<endl;
    return detUp*detDown*correlation;
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
    double ddr = 0;

    if (m_system->getAnalyticalLaplacian()){
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
                crossTerm += 2 * correlationGrad(particles, i, d) * slaterGrad(particles, i, d) * m_RSD;        // FILIP
            }

        }

        for (int i=m_npHalf; i<m_np; i++){

            for (int l=m_npHalf; l<m_np; l++){
                slaterLap += ddSingleParticleWF(i,l)*m_slaterDownInverse(l-m_npHalf,i-m_npHalf);
            }

            corrLap += correlationLap(particles, i);

            for (int d=0; d<2; d++){
                crossTerm += 2 * correlationGrad(particles, i, d) * slaterGrad(particles, i, d) * m_RSD;        // FILIP
            }
        }

        ddr = corrLap + slaterLap + crossTerm;
    }

    else{

        m_derivativeStepLength = 1e-5;
        double psi = evaluate( particles );
        for (int i=0; i<m_system->getNumberOfParticles(); i++){
            for (int j=0; j<m_system->getNumberOfDimensions(); j++){
                particles[i]->adjustNewPosition( m_derivativeStepLength, j );      // +
                updateDistances(i);
                double psiPlus  =   evaluate( particles );
                particles[i]->adjustNewPosition( -2 * m_derivativeStepLength, j ); // -
                updateDistances(i);
                double psiMinus =   evaluate( particles );
                particles[i]->adjustNewPosition( m_derivativeStepLength, j );      // reset
                updateDistances(i);
                ddr += psiPlus - 2*psi + psiMinus;
            }
        }
        ddr = ddr / (m_derivativeStepLength * m_derivativeStepLength);
    }
    return ddr;
}






double ManyBodyQuantumDotWaveFunction::computeGradient(std::vector<Particle *>& particles, int particle, int dimension)
{
//    double dr=0;
//    m_derivativeStepLength = 1e-5;
//    for (int i=0; i<m_system->getNumberOfParticles(); i++){
//        for (int j=0; j<m_system->getNumberOfDimensions(); j++){
//            particles[i]->adjustNewPosition( m_derivativeStepLength, j );      // +
//            double psiPlus  =   evaluate( particles );
//            particles[i]->adjustNewPosition( -2 * m_derivativeStepLength, j ); // -
//            double psiMinus =   evaluate( particles );
//            particles[i]->adjustNewPosition( m_derivativeStepLength, j );      // reset
//            dr += psiPlus - psiMinus;
//        }
//    }
//    return dr/(2*m_derivativeStepLength);
    return slaterGrad(particles, particle, dimension) + correlationGrad(particles, particle, dimension);
}






double ManyBodyQuantumDotWaveFunction::SingleParticleWF(int nx, int ny, double x, double y)
{
    return  hermite(nx, x) *
            hermite(ny, y) *
            exp(-m_oa*(x*x + y*y)*0.5);
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

    double   Hx = hermite                       (nx, m_system->getParticles()[i]->getNewPosition()[0]);
    double   Hy = hermite                       (ny, m_system->getParticles()[i]->getNewPosition()[1]);
    double  dHx = hermiteDerivative             (nx, m_system->getParticles()[i]->getNewPosition()[0]);
    double  dHy = hermiteDerivative             (ny, m_system->getParticles()[i]->getNewPosition()[1]);
    double ddHx = hermiteDoubleDerivative       (nx, m_system->getParticles()[i]->getNewPosition()[0]);
    double ddHy = hermiteDoubleDerivative       (ny, m_system->getParticles()[i]->getNewPosition()[1]);


    return (Hy*ddHx + Hx*ddHy -
            2*m_oa*(x*Hy*dHx + y*Hx*dHy) -
            m_oa*Hx*Hy*(m_system->getNumberOfDimensions()-m_oa*ri2))*
            exp(-m_oa*ri2*0.5);
}




double ManyBodyQuantumDotWaveFunction::slaterGrad(std::vector<Particle *>& particles, int k, int j)
{
    double element = 0;
    double slater = 0;

    double expRK = exp(-0.5*m_oa*( particles[k]->getNewPosition()[0] * particles[k]->getNewPosition()[0]+
                                      particles[k]->getNewPosition()[1] * particles[k]->getNewPosition()[1]) );

    if (k < m_npHalf){
        for (int i=0; i<m_npHalf; i++){
            int nx = m_quantumNumbers(i,0);
            int ny = m_quantumNumbers(i,1);

            if (j == 0){
                element = hermite(ny, particles[k]->getNewPosition()[1])
                         *hermiteDerivative(nx, particles[k]->getNewPosition()[0]);
            }
            else{
                element = hermite(nx, particles[k]->getNewPosition()[0])
                         *hermiteDerivative(ny, particles[k]->getNewPosition()[1]);
            }
            element *= expRK;
            slater += (element - m_oa*particles[k]->getNewPosition()[j]
                       *SingleParticleWF(m_quantumNumbers(i,0),
                                         m_quantumNumbers(i,1),
                                         particles[k]->getNewPosition()[0],
                                         particles[k]->getNewPosition()[1]))
                                         *m_slaterUpInverse(i,k);
        }
        /* element *= expRK;
        slater += element - m_oa*particles[k]->getNewPosition()[j]
                *SingleParticleWF(m_quantumNumbers(k,0),m_quantumNumbers(k,1),particles[k]->getNewPosition()[0],particles[k]->getNewPosition()[1]);*/
    }
    else{
        for (int i=m_npHalf; i<m_system->getNumberOfParticles(); i++){
            int nx = m_quantumNumbers(i-m_npHalf,0);
            int ny = m_quantumNumbers(i-m_npHalf,1);

            if (j == 0){
                element = hermite(ny, particles[k]->getNewPosition()[1])
                         *hermiteDerivative(nx, particles[k]->getNewPosition()[0]);
            }
            else{
                element = hermite(nx, particles[k]->getNewPosition()[0])
                         *hermiteDerivative(ny, particles[k]->getNewPosition()[1]);
            }
            element *= expRK;
            slater += (element - m_oa*particles[k]->getNewPosition()[j]
                       *SingleParticleWF(m_quantumNumbers(i-m_npHalf,0),
                                         m_quantumNumbers(i-m_npHalf,1),
                                         particles[k]->getNewPosition()[0],
                                         particles[k]->getNewPosition()[1]))
                                         *m_slaterDownInverse(i-m_npHalf, k-m_npHalf);

        }
        /*
        element *= expRK;
        slater += element - m_oa*particles[k]->getNewPosition()[j]
                *SingleParticleWF(m_quantumNumbers(k-m_npHalf,0),m_quantumNumbers(k-m_npHalf,1),particles[k]->getNewPosition()[0],particles[k]->getNewPosition()[1]);
*/
    }
    //cout << "element: "<< setw(8) << element << "    slater: "<<setw(8)<<setprecision(5)<<slater<<endl;
    //cout <<slater<< endl;
    return slater/m_RSD;    // FILIP
}


double ManyBodyQuantumDotWaveFunction::correlationGrad(std::vector<Particle *>& particles, int k, int d){
    //There might be an error in sign here, PLS CHECK OMGLOL
    double correlation = 0;
    for (int i=0; i<k; i++){
        double betaFrac = m_a(i,k)/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        correlation += (particles[k]->getNewPosition()[d] -
                        particles[i]->getNewPosition()[d])*
                        (betaFrac/m_distances(i,k));
    }
    for (int j=k+1; j<m_np; j++){
        double betaFrac = m_a(k,j)/((1+m_beta*m_distances(k,j))*(1+m_beta*m_distances(k,j)));
        correlation -= (particles[j]->getNewPosition()[d] -
                        particles[k]->getNewPosition()[d])*
                        (betaFrac/m_distances(k,j));
    }
    //cout << "correlation: "<<correlation<<endl;
    return correlation;
}

double ManyBodyQuantumDotWaveFunction::correlationLap(std::vector<Particle *>& particles, int k)
{
    double correlation = 0;

    for (int d; d<2; d++){
        correlation += correlationGrad(particles,k,d)*correlationGrad(particles,k,d);
    }

    for (int i=0; i<k; i++){
        double betaFrac2 = m_a(i,k)/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        double betaFrac3 = -2*m_a(i,k)*m_beta/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        correlation += (1/m_distances(i,k))*betaFrac2 + betaFrac3;
    }

    for (int i=k+1; i<m_np; i++){
        double betaFrac2 = m_a(i,k)/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        double betaFrac3 = -2*m_a(i,k)*m_beta/((1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k))*(1+m_beta*m_distances(i,k)));
        correlation += (1/m_distances(i,k))*betaFrac2 + betaFrac3;
    }

    return correlation;
}



double ManyBodyQuantumDotWaveFunction::hermite(int energyLevel, double position)
{
    double x = position;
    double sW = sqrt(m_oa);
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
    double sW = sqrt(m_oa);
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
    double sW = sqrt(m_oa);
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
    particles[i]->adjustNewPosition(change, j);
    updateDistances(i);

    double RSD = 0;
    double RC = 0;

    double x =  particles[i]->getNewPosition()[0];
    double y =  particles[i]->getNewPosition()[1];

    if (i<m_npHalf){
        for (int j=0; j<m_npHalf; j++){
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            RSD += m_slaterUpInverse(j,i)*SingleParticleWF(nx, ny, x, y);
        }
    }
    else{
        for (int j=0; j<m_npHalf; j++){
            int nx = m_quantumNumbers(j,0);
            int ny = m_quantumNumbers(j,1);
            RSD += m_slaterDownInverse(j,i-m_npHalf)*SingleParticleWF(nx, ny, x, y);
        }

    }

    double fn = 0;
    double fo = 0;

    for (int k=0; k<i; k++){
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

    for (int k=i+1; k<2*m_npHalf; k++){
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
    //cout << m_a << endl;

    RC = exp(fn - fo);
    m_RSD = RSD;
    m_R = RSD*RC;
    //cout << "RC: "<<RC<<"   RSD: "<<RSD<<"   R: "<<m_R<<endl;
    psiBeta = 0;
    for (int i=0; i<m_np; i++){
        for (int j=i+1; j<m_np; j++){
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
                    m_slaterUpInverse(k,j) = upOldInverse(k,j) - (upOldInverse(k,i)/m_RSD)*sum;
                }
                else{
                    m_slaterUpInverse(k,j) = upOldInverse(k,i)/m_RSD;
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
                    m_slaterDownInverse(k,j) = downOldInverse(k,j) - (downOldInverse(k,i-m_npHalf)/m_RSD)*sum;
                }
                else{
                    m_slaterDownInverse(k,j) = downOldInverse(k,i-m_npHalf)/m_RSD;
                }
            }
        }
    }


    if (m_system->OptimizingParameters()){
        psiAlpha = 0;
        psiBeta  = 0;
        double dAijUp   = 0;
        double dAijDown = 0;
        for (int i=0; i<m_npHalf; i++){
            double xiUp = m_system->getParticles()[i]->getNewPosition()[0];
            double yiUp = m_system->getParticles()[i]->getNewPosition()[1];
            double xiDown = m_system->getParticles()[i+m_npHalf]->getNewPosition()[0];
            double yiDown = m_system->getParticles()[i+m_npHalf]->getNewPosition()[1];
            for (int j=0; j<m_npHalf; j++){

                int nx = m_quantumNumbers(j,0);
                int ny = m_quantumNumbers(j,1);

                dAijUp = (hermiteDerivativeAlpha(nx,xiUp)*hermite(ny,yiUp)
                          + hermiteDerivativeAlpha(ny,yiUp)*hermite(nx,xiUp)
                          - hermite(nx,xiUp)*hermite(ny,yiUp)*m_oa*(xiUp*xiUp + yiUp*yiUp)/(2*m_alpha))
                         * m_slaterUpInverse(i,j);

                dAijDown = (hermiteDerivativeAlpha(nx,xiDown)*hermite(ny,yiDown)
                            + hermiteDerivativeAlpha(ny,yiDown)*hermite(nx,xiDown)
                            - hermite(nx,xiDown)*hermite(ny,yiDown)*m_oa*(xiDown*xiDown + yiDown*yiDown)/(2*m_alpha))
                           * m_slaterDownInverse(i,j);

                // Dividing by alpha because we defined m_oa=omega*alpha in constructer.

            psiAlpha += dAijUp*exp(-m_oa*(xiUp*xiUp + yiUp*yiUp));
            psiAlpha += dAijDown*exp(-m_oa*(xiDown*xiDown + yiDown*yiDown));
            }
        }


        for (int i=0; i<m_np; i++){
            for (int j=i+1; j<m_np;j++){
                psiBeta -= m_a(i,j) * (m_distances(i,j)/(1+m_beta*m_distances(i,j))) * (m_distances(i,j)/(1+m_beta*m_distances(i,j)));
            }
        }
    }
}


double ManyBodyQuantumDotWaveFunction::hermiteDerivativeAlpha(int energyLevel, double position){
    double x = position;
    double sW = sqrt(m_oa);  // m_oa = omega*alpha. See constructor.
    if (energyLevel == 0){
        return 0;
    }
    else if (energyLevel == 1) {
        return x*sW/m_alpha;
    }
    else if (energyLevel == 2) {
        return 4*x*x*sW*sW/m_alpha;
    }
    else if (energyLevel == 3){
        return 12*x*x*x*m_oa*m_oa/(sW*m_alpha) - 6*sW*x/m_alpha;
    }
    else {
        cout << energyLevel<<endl;
        cout << "Way too high energy level!"<<endl;
        exit(0);
    }
}

void ManyBodyQuantumDotWaveFunction::updateDistances(int i)
{
    for (int k=0; k<i;k++){
        double rki = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            rki +=   (m_system->getParticles()[k]->getNewPosition()[d] - m_system->getParticles()[i]->getNewPosition()[d])
                    *(m_system->getParticles()[k]->getNewPosition()[d] - m_system->getParticles()[i]->getNewPosition()[d]);
        }
        m_distances(k,i) = sqrt(rki);
        m_distances(i,k) = sqrt(rki);   // can be removed
    }
    for (int k=i+1; k<m_np;k++){
        double rki = 0;
        for (int d=0; d<m_system->getNumberOfDimensions(); d++){
            rki +=   (m_system->getParticles()[k]->getNewPosition()[d] - m_system->getParticles()[i]->getNewPosition()[d])
                    *(m_system->getParticles()[k]->getNewPosition()[d] - m_system->getParticles()[i]->getNewPosition()[d]);
        }
        m_distances(i,k) = sqrt(rki);
        m_distances(k,i) = sqrt(rki);   // can be removed
    }
}












void ManyBodyQuantumDotWaveFunction::printParameters()
{
    cout << setw(25) << setprecision(10) << " Alpha :      "<< m_alpha<<endl;
    cout << setw(25) << setprecision(10) << " Beta  :      "<< m_beta <<endl;
    cout << setw(25) << setprecision(5)  << " Omega :      "<< m_omega<<endl;
    //cout << setw(25) << setprecision(5) << " a     :      "<< m_parameters[2]<<endl;
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

    for (int i=0; i<m_np; i++){
        for (int j=i+1; j<m_np; j++){
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
    /*
    cout << "Initial slater up:\n"         << m_slaterUp << endl;
    cout << "Initial slater up inverse:\n" << m_slaterUpInverse << endl;
    cout << "Initial slater down:\n"         << m_slaterDown << endl;
    cout << "Initial slater down inverse:\n" << m_slaterDownInverse << endl;
    */
}




