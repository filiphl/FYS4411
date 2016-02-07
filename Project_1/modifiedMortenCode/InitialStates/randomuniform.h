#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
private:
    double m_boxSize = 1;
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles);
    void setupInitialState();

    void setBoxSize(double boxSize){m_boxSize = boxSize;}
};

