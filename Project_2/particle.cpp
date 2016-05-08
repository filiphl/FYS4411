#include "particle.h"
#include <cassert>
#include <iostream>

using namespace std;
Particle::Particle() {
}

void Particle::setOldPosition(const std::vector<double> &position) {
    assert(position.size() == m_numberOfDimensions);
    m_oldPosition = position;
}

void Particle::setNewPosition(const std::vector<double> &position)
{
    assert(position.size() == m_numberOfDimensions);
    m_newPosition = position;
}

void Particle::adjustOldPosition(double change, int dimension) {
    m_oldPosition[dimension] += change;
}

void Particle::adjustNewPosition(double change, int dimension) {
    m_newPosition[dimension] += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
