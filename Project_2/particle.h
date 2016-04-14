#pragma once
#include <vector>

class Particle {

private:
    int     m_numberOfDimensions = 0;
    std::vector<double> m_oldPosition = std::vector<double>();
    std::vector<double> m_newPosition = std::vector<double>();

public:
    Particle();
    void setNumberOfDimensions(int numberOfDimensions);
    void setOldPosition(const std::vector<double> &position);
    void setNewPosition(const std::vector<double> &position);
    void adjustOldPosition(double change, int dimension);
    void adjustNewPosition(double change, int dimension);

    std::vector<double>& getNewPosition() { return m_newPosition; }
    std::vector<double>& getOldPosition() { return m_oldPosition; }
};

