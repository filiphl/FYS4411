#ifndef SLATER_H
#define SLATER_H
#include "system.h"

using namespace std;

class Slater
{
private:
    double m_R = 0;
    System* m_system = nullptr;
    vector<vector<double>> m_up;
    vector<vector<double>> m_down;
    vector<vector<double>> m_upInverse;
    vector<vector<double>> m_downInverse;

public:
    Slater(System* system);
    void update(int row, double value);


    // get / set

    vector<vector<double>>& getUp()          {return m_up;}
    vector<vector<double>>& getDown()        {return m_down;}
    vector<vector<double>>& getUpInverse()   {return m_upInverse;}
    vector<vector<double>>& getDownInverse() {return m_downInverse;}

    double getR() const;
    void setR(double R);
};

#endif // SLATER_H
