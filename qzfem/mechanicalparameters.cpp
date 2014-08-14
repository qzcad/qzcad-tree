#include "mechanicalparameters.h"

MechanicalParameters::MechanicalParameters()
{
    E_ = 70000.0;
    nu_ = 0.3;
}

MechanicalParameters::MechanicalParameters(double E, double nu)
{
    E_ = E;
    nu_ = nu;
}

double MechanicalParameters::E() const
{
    return E_;
}

void MechanicalParameters::setE(double E)
{
    E_ = E;
}
double MechanicalParameters::nu() const
{
    return nu_;
}

void MechanicalParameters::setNu(double nu)
{
    nu_ = nu;
}

void MechanicalParameters::set(double E, double nu)
{
    E_ = E;
    nu_ = nu;
}

void MechanicalParameters::setByGK(double G, double K)
{
    nu_ = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G));
    E_ = 9.0 * K * G / (3.0 * K + G);
}



