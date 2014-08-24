#include "mechanicalparameters3d.h"

MechanicalParameters3D::MechanicalParameters3D(double E, double nu)
{
    E1_ = E2_ = E3_ = E;
    nu12_ = nu21_ =
            nu13_ = nu31_ =
            nu23_ = nu32_ = nu;
    G12_ = G13_ = G23_ = E / (2.0 * (1.0 + nu));
}

MechanicalParameters3D::MechanicalParameters3D(double E, double nu, double G)
{
    E1_ = E2_ = E3_ = E;
    nu12_ = nu21_ =
            nu13_ = nu31_ =
            nu23_ = nu32_ = nu;
    G12_ = G13_ = G23_ = G;
}

MechanicalParameters3D::MechanicalParameters3D(double E1, double E2, double E3, double nu12, double nu13, double nu23, double G12, double G13, double G23)
{
    E1_ = E1;
    E2_ = E2;
    E3_ = E3;
    nu12_ = nu12;   nu21_ = nu12 * E2 / E1;
    nu13_ = nu13;   nu31_ = nu13 * E3 / E1;
    nu23_ = nu23;   nu32_ = nu23 * E3 / E2;
    G12_ = G12;
    G13_ = G13;
    G23_ = G23;
}
double MechanicalParameters3D::E1() const
{
    return E1_;
}

void MechanicalParameters3D::setE1(double E1)
{
    E1_ = E1;
}
double MechanicalParameters3D::E2() const
{
    return E2_;
}

void MechanicalParameters3D::setE2(double E2)
{
    E2_ = E2;
}
double MechanicalParameters3D::E3() const
{
    return E3_;
}

void MechanicalParameters3D::setE3(double E3)
{
    E3_ = E3;
}
double MechanicalParameters3D::nu21() const
{
    return nu21_;
}

void MechanicalParameters3D::setNu21(double nu21)
{
    nu21_ = nu21;
}
double MechanicalParameters3D::nu31() const
{
    return nu31_;
}

void MechanicalParameters3D::setNu31(double nu31)
{
    nu31_ = nu31;
}
double MechanicalParameters3D::nu12() const
{
    return nu12_;
}

void MechanicalParameters3D::setNu12(double nu12)
{
    nu12_ = nu12;
}
double MechanicalParameters3D::nu32() const
{
    return nu32_;
}

void MechanicalParameters3D::setNu32(double nu32)
{
    nu32_ = nu32;
}
double MechanicalParameters3D::nu13() const
{
    return nu13_;
}

void MechanicalParameters3D::setNu13(double nu13)
{
    nu13_ = nu13;
}
double MechanicalParameters3D::nu23() const
{
    return nu23_;
}

void MechanicalParameters3D::setNu23(double nu23)
{
    nu23_ = nu23;
}
double MechanicalParameters3D::G23() const
{
    return G23_;
}

void MechanicalParameters3D::setG23(double G23)
{
    G23_ = G23;
}
double MechanicalParameters3D::G13() const
{
    return G13_;
}

void MechanicalParameters3D::setG13(double G13)
{
    G13_ = G13;
}
double MechanicalParameters3D::G12() const
{
    return G12_;
}

void MechanicalParameters3D::setG12(double G12)
{
    G12_ = G12;
}















