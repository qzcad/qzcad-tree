#include "fem2d.h"

#include "doublematrix.h"

Fem2D::Fem2D(Mesh *mesh, UInteger freedom, const std::list<FemCondition *> &conditions) :
    Fem(mesh, freedom, conditions)
{
}

DoubleMatrix Fem2D::evalPlaneStrainMatrix(const double &E, const double &nu)
{
    DoubleMatrix d(3);
    d(0, 0) = 1.0 - nu; d(0, 1) = nu;       d(0, 2) = 0.0;
    d(1, 0) = nu;       d(1, 1) = 1.0 - nu; d(1, 2) = 0.0;
    d(2, 0) = 0.0;      d(2, 1) = 0.0;      d(2, 2) = (1.0 - 2.0 * nu) / 2.0;
    return (E / ((1.0 + nu) * (1.0 - 2.0 * nu))) * d;
}

DoubleMatrix Fem2D::evalPlaneStrainMatrix(const double &E, const double &nu, const double &G)
{
    DoubleMatrix d = evalPlaneStrainMatrix(E, nu);
    d(2, 2) = G;
    return d;
}

DoubleMatrix Fem2D::evalPlaneStressMatrix(const double &E, const double &nu)
{
    DoubleMatrix d(3);
    d(0, 0) = 1.0;  d(0, 1) = nu;   d(0, 2) = 0.0;
    d(1, 0) = nu;   d(1, 1) = 1.0;  d(1, 2) = 0.0;
    d(2, 0) = 0.0;  d(2, 1) = 0.0;  d(2, 2) = (1 - nu) / 2.0;
    return (E / (1.0 - nu*nu)) * d;
}

DoubleMatrix Fem2D::evalPlaneStressMatrix(const double &E, const double &nu, const double &G)
{
    DoubleMatrix d = evalPlaneStressMatrix(E, nu);
    d(2, 2) = G;
    return d;
}

double Fem2D::isoQuad4(const double &xi, const double &eta, double x[], double y[], DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY)
{
    const unsigned int elementNodes = 4;
    DoubleVector dNdXi(elementNodes);
    DoubleVector dNdEta(elementNodes);
    DoubleMatrix Jacobi(2, 0.0); // матрица Якоби
    // билинейные функции формы
    N(0) = (1.0 - xi) * (1.0 - eta) / 4.0;
    N(1) = (1.0 + xi) * (1.0 - eta) / 4.0;
    N(2) = (1.0 + xi) * (1.0 + eta) / 4.0;
    N(3) = (1.0 - xi) * (1.0 + eta) / 4.0;
    // производные функций формы
    dNdXi(0) = -(1.0 - eta) / 4.0;
    dNdXi(1) = (1.0 - eta) / 4.0;
    dNdXi(2) = (1.0 + eta) / 4.0;
    dNdXi(3) = -(1.0 + eta) / 4.0;

    dNdEta(0) = -(1.0 - xi) / 4.0;
    dNdEta(1) = -(1.0 + xi) / 4.0;
    dNdEta(2) = (1.0 + xi) / 4.0;
    dNdEta(3) = (1.0 - xi) / 4.0;
    // матрица Якоби
    for (unsigned int i = 0; i < elementNodes; i++)
    {
        Jacobi(0, 0) += dNdXi(i) * x[i];  Jacobi(0, 1) += dNdXi(i) * y[i];
        Jacobi(1, 0) += dNdEta(i) * x[i]; Jacobi(1, 1) += dNdEta(i) * y[i];
    }
    // якобиан
    double jacobian = Jacobi.det2x2();
    // обратный якобиан
    DoubleMatrix invJacobi = Jacobi.inverted2x2();
    // производные в глобальных координатах
    for (unsigned int i = 0; i < elementNodes; i++)
    {
        dNdX(i) = invJacobi(0, 0) * dNdXi(i) + invJacobi(0, 1) * dNdEta(i);
        dNdY(i) = invJacobi(1, 0) * dNdXi(i) + invJacobi(1, 1) * dNdEta(i);
    }
    return jacobian;
}

double Fem2D::isoTriangle3(const double &xi, const double &eta, double x[], double y[], DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY)
{
    const unsigned int elementNodes = 3;
    DoubleVector dNdXi(elementNodes);
    DoubleVector dNdEta(elementNodes);
    DoubleMatrix Jacobi(2, 0.0); // матрица Якоби
    // линейные функции формы
    N(0) = 1.0 - xi - eta;
    N(1) = xi;
    N(2) = eta;
    // производные функций формы
    dNdXi(0) = -1.0;
    dNdXi(1) = 1.0;
    dNdXi(2) = 0.0;

    dNdEta(0) = -1.0;
    dNdEta(1) = 0.0;
    dNdEta(2) = 1.0;
    for (unsigned int i = 0; i < elementNodes; i++)
    {
        Jacobi(0, 0) += dNdXi(i) * x[i];  Jacobi(0, 1) += dNdXi(i) * y[i];
        Jacobi(1, 0) += dNdEta(i) * x[i]; Jacobi(1, 1) += dNdEta(i) * y[i];
    }
    // якобиан
    double jacobian = Jacobi.det2x2();
    // обратный якобиан
    DoubleMatrix invJacobi = Jacobi.inverted2x2();
    // производные в глобальных координатах
    for (unsigned int i = 0; i < elementNodes; i++)
    {
        dNdX(i) = invJacobi(0, 0) * dNdXi(i) + invJacobi(0, 1) * dNdEta(i);
        dNdY(i) = invJacobi(1, 0) * dNdXi(i) + invJacobi(1, 1) * dNdEta(i);
    }
    return jacobian;
}
