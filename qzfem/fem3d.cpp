#include "fem3d.h"

Fem3D::Fem3D(Mesh *mesh, UInteger freedom, const std::list<FemCondition *> &conditions) : Fem(mesh, freedom, conditions)
{
}

double Fem3D::isoHex8(const double &xi, const double &eta, const double &mu, const DoubleVector &x, const DoubleVector &y, const DoubleVector &z, DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY, DoubleVector &dNdZ)
{
    const mtx::size_type elementNodes = 8;
    // значения производных функций формы в местных координатах
    DoubleVector dNdXi(elementNodes);
    DoubleVector dNdEta(elementNodes);
    DoubleVector dNdMu(elementNodes);
    // Якобиан
    DoubleMatrix Jacobi(3, 0.0);
    DoubleMatrix invJacobian(3);
    double jacobian;
    N(0) = (1 - xi) * (1 - eta) * (1 - mu) / 8.0;
    N(1) = (1 - xi) * (1 - eta) * (1 + mu) / 8.0;
    N(2) = (1 + xi) * (1 - eta) * (1 + mu) / 8.0;
    N(3) = (1 + xi) * (1 - eta) * (1 - mu) / 8.0;
    N(4) = (1 - xi) * (1 + eta) * (1 - mu) / 8.0;
    N(5) = (1 - xi) * (1 + eta) * (1 + mu) / 8.0;
    N(6) = (1 + xi) * (1 + eta) * (1 + mu) / 8.0;
    N(7) = (1 + xi) * (1 + eta) * (1 - mu) / 8.0;

    dNdXi(0) = -(1 - eta) * (1 - mu) / 8.0;
    dNdXi(1) = -(1 - eta) * (1 + mu) / 8.0;
    dNdXi(2) = (1 - eta) * (1 + mu) / 8.0;
    dNdXi(3) = (1 - eta) * (1 - mu) / 8.0;
    dNdXi(4) = -(1 + eta) * (1 - mu) / 8.0;
    dNdXi(5) = -(1 + eta) * (1 + mu) / 8.0;
    dNdXi(6) = (1 + eta) * (1 + mu) / 8.0;
    dNdXi(7) = (1 + eta) * (1 - mu) / 8.0;

    dNdEta(0) = -(1 - xi) * (1 - mu) / 8.0;
    dNdEta(1) = -(1 - xi) * (1 + mu) / 8.0;
    dNdEta(2) = -(1 + xi) * (1 + mu) / 8.0;
    dNdEta(3) = -(1 + xi) * (1 - mu) / 8.0;
    dNdEta(4) = (1 - xi) * (1 - mu) / 8.0;
    dNdEta(5) = (1 - xi) * (1 + mu) / 8.0;
    dNdEta(6) = (1 + xi) * (1 + mu) / 8.0;
    dNdEta(7) = (1 + xi) * (1 - mu) / 8.0;

    dNdMu(0) = -(1 - xi) * (1 - eta) / 8.0;
    dNdMu(1) = (1 - xi) * (1 - eta) / 8.0;
    dNdMu(2) = (1 + xi) * (1 - eta) / 8.0;
    dNdMu(3) = -(1 + xi) * (1 - eta) / 8.0;
    dNdMu(4) = -(1 - xi) * (1 + eta) / 8.0;
    dNdMu(5) = (1 - xi) * (1 + eta) / 8.0;
    dNdMu(6) = (1 + xi) * (1 + eta) / 8.0;
    dNdMu(7) = -(1 + xi) * (1 + eta) / 8.0;

    Jacobi(0, 0) = 0.0; Jacobi(0, 1) = 0.0; Jacobi(0, 2) = 0.0;
    Jacobi(1, 0) = 0.0; Jacobi(1, 1) = 0.0; Jacobi(1, 2) = 0.0;
    Jacobi(2, 0) = 0.0; Jacobi(2, 1) = 0.0; Jacobi(2, 2) = 0.0;

    for (mtx::size_type i = 0; i < elementNodes; i++)
    {
        Jacobi(0, 0) += dNdXi(i) * x[i]; Jacobi(0, 1) += dNdXi(i) * y[i]; Jacobi(0, 2) += dNdXi(i) * z[i];
        Jacobi(1, 0) += dNdEta(i) * x[i]; Jacobi(1, 1) += dNdEta(i) * y[i]; Jacobi(1, 2) += dNdEta(i) * z[i];
        Jacobi(2, 0) += dNdMu(i) * x[i]; Jacobi(2, 1) += dNdMu(i) * y[i]; Jacobi(2, 2) += dNdMu(i) * z[i];
    }
    invJacobian = Jacobi.inverted3x3();
    jacobian = Jacobi.det3x3();

    for (mtx::size_type i = 0; i < elementNodes; i++)
    {
        dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
        dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
        dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
    }
    return jacobian;
}

double Fem3D::isoTet4(const double &xi, const double &eta, const double &mu, const DoubleVector &x, const DoubleVector &y, const DoubleVector &z, DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY, DoubleVector &dNdZ)
{
    const mtx::size_type elementNodes = 4;
    // значения производных функций формы в местных координатах
    DoubleVector dNdXi(elementNodes);
    DoubleVector dNdEta(elementNodes);
    DoubleVector dNdMu(elementNodes);
    // Якобиан
    DoubleMatrix Jacobi(3, 0.0);
    DoubleMatrix invJacobian(3);
    double jacobian;
    N(0) = 1.0 - xi - eta - mu;
    N(1) = xi;
    N(2) = eta;
    N(3) = mu;


    dNdXi(0) = -1.0;
    dNdXi(1) = 1.0;
    dNdXi(2) = 0.0;
    dNdXi(3) = 0.0;


    dNdEta(0) = -1.0;
    dNdEta(1) = 0.0;
    dNdEta(2) = 1.0;
    dNdEta(3) = 0.0;

    dNdMu(0) = -1.0;
    dNdMu(1) = 0.0;
    dNdMu(2) = 0.0;
    dNdMu(3) = 1.0;

    Jacobi(0, 0) = 0.0; Jacobi(0, 1) = 0.0; Jacobi(0, 2) = 0.0;
    Jacobi(1, 0) = 0.0; Jacobi(1, 1) = 0.0; Jacobi(1, 2) = 0.0;
    Jacobi(2, 0) = 0.0; Jacobi(2, 1) = 0.0; Jacobi(2, 2) = 0.0;

    for (mtx::size_type i = 0; i < elementNodes; i++)
    {
        Jacobi(0, 0) += dNdXi(i) * x[i]; Jacobi(0, 1) += dNdXi(i) * y[i]; Jacobi(0, 2) += dNdXi(i) * z[i];
        Jacobi(1, 0) += dNdEta(i) * x[i]; Jacobi(1, 1) += dNdEta(i) * y[i]; Jacobi(1, 2) += dNdEta(i) * z[i];
        Jacobi(2, 0) += dNdMu(i) * x[i]; Jacobi(2, 1) += dNdMu(i) * y[i]; Jacobi(2, 2) += dNdMu(i) * z[i];
    }
    invJacobian = Jacobi.inverted3x3();
    jacobian = Jacobi.det3x3();

    for (mtx::size_type i = 0; i < elementNodes; i++)
    {
        dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
        dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
        dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
    }
    return jacobian;
}
