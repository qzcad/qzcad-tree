#include "hexahedralfem.h"

#include "floatingvector.h"
#include "invertmatrix.hpp"
#include "determinant.hpp"
#include "conjugategradient.hpp"
#include <boost/progress.hpp>
#include <iostream>
#include <float.h>

HexahedralFEM::HexahedralFEM(HexahedralMesh3D *mesh, const MechanicalParameters3D &parameters, FEMCondition3DPointer boundaryForce, const std::vector<FEMCondition3DPointer> &boundaryConditions)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger systemDimension = nodesCount * freedom;

    FloatingMatrix D(6, 6);

    GlobalMatrix globalMatrix (systemDimension, systemDimension);
    FloatingVector force(systemDimension);

    FloatingVector displacement(systemDimension);

    buildElasticMatrix(parameters, D);

    std::cout << "D:" << std::endl;
    for (int i = 0; i < 6; i++)
    {
        std::cout << "[\t";
        for (int j = 0; j < 6; j++)
        {
            std::cout << D(i, j) << '\t';
        }
        std::cout << "\t]" << std::endl;
    }

    assebly(mesh, D, globalMatrix);

    std::cout << "Построение вектора нагрузок... " << std::endl;
    for (UInteger i = 0; i < systemDimension; i++)
        force(i) = 0.0;

    processForce(mesh, boundaryForce, force);

    processBoundaryConditions(mesh, boundaryConditions, globalMatrix, force);

    std::cout << "Ассоциация и частичная очистка памяти..." << std::endl;
    compressed_matrix<double> CGM(systemDimension, systemDimension);
    CGM.assign(globalMatrix);
    globalMatrix.clear();
    std::cout << "Решение СЛАУ..." << std::endl;
    int res = conjugateGradient(CGM, displacement, force, row_major_tag());
    std::cout << "Выполнено итераций в метода сопряженных градиентов: " << res << std::endl;

    displacementToUVW(displacement, nodesCount);

    printDisplacementExtremum();

    recoverStress(mesh, D);
}

HexahedralFEM::HexahedralFEM(HexahedralMesh3D *mesh, const MechanicalParameters3D &parameters, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger systemDimension = nodesCount * freedom;

    FloatingMatrix D(6, 6);

    GlobalMatrix globalMatrix (systemDimension, systemDimension);
    FloatingVector force(systemDimension);

    FloatingVector displacement(systemDimension);

    buildElasticMatrix(parameters, D);

    std::cout << "D:" << std::endl;
    for (int i = 0; i < 6; i++)
    {
        std::cout << "[\t";
        for (int j = 0; j < 6; j++)
        {
            std::cout << D(i, j) << '\t';
        }
        std::cout << "\t]" << std::endl;
    }

    assebly(mesh, D, globalMatrix);

    std::cout << "Построение вектора нагрузок... " << std::endl;
    for (UInteger i = 0; i < systemDimension; i++)
        force(i) = 0.0;

    for(UInteger bf = 0; bf < boundaryForces.size(); bf++)
    {
        FEMCondition3DPointer forcePointer = boundaryForces[bf];
        processForce(mesh, forcePointer, force);
    } // for bf

    processBoundaryConditions(mesh, boundaryConditions, globalMatrix, force);

    std::cout << "Ассоциация и частичная очистка памяти..." << std::endl;
    compressed_matrix<double> CGM(systemDimension, systemDimension);
    CGM.assign(globalMatrix);
    globalMatrix.clear();
    std::cout << "Решение СЛАУ..." << std::endl;
    int res = conjugateGradient(CGM, displacement, force, row_major_tag());
    std::cout << "Выполнено итераций в метода сопряженных градиентов: " << res << std::endl;

    displacementToUVW(displacement, nodesCount);

    printDisplacementExtremum();

    recoverStress(mesh, D);
}

HexahedralFEM::HexahedralFEM(HexahedralMesh3D *mesh, const std::vector<MechanicalParameters3D> &parameters, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger systemDimension = nodesCount * freedom;

    UInteger layersCount = parameters.size(); // количество слоев

    FloatingMatrix D[layersCount];
    for (UInteger i = 0; i< layersCount; i++) D[i].resize(6, 6);

    GlobalMatrix globalMatrix (systemDimension, systemDimension);
    FloatingVector force(systemDimension);

    FloatingVector displacement(systemDimension);


    buildElasticMatrix(parameters, D);

    std::cout << "D:" << std::endl;
    for (UInteger p = 0; p < layersCount; p++)
        for (int i = 0; i < 6; i++)
        {
            std::cout << "[\t";
            for (int j = 0; j < 6; j++)
            {
                std::cout << D[p](i, j) << '\t';
            }
            std::cout << "\t]" << std::endl;
        }

    assebly(mesh, D, globalMatrix);

    std::cout << "Построение вектора нагрузок... " << std::endl;
    for (UInteger i = 0; i < systemDimension; i++)
        force(i) = 0.0;

    for(UInteger bf = 0; bf < boundaryForces.size(); bf++)
    {
        FEMCondition3DPointer forcePointer = boundaryForces[bf];
        processForce(mesh, forcePointer, force);
    } // for bf

    processBoundaryConditions(mesh, boundaryConditions, globalMatrix, force);

    std::cout << "Ассоциация и частичная очистка памяти..." << std::endl;
    compressed_matrix<double> CGM(systemDimension, systemDimension);
    CGM.assign(globalMatrix);
    globalMatrix.clear();
    std::cout << "Решение СЛАУ..." << std::endl;
    int res = conjugateGradient(CGM, displacement,force,row_major_tag());
    std::cout << "Выполнено итераций в метода сопряженных градиентов: " << res << std::endl;

    displacementToUVW(displacement, nodesCount);

    printDisplacementExtremum();

    recoverStress(mesh, D);
}

std::vector<double> HexahedralFEM::u() const
{
    return u_;
}

std::vector<double> HexahedralFEM::v() const
{
    return v_;
}

std::vector<double> HexahedralFEM::w() const
{
    return w_;
}

std::vector<double> HexahedralFEM::sigmaX() const
{
    return sigmaX_;
}

std::vector<double> HexahedralFEM::sigmaY() const
{
    return sigmaY_;
}

std::vector<double> HexahedralFEM::sigmaZ() const
{
    return sigmaZ_;
}

std::vector<double> HexahedralFEM::tauXY() const
{
    return tauXY_;
}

std::vector<double> HexahedralFEM::tauYZ() const
{
    return tauYZ_;
}

std::vector<double> HexahedralFEM::tauZX() const
{
    return tauZX_;
}

void HexahedralFEM::buildElasticMatrix(const MechanicalParameters3D &params, FloatingMatrix &D)
{
    // матрица податливости
    FloatingMatrix S(6, 6);

    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            S(i, j) = 0.0;

    S(0, 0) = 1.0 / params.E1();    S(0, 1) = - params.nu21() / params.E2();    S(0, 2) = - params.nu31() / params.E3();
    S(1, 0) = - params.nu12() / params.E1();    S(1, 1) = 1.0 / params.E2();    S(1, 2) = - params.nu32() / params.E3();
    S(2, 0) = - params.nu13() / params.E1();    S(2, 1) = - params.nu23() / params.E2();    S(2, 2) = 1.0 / params.E3();
    S(3, 3) = 1.0 / params.G23();
    S(4, 4) = 1.0 / params.G13();
    S(5, 5) = 1.0 / params.G12();

    invertMatrix(S, D);
}

void HexahedralFEM::buildElasticMatrix(const std::vector<MechanicalParameters3D> &params, FloatingMatrix D[])
{
    for (UInteger p = 0; p < params.size(); p++)
    {
        buildElasticMatrix(params[p], D[p]);
    }
}

void HexahedralFEM::assebly(HexahedralMesh3D *mesh, const FloatingMatrix &D, GlobalMatrix &globalMatrix)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger elementsCount = mesh->elementsCount();
    const int elementNodes = 8;
    const int gaussCount = 2; // количество точек в квадратурах для интегрирования
    double gaussPoint[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}; // координаты точек квадратуры
    double gaussWeight[] = {1.0, 1.0}; // веса точек квадратуры

    std::cout << "Построение глобальной матрицы..." << std::endl;

    boost::progress_display progressBar(elementsCount);
    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;

        double x[elementNodes];
        double y[elementNodes];
        double z[elementNodes];
        FloatingMatrix localMatrix(elementNodes * freedom, elementNodes * freedom);
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            for (int j = 0; j < elementNodes * freedom; j++)
            {
                localMatrix(i, j) = 0.0;
            }
        }
        // извлечение координат узлов
        ElementPointer element = mesh->element(elementNumber);
        for (int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }
        // интегрирование квадратурами
        for (int iXi = 0; iXi < gaussCount; iXi++)
        {
            double xi = gaussPoint[iXi];
            double xiWeight = gaussWeight[iXi];
            for (int iEta = 0; iEta < gaussCount; iEta++)
            {
                double eta = gaussPoint[iEta];
                double etaWeight = gaussWeight[iEta];
                for (int iMu = 0; iMu < gaussCount; iMu++)
                {
                    double mu = gaussPoint[iMu];
                    double muWeight = gaussWeight[iMu];
                    // значения функций формы
                    FloatingVector N(elementNodes);
                    // значения производных функций формы в местных координатах
                    FloatingVector dNdXi(elementNodes);
                    FloatingVector dNdEta(elementNodes);
                    FloatingVector dNdMu(elementNodes);
                    // Якобиан
                    FloatingMatrix Jacobian(freedom, freedom);
                    FloatingMatrix invJacobian(freedom, freedom);
                    double detJacobian;
                    // значения производных функций формы в глобальных координатах
                    FloatingVector dNdX(elementNodes);
                    FloatingVector dNdY(elementNodes);
                    FloatingVector dNdZ(elementNodes);
                    // матрицы вариационной постановки
                    FloatingMatrix B(6, 24);
                    FloatingMatrix BT(24, 6);
                    FloatingMatrix BTD(24, 6);

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

                    Jacobian(0, 0) = 0.0; Jacobian(0, 1) = 0.0; Jacobian(0, 2) = 0.0;
                    Jacobian(1, 0) = 0.0; Jacobian(1, 1) = 0.0; Jacobian(1, 2) = 0.0;
                    Jacobian(2, 0) = 0.0; Jacobian(2, 1) = 0.0; Jacobian(2, 2) = 0.0;

                    for (int i = 0; i < elementNodes; i++)
                    {
                        Jacobian(0, 0) += dNdXi(i) * x[i]; Jacobian(0, 1) += dNdXi(i) * y[i]; Jacobian(0, 2) += dNdXi(i) * z[i];
                        Jacobian(1, 0) += dNdEta(i) * x[i]; Jacobian(1, 1) += dNdEta(i) * y[i]; Jacobian(1, 2) += dNdEta(i) * z[i];
                        Jacobian(2, 0) += dNdMu(i) * x[i]; Jacobian(2, 1) += dNdMu(i) * y[i]; Jacobian(2, 2) += dNdMu(i) * z[i];
                    }
                    invertMatrix(Jacobian, invJacobian);
                    detJacobian = determinant(Jacobian);

                    for (int i = 0; i < elementNodes; i++)
                    {
                        dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
                        dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
                        dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
                    }

                    for (int i = 0; i < 6; i++)
                    {
                        for (int j = 0; j < 24; j++)
                        {
                            B(i, j) = 0.0;
                        }
                    }
                    for (int i = 0; i < 8; i++)
                    {
                        B(0, i) = dNdX(i);
                        B(1, i + 8) = dNdY(i);
                        B(2, i + 16) = dNdZ(i);
                        B(3, i) = dNdY(i); B(3, i + 8) = dNdX(i);
                        B(4, i + 8) = dNdZ(i); B(4, i + 16) = dNdY(i);
                        B(5, i) = dNdZ(i); B(5, i + 16) = dNdX(i);
                    }
                    BT = trans(B);

                    BTD = prod(BT, D);

                    localMatrix += prod(BTD, B) * detJacobian * xiWeight * etaWeight * muWeight;

                } // for iMu
            } // for iEta
        } // for iXi
        // Ансамблирование
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }
            else
            {
                index_i = element->vertexNode(i - 2L * elementNodes) + 2L * nodesCount;
            }

            for (int j = 0; j < elementNodes * freedom; j++)
            {
                if (j < elementNodes)
                {
                    index_j = element->vertexNode(j);
                }
                else if (j < 2 * elementNodes)
                {
                    index_j = element->vertexNode(j - elementNodes) + nodesCount;
                }
                else
                {
                    index_j = element->vertexNode(j - 2L * elementNodes) + 2L * nodesCount;
                }
                globalMatrix(index_i, index_j) += localMatrix(i, j);
            } // for j
        } // for i
    } // for elementNumber
}

void HexahedralFEM::assebly(HexahedralMesh3D *mesh, FloatingMatrix D[], GlobalMatrix &globalMatrix)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger elementsCount = mesh->elementsCount();
    const int elementNodes = 8;
    const int gaussCount = 2; // количество точек в квадратурах для интегрирования
    double gaussPoint[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}; // координаты точек квадратуры
    double gaussWeight[] = {1.0, 1.0}; // веса точек квадратуры

    std::cout << "Построение глобальной матрицы..." << std::endl;

    boost::progress_display progressBar(elementsCount);
    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;

        double x[elementNodes];
        double y[elementNodes];
        double z[elementNodes];
        FloatingMatrix localMatrix(elementNodes * freedom, elementNodes * freedom);
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            for (int j = 0; j < elementNodes * freedom; j++)
            {
                localMatrix(i, j) = 0.0;
            }
        }
        // извлечение координат узлов
        ElementPointer element = mesh->element(elementNumber);
        for (int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }
        // интегрирование квадратурами
        for (int iXi = 0; iXi < gaussCount; iXi++)
        {
            double xi = gaussPoint[iXi];
            double xiWeight = gaussWeight[iXi];
            for (int iEta = 0; iEta < gaussCount; iEta++)
            {
                double eta = gaussPoint[iEta];
                double etaWeight = gaussWeight[iEta];
                for (int iMu = 0; iMu < gaussCount; iMu++)
                {
                    double mu = gaussPoint[iMu];
                    double muWeight = gaussWeight[iMu];
                    // значения функций формы
                    FloatingVector N(elementNodes);
                    // значения производных функций формы в местных координатах
                    FloatingVector dNdXi(elementNodes);
                    FloatingVector dNdEta(elementNodes);
                    FloatingVector dNdMu(elementNodes);
                    // Якобиан
                    FloatingMatrix Jacobian(freedom, freedom);
                    FloatingMatrix invJacobian(freedom, freedom);
                    double detJacobian;
                    // значения производных функций формы в глобальных координатах
                    FloatingVector dNdX(elementNodes);
                    FloatingVector dNdY(elementNodes);
                    FloatingVector dNdZ(elementNodes);
                    // матрицы вариационной постановки
                    FloatingMatrix B(6, 24);
                    FloatingMatrix BT(24, 6);
                    FloatingMatrix BTD(24, 6);

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

                    Jacobian(0, 0) = 0.0; Jacobian(0, 1) = 0.0; Jacobian(0, 2) = 0.0;
                    Jacobian(1, 0) = 0.0; Jacobian(1, 1) = 0.0; Jacobian(1, 2) = 0.0;
                    Jacobian(2, 0) = 0.0; Jacobian(2, 1) = 0.0; Jacobian(2, 2) = 0.0;

                    for (int i = 0; i < elementNodes; i++)
                    {
                        Jacobian(0, 0) += dNdXi(i) * x[i]; Jacobian(0, 1) += dNdXi(i) * y[i]; Jacobian(0, 2) += dNdXi(i) * z[i];
                        Jacobian(1, 0) += dNdEta(i) * x[i]; Jacobian(1, 1) += dNdEta(i) * y[i]; Jacobian(1, 2) += dNdEta(i) * z[i];
                        Jacobian(2, 0) += dNdMu(i) * x[i]; Jacobian(2, 1) += dNdMu(i) * y[i]; Jacobian(2, 2) += dNdMu(i) * z[i];
                    }
                    invertMatrix(Jacobian, invJacobian);
                    detJacobian = determinant(Jacobian);

                    for (int i = 0; i < elementNodes; i++)
                    {
                        dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
                        dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
                        dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
                    }

                    for (int i = 0; i < 6; i++)
                    {
                        for (int j = 0; j < 24; j++)
                        {
                            B(i, j) = 0.0;
                        }
                    }
                    for (int i = 0; i < 8; i++)
                    {
                        B(0, i) = dNdX(i);
                        B(1, i + 8) = dNdY(i);
                        B(2, i + 16) = dNdZ(i);
                        B(3, i) = dNdY(i); B(3, i + 8) = dNdX(i);
                        B(4, i + 8) = dNdZ(i); B(4, i + 16) = dNdY(i);
                        B(5, i) = dNdZ(i); B(5, i + 16) = dNdX(i);
                    }
                    BT = trans(B);

                    int dNum = (int)mesh->elementValue(elementNumber);
                    BTD = prod(BT, D[dNum]);

                    localMatrix += prod(BTD, B) * detJacobian * xiWeight * etaWeight * muWeight;

                } // for iMu
            } // for iEta
        } // for iXi
        // Ансамблирование
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }
            else
            {
                index_i = element->vertexNode(i - 2L * elementNodes) + 2L * nodesCount;
            }

            for (int j = 0; j < elementNodes * freedom; j++)
            {
                if (j < elementNodes)
                {
                    index_j = element->vertexNode(j);
                }
                else if (j < 2 * elementNodes)
                {
                    index_j = element->vertexNode(j - elementNodes) + nodesCount;
                }
                else
                {
                    index_j = element->vertexNode(j - 2L * elementNodes) + 2L * nodesCount;
                }
                globalMatrix(index_i, index_j) += localMatrix(i, j);
            } // for j
        } // for i
    } // for elementNumber
}

void HexahedralFEM::processForce(HexahedralMesh3D *mesh, FEMCondition3DPointer boundaryForce, FloatingVector &force)
{
    boost::progress_display progressBar(mesh->elementsCount());
    const UInteger nodesCount = mesh->nodesCount();
    double area = 0.0;
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        ++progressBar;
        if (mesh->isBorderElement(i))
        {
            ElementPointer element = mesh->element(i);
            for (int j = 0; j < element->facesCount(); j++)
            {
                UIntegerVector face = element->face(j);
                bool isBorderFace = true;
                for (int k = 0; k < 4; k++)
                {
                    if (mesh->nodeType(face[k]) == INNER || !boundaryForce->isApplied(mesh->node(face[k])))
                    {
                        isBorderFace = false;
                    }
                } // for k
                if (isBorderFace)
                {
                    double faceArea = mesh->faceArea(face);
                    area += faceArea;
                    for (int k = 0; k < 4; k++)
                    {
                        force(face[k]) = force(face[k]) + boundaryForce->u() * faceArea / 4.0;
                        force(face[k] + nodesCount) = force(face[k] + nodesCount) + boundaryForce->v() * faceArea / 4.0;
                        force(face[k] + 2L * nodesCount) = force(face[k] + 2L * nodesCount) + boundaryForce->w() * faceArea / 4.0;
                    } // for k
                }// if
            } // for j
        } // if
    } // for i
    std::cout << "Нагрузка: " << boundaryForce->u() << "; " << boundaryForce->v() << "; " << boundaryForce->w() << ". Площадь нагруженной поверхноти: " << area << std::endl;
}

void HexahedralFEM::processBoundaryConditions(HexahedralMesh3D *mesh, const std::vector<FEMCondition3DPointer> &boundaryConditions, GlobalMatrix &globalMatrix, FloatingVector &force)
{
    const UInteger nodesCount = mesh->nodesCount();
    boost::progress_display progressBar(nodesCount);
    const UInteger freedom = 3;
    const UInteger systemDimension = nodesCount * freedom;
    UInteger boundaryNodes = 0;

    std::cout << "Учет граничных условий..." << std::endl;

    for (UInteger i = 0; i < nodesCount; i++)
    {
        ++progressBar;
        PointPointer point = mesh->node(i);
        for (UInteger b = 0; b < boundaryConditions.size(); b++)
        {
            FEMCondition3DPointer condition = boundaryConditions[b];
            if (condition->isApplied(point))
            {
                boundaryNodes++;
                if (condition->isU())
                {
                    for (UInteger j = 0; j < systemDimension; j++)
                    {
                        if (i != j && globalMatrix(i, j) != 0.0)
                        { // см. Зенкевич, стр. 485
                            force(j) = force(j) - globalMatrix(i, j) * condition->u();
                            globalMatrix(i, j) = 0.0; // обнуление строки/столбца
                            globalMatrix(j, i) = 0.0;
                        }
                    }
                    force(i) = condition->u();
                    globalMatrix(i, i) = 1.0;
                } // if
                if (condition->isV())
                {
                    UInteger rowNumber = i + nodesCount;
                    for (UInteger j = 0; j < systemDimension; j++)
                    {
                        if (rowNumber != j && globalMatrix(rowNumber, j) != 0.0)
                        { // см. Зенкевич, стр. 485
                            force(j) = force(j) - globalMatrix(rowNumber, j) * condition->v();
                            globalMatrix(rowNumber, j) = 0.0; // обнуление строки/столбца
                            globalMatrix(j, rowNumber) = 0.0;
                        }
                    }
                    force(rowNumber) = condition->v();
                    globalMatrix(rowNumber, rowNumber) = 1.0;
                } // if
                if (condition->isW())
                {
                    UInteger rowNumber = i + 2L * nodesCount;
                    for (UInteger j = 0; j < systemDimension; j++)
                    {
                        if (rowNumber != j && globalMatrix(rowNumber, j) != 0.0)
                        { // см. Зенкевич, стр. 485
                            force(j) = force(j) - globalMatrix(rowNumber, j) * condition->w();
                            globalMatrix(rowNumber, j) = 0.0; // обнуление строки/столбца
                            globalMatrix(j, rowNumber) = 0.0;
                        }
                    }
                    force(rowNumber) = condition->w();
                    globalMatrix(rowNumber, rowNumber) = 1.0;
                } // if
            } // if
        } // for b
    } // for i
    std::cout << "Обработано " << boundaryNodes << " узлов." << std::endl;
}

void HexahedralFEM::printDisplacementExtremum()
{
    std::cout << "Обработка вектора перемещений..." << std::endl;
    boost::progress_display progressBar(u_.size());
    double maxU = u_[0];
    double maxV = v_[0];
    double maxW = w_[0];
    double minU = u_[0];
    double minV = v_[0];
    double minW = w_[0];
    for (UInteger i = 0; i < u_.size(); i++)
    {
        ++progressBar;
        double u = u_[i], v = v_[i], w = w_[i];
        // max
        if (maxU < u)
        {
            maxU = u;
        }
        if (maxV < v)
        {
            maxV = v;
        }
        if (maxW < w)
        {
            maxW = w;
        }
        // min
        if (minU > u)
        {
            minU = u;
        }
        if (minV > v)
        {
            minV = v;
        }
        if (minW > w)
        {
            minW = w;
        }
    } // for i
    std::cout << "Первое (x) направление:\t" << minU << " <= U <= " << maxU << std::endl;
    std::cout << "Второе (y) направление:\t" << minV << " <= V <= " << maxV << std::endl;
    std::cout << "Третье (z) направление:\t" << minW << " <= W <= " << maxW << std::endl;
}

void HexahedralFEM::recoverStress(HexahedralMesh3D *mesh, const FloatingMatrix &D)
{
    const int freedom = 3;
    const UInteger elementsCount = mesh->elementsCount();
    const int elementNodes = 8;

    double minSigmaX = DBL_MAX;
    double minSigmaY = DBL_MAX;
    double minSigmaZ = DBL_MAX;
    double maxSigmaX = -DBL_MAX;
    double maxSigmaY = -DBL_MAX;
    double maxSigmaZ = -DBL_MAX;
    double minTauXY = DBL_MAX;
    double minTauYZ = DBL_MAX;
    double minTauZX = DBL_MAX;
    double maxTauXY = -DBL_MAX;
    double maxTauYZ = -DBL_MAX;
    double maxTauZX = -DBL_MAX;

    std::cout << "Вычисление напряжений..." << std::endl;

    boost::progress_display progressBar(elementsCount);
    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;

        double x[elementNodes];
        double y[elementNodes];
        double z[elementNodes];
        // извлечение координат узлов
        ElementPointer element = mesh->element(elementNumber);
        for (int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }
        const double xi = 0.0;
        const double eta = 0.0;
        const double mu = 0.0;

        // значения функций формы
        FloatingVector N(elementNodes);
        // значения производных функций формы в местных координатах
        FloatingVector dNdXi(elementNodes);
        FloatingVector dNdEta(elementNodes);
        FloatingVector dNdMu(elementNodes);
        // Якобиан
        FloatingMatrix Jacobian(freedom, freedom);
        FloatingMatrix invJacobian(freedom, freedom);
        // значения производных функций формы в глобальных координатах
        FloatingVector dNdX(elementNodes);
        FloatingVector dNdY(elementNodes);
        FloatingVector dNdZ(elementNodes);
        // матрицы вариационной постановки
        FloatingMatrix B(6, 24);
        FloatingMatrix DB(6, 24);
        FloatingMatrix U(24, 1);
        FloatingMatrix S(6, 1);

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

        Jacobian(0, 0) = 0.0; Jacobian(0, 1) = 0.0; Jacobian(0, 2) = 0.0;
        Jacobian(1, 0) = 0.0; Jacobian(1, 1) = 0.0; Jacobian(1, 2) = 0.0;
        Jacobian(2, 0) = 0.0; Jacobian(2, 1) = 0.0; Jacobian(2, 2) = 0.0;

        for (int i = 0; i < elementNodes; i++)
        {
            Jacobian(0, 0) += dNdXi(i) * x[i]; Jacobian(0, 1) += dNdXi(i) * y[i]; Jacobian(0, 2) += dNdXi(i) * z[i];
            Jacobian(1, 0) += dNdEta(i) * x[i]; Jacobian(1, 1) += dNdEta(i) * y[i]; Jacobian(1, 2) += dNdEta(i) * z[i];
            Jacobian(2, 0) += dNdMu(i) * x[i]; Jacobian(2, 1) += dNdMu(i) * y[i]; Jacobian(2, 2) += dNdMu(i) * z[i];
        }
        invertMatrix(Jacobian, invJacobian);

        for (int i = 0; i < elementNodes; i++)
        {
            dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
            dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
            dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
        }

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 24; j++)
            {
                B(i, j) = 0.0;
            }
        }
        for (int i = 0; i < 8; i++)
        {
            B(0, i) = dNdX(i);
            B(1, i + 8) = dNdY(i);
            B(2, i + 16) = dNdZ(i);
            B(3, i) = dNdY(i); B(3, i + 8) = dNdX(i);
            B(4, i + 8) = dNdZ(i); B(4, i + 16) = dNdY(i);
            B(5, i) = dNdZ(i); B(5, i + 16) = dNdX(i);

            U(i, 0) = u_[element->vertexNode(i)];
            U(i + 8, 0) = v_[element->vertexNode(i)];
            U(i + 16, 0) = w_[element->vertexNode(i)];
        }

        DB = prod(D, B);

        S = prod(DB, U);

        sigmaX_.push_back(S(0, 0));
        sigmaY_.push_back(S(1, 0));
        sigmaZ_.push_back(S(2, 0));
        tauXY_.push_back(S(3, 0));
        tauYZ_.push_back(S(4, 0));
        tauZX_.push_back(S(5, 0));

        if (S(0, 0) > maxSigmaX) maxSigmaX = S(0, 0);
        if (S(0, 0) < minSigmaX) minSigmaX = S(0, 0);

        if (S(1, 0) > maxSigmaY) maxSigmaY = S(1, 0);
        if (S(1, 0) < minSigmaY) minSigmaY = S(1, 0);

        if (S(2, 0) > maxSigmaZ) maxSigmaZ = S(2, 0);
        if (S(2, 0) < minSigmaZ) minSigmaZ = S(2, 0);

        if (S(3, 0) > maxTauXY) maxTauXY = S(3, 0);
        if (S(3, 0) < minTauXY) minTauXY = S(3, 0);

        if (S(4, 0) > maxTauYZ) maxTauYZ = S(4, 0);
        if (S(4, 0) < minTauYZ) minTauYZ = S(4, 0);

        if (S(5, 0) > maxTauZX) maxTauZX = S(5, 0);
        if (S(5, 0) < minTauZX) minTauZX = S(5, 0);
    } // for elementNumber
    std::cout << minSigmaX << " <= SigmaX <= " << maxSigmaX << std::endl;
    std::cout << minSigmaY << " <= SigmaY <= " << maxSigmaY << std::endl;
    std::cout << minSigmaZ << " <= SigmaZ <= " << maxSigmaZ << std::endl;
    std::cout << minTauXY << " <= TauXY <= " << maxTauXY << std::endl;
    std::cout << minTauYZ << " <= TauYZ <= " << maxTauYZ << std::endl;
    std::cout << minTauZX << " <= TauZX <= " << maxTauZX << std::endl;
}

void HexahedralFEM::recoverStress(HexahedralMesh3D *mesh, const FloatingMatrix D[])
{
    const int freedom = 3;
    const UInteger elementsCount = mesh->elementsCount();
    const int elementNodes = 8;

    double minSigmaX = DBL_MAX;
    double minSigmaY = DBL_MAX;
    double minSigmaZ = DBL_MAX;
    double maxSigmaX = -DBL_MAX;
    double maxSigmaY = -DBL_MAX;
    double maxSigmaZ = -DBL_MAX;
    double minTauXY = DBL_MAX;
    double minTauYZ = DBL_MAX;
    double minTauZX = DBL_MAX;
    double maxTauXY = -DBL_MAX;
    double maxTauYZ = -DBL_MAX;
    double maxTauZX = -DBL_MAX;

    std::cout << "Вычисление напряжений..." << std::endl;

    boost::progress_display progressBar(elementsCount);
    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;

        double x[elementNodes];
        double y[elementNodes];
        double z[elementNodes];
        // извлечение координат узлов
        ElementPointer element = mesh->element(elementNumber);
        for (int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }
        const double xi = 0.0;
        const double eta = 0.0;
        const double mu = 0.0;

        // значения функций формы
        FloatingVector N(elementNodes);
        // значения производных функций формы в местных координатах
        FloatingVector dNdXi(elementNodes);
        FloatingVector dNdEta(elementNodes);
        FloatingVector dNdMu(elementNodes);
        // Якобиан
        FloatingMatrix Jacobian(freedom, freedom);
        FloatingMatrix invJacobian(freedom, freedom);
        // значения производных функций формы в глобальных координатах
        FloatingVector dNdX(elementNodes);
        FloatingVector dNdY(elementNodes);
        FloatingVector dNdZ(elementNodes);
        // матрицы вариационной постановки
        FloatingMatrix B(6, 24);
        FloatingMatrix DB(6, 24);
        FloatingMatrix U(24, 1);
        FloatingMatrix S(6, 1);

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

        Jacobian(0, 0) = 0.0; Jacobian(0, 1) = 0.0; Jacobian(0, 2) = 0.0;
        Jacobian(1, 0) = 0.0; Jacobian(1, 1) = 0.0; Jacobian(1, 2) = 0.0;
        Jacobian(2, 0) = 0.0; Jacobian(2, 1) = 0.0; Jacobian(2, 2) = 0.0;

        for (int i = 0; i < elementNodes; i++)
        {
            Jacobian(0, 0) += dNdXi(i) * x[i]; Jacobian(0, 1) += dNdXi(i) * y[i]; Jacobian(0, 2) += dNdXi(i) * z[i];
            Jacobian(1, 0) += dNdEta(i) * x[i]; Jacobian(1, 1) += dNdEta(i) * y[i]; Jacobian(1, 2) += dNdEta(i) * z[i];
            Jacobian(2, 0) += dNdMu(i) * x[i]; Jacobian(2, 1) += dNdMu(i) * y[i]; Jacobian(2, 2) += dNdMu(i) * z[i];
        }
        invertMatrix(Jacobian, invJacobian);

        for (int i = 0; i < elementNodes; i++)
        {
            dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
            dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
            dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
        }

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 24; j++)
            {
                B(i, j) = 0.0;
            }
        }
        for (int i = 0; i < 8; i++)
        {
            B(0, i) = dNdX(i);
            B(1, i + 8) = dNdY(i);
            B(2, i + 16) = dNdZ(i);
            B(3, i) = dNdY(i); B(3, i + 8) = dNdX(i);
            B(4, i + 8) = dNdZ(i); B(4, i + 16) = dNdY(i);
            B(5, i) = dNdZ(i); B(5, i + 16) = dNdX(i);

            U(i, 0) = u_[element->vertexNode(i)];
            U(i + 8, 0) = v_[element->vertexNode(i)];
            U(i + 16, 0) = w_[element->vertexNode(i)];
        }

        int dNum = (int)mesh->elementValue(elementNumber);

        DB = prod(D[dNum], B);

        S = prod(DB, U);

        sigmaX_.push_back(S(0, 0));
        sigmaY_.push_back(S(1, 0));
        sigmaZ_.push_back(S(2, 0));
        tauXY_.push_back(S(3, 0));
        tauYZ_.push_back(S(4, 0));
        tauZX_.push_back(S(5, 0));

        if (S(0, 0) > maxSigmaX) maxSigmaX = S(0, 0);
        if (S(0, 0) < minSigmaX) minSigmaX = S(0, 0);

        if (S(1, 0) > maxSigmaY) maxSigmaY = S(1, 0);
        if (S(1, 0) < minSigmaY) minSigmaY = S(1, 0);

        if (S(2, 0) > maxSigmaZ) maxSigmaZ = S(2, 0);
        if (S(2, 0) < minSigmaZ) minSigmaZ = S(2, 0);

        if (S(3, 0) > maxTauXY) maxTauXY = S(3, 0);
        if (S(3, 0) < minTauXY) minTauXY = S(3, 0);

        if (S(4, 0) > maxTauYZ) maxTauYZ = S(4, 0);
        if (S(4, 0) < minTauYZ) minTauYZ = S(4, 0);

        if (S(5, 0) > maxTauZX) maxTauZX = S(5, 0);
        if (S(5, 0) < minTauZX) minTauZX = S(5, 0);
    } // for elementNumber
    std::cout << minSigmaX << " <= SigmaX <= " << maxSigmaX << std::endl;
    std::cout << minSigmaY << " <= SigmaY <= " << maxSigmaY << std::endl;
    std::cout << minSigmaZ << " <= SigmaZ <= " << maxSigmaZ << std::endl;
    std::cout << minTauXY << " <= TauXY <= " << maxTauXY << std::endl;
    std::cout << minTauYZ << " <= TauYZ <= " << maxTauYZ << std::endl;
    std::cout << minTauZX << " <= TauZX <= " << maxTauZX << std::endl;
}

void HexahedralFEM::displacementToUVW(const FloatingVector &displacement, const UInteger &nodesCount)
{
    for (UInteger i = 0; i < nodesCount; i++)
    {
        u_.push_back(displacement(i));
        v_.push_back(displacement(i + nodesCount));
        w_.push_back(displacement(i + 2UL * nodesCount));
    }
}
