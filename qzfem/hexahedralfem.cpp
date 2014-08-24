#include "hexahedralfem.h"

#include "floatingvector.h"
#include "invertmatrix.hpp"
#include "determinant.hpp"
#include "conjugategradient.hpp"
#include <boost/progress.hpp>
#include <iostream>

HexahedralFEM::HexahedralFEM(HexahedralMesh3D *mesh, const MechanicalParameters3D &parameters, FEMCondition3DPointer boundaryForce, std::vector<FEMCondition3DPointer> boundaryConditions)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger systemDimension = nodesCount * freedom;

    FloatingMatrix D(6, 6);

    GlobalMatrix globalMatrix (systemDimension, systemDimension);
    FloatingVector force(systemDimension);

    displacement.resize(systemDimension);

//    std::cout << "E = " << E << " nu = " << nu << std::endl;

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
    compressed_matrix<Floating> CGM(systemDimension, systemDimension);
    CGM.assign(globalMatrix);
    globalMatrix.clear();
    std::cout << "Решение СЛАУ..." << std::endl;
    int res = conjugateGradient(CGM, displacement,force,row_major_tag());
    std::cout << "Выполнено итераций в метода сопряженных градиентов: " << res << std::endl;

    printDisplacementExtremum(nodesCount);
}

HexahedralFEM::HexahedralFEM(HexahedralMesh3D *mesh, const MechanicalParameters3D &parameters, std::vector<FEMCondition3DPointer> boundaryForces, std::vector<FEMCondition3DPointer> boundaryConditions)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger systemDimension = nodesCount * freedom;

    FloatingMatrix D(6, 6);

    GlobalMatrix globalMatrix (systemDimension, systemDimension);
    FloatingVector force(systemDimension);

    displacement.resize(systemDimension);

//    std::cout << "E = " << E << " nu = " << nu << std::endl;

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
    compressed_matrix<Floating> CGM(systemDimension, systemDimension);
    CGM.assign(globalMatrix);
    globalMatrix.clear();
    std::cout << "Решение СЛАУ..." << std::endl;
    int res = conjugateGradient(CGM, displacement,force,row_major_tag());
    std::cout << "Выполнено итераций в метода сопряженных градиентов: " << res << std::endl;

    printDisplacementExtremum(nodesCount);
}

void HexahedralFEM::setNodeDisplacement(HexahedralMesh3D *mesh, const UInteger &direction)
{
    const UInteger nodesCount = mesh->nodesCount();
    mesh->clearNodeValues();
    for (UInteger i = 0; i < nodesCount; i++)
    {
        mesh->pushNodeValue(displacement[i + direction * nodesCount]);
    }
}

Point3D HexahedralFEM::getDisplacemementVector(const UInteger &i, const UInteger &nodesCount)
{
    Point3D point(displacement[i], displacement[i + nodesCount], displacement[i + 2L * nodesCount]);
    return point;
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

//    // Инициализация матрицы D
//    for (int i = 0; i < 6; i++)
//        for (int j = 0; j < 6; j++)
//            D(i, j) = 0.0;

//    D(0, 0) = 1.0; D(0, 1) = nu / (1.0 - nu); D(0, 2) = nu / (1.0 - nu);
//    D(1, 0) = nu / (1.0 - nu); D(1, 1) = 1.0; D(1, 2) = nu / (1.0 - nu);
//    D(2, 0) = nu / (1.0 - nu); D(2, 1) = nu / (1.0 - nu); D(2, 2) = 1.0;
//    D(3, 3) = (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
//    D(4, 4) = (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
//    D(5, 5) = (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
//    D = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) * D;
}

void HexahedralFEM::assebly(HexahedralMesh3D *mesh, const FloatingMatrix &D, GlobalMatrix &globalMatrix)
{
    const int freedom = 3;
    const UInteger nodesCount = mesh->nodesCount();
    const UInteger elementsCount = mesh->elementsCount();
    const int elementNodes = 8;
    const int gaussCount = 2; // количество точек в квадратурах для интегрирования
    Floating gaussPoint[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}; // координаты точек квадратуры
    Floating gaussWeight[] = {1.0, 1.0}; // веса точек квадратуры

    std::cout << "Построение глобальной матрицы..." << std::endl;

    boost::progress_display progressBar(elementsCount);
    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;

        Floating x[elementNodes];
        Floating y[elementNodes];
        Floating z[elementNodes];
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
            Floating xi = gaussPoint[iXi];
            Floating xiWeight = gaussWeight[iXi];
            for (int iEta = 0; iEta < gaussCount; iEta++)
            {
                Floating eta = gaussPoint[iEta];
                Floating etaWeight = gaussWeight[iEta];
                for (int iMu = 0; iMu < gaussCount; iMu++)
                {
                    Floating mu = gaussPoint[iMu];
                    Floating muWeight = gaussWeight[iMu];
                    // значения функций формы
                    FloatingVector N(elementNodes);
                    // значения производных функций формы в местных координатах
                    FloatingVector dNdXi(elementNodes);
                    FloatingVector dNdEta(elementNodes);
                    FloatingVector dNdMu(elementNodes);
                    // Якобиан
                    FloatingMatrix Jacobian(freedom, freedom);
                    FloatingMatrix invJacobian(freedom, freedom);
                    Floating detJacobian;
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

void HexahedralFEM::processForce(HexahedralMesh3D *mesh, FEMCondition3DPointer boundaryForce, FloatingVector &force)
{
    boost::progress_display progressBar(mesh->elementsCount());
    const UInteger nodesCount = mesh->nodesCount();
    Floating area = 0.0;
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
                    Floating faceArea = mesh->faceArea(face);
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

void HexahedralFEM::processBoundaryConditions(HexahedralMesh3D *mesh, std::vector<FEMCondition3DPointer> boundaryConditions, GlobalMatrix &globalMatrix, FloatingVector &force)
{
    const UInteger nodesCount = mesh->nodesCount();
    boost::progress_display progressBar(nodesCount);
    const UInteger freedom = 3;
    const UInteger systemDimension = nodesCount * freedom;

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
}

void HexahedralFEM::printDisplacementExtremum(const UInteger &nodesCount)
{
    std::cout << "Обработка вектора перемещений..." << std::endl;
    boost::progress_display progressBar(nodesCount);
    Floating maxU = displacement(0);
    Floating maxV = displacement(nodesCount);
    Floating maxW = displacement(2UL * nodesCount);
    Floating minU = displacement(0);
    Floating minV = displacement(nodesCount);
    Floating minW = displacement(2UL * nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        ++progressBar;
        // max
        if (maxU < displacement(i))
        {
            maxU = displacement(i);
        }
        if (maxV < displacement(i + nodesCount))
        {
            maxV = displacement(i + nodesCount);
        }
        if (maxW < displacement(i + 2UL * nodesCount))
        {
            maxW = displacement(i + 2UL * nodesCount);
        }
        // min
        if (minU > displacement(i))
        {
            minU = displacement(i);
        }
        if (minV > displacement(i + nodesCount))
        {
            minV = displacement(i + nodesCount);
        }
        if (minW > displacement(i + 2UL * nodesCount))
        {
            minW = displacement(i + 2UL * nodesCount);
        }
    } // for i
    std::cout << "Первое (x) направление:\t" << minU << " <= U <= " << maxU << std::endl;
    std::cout << "Второе (y) направление:\t" << minV << " <= V <= " << maxV << std::endl;
    std::cout << "Третье (z) направление:\t" << minW << " <= W <= " << maxW << std::endl;
}
