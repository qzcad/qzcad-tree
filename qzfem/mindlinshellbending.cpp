#include "mindlinshellbending.h"

#include <iostream>
#include "consoleprogress.h"
#include "rowdoublematrix.h"
#include <math.h>

MindlinShellBending::MindlinShellBending(Mesh3D *mesh, double thickness, const ElasticMatrix &elasticMatrix, std::list<FemCondition *> conditions) :
    Fem2D(mesh, 6)
{
    const double kappa = 5.0 / 6.0;

    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 5; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
    {
        gaussPoints = line_count * line_count;
        gxi.resize(gaussPoints);
        geta.resize(gaussPoints);
        gweight.resize(gaussPoints);
        elementNodes = 4;
        for (int i = 0; i < line_count; i++)
        {
            for (int j = 0; j < line_count; j++)
            {
                gxi[i * line_count + j] = line_points[i];
                geta[i * line_count + j] = line_points[j];
                gweight[i * line_count + j] = line_weights[i] * line_weights[j];
            }
        }
    }

    DoubleMatrix D = elasticMatrix.D(); // матрица упругости
    DoubleMatrix Dc(2, 0.0);
    Dc(0, 0) = D(2, 2); Dc(1, 1) = D(2, 2);

    std::cout << "D:" << std::endl;
    D.print();
    std::cout << "Dc:" << std::endl;
    Dc.print();
    std::cout << "Thickness: " << thickness << std::endl;

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    MappedDoubleMatrix global (dimension_); // глобальная матрица жесткости
    DoubleVector force(dimension_, 0.0); // вектор сил

    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        ElementPointer element = mesh->element(elNum);
        Point3D A = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(0))));
        Point3D B = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(1))));
        Point3D C = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(2))));
        DoubleMatrix lambda = cosinuses(A, B, C);
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов

        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
        }

        for (int ig = 0; ig < gaussPoints; ig++)
        {
            double xi = gxi(ig);
            double eta = geta(ig);
            double w = gweight(ig);
            // значения функций формы
            DoubleVector N(elementNodes);
            // значения производных функций формы
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            // якобиан
            double jacobian = 1.0;
            if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
            {
                jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
            {
                jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bm(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);

            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bm(0, i * freedom_) = dNdX(i);
                                                    Bm(1, i * freedom_ + 1) = dNdY(i);
                Bm(2, i * freedom_) = dNdY(i);      Bm(2, i * freedom_ + 1) = dNdX(i);

                Bf(0, i * freedom_ + 3) = dNdX(i);
                                                    Bf(1, i * freedom_ + 4) = dNdY(i);
                Bf(2, i * freedom_ + 3) = dNdY(i);  Bf(2, i * freedom_ + 4) = dNdX(i);

                Bc(0, i * freedom_ + 2) = dNdX(i);  Bc(0, i * freedom_ + 3) = N(i);
                Bc(1, i * freedom_ + 2) = dNdY(i);  Bc(1, i * freedom_ + 4) = N(i);
            }

            local += jacobian * w * thickness * (Bm.transpose() * D * Bm);
            local += jacobian * w * thickness*thickness*thickness / 12.0 * (Bf.transpose() * D * Bf);
            local += jacobian * w * kappa * thickness * (Bc.transpose() * Dc * Bc);
        } // ig

        DoubleMatrix surf = T.transpose() * local * T;
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (UInteger i = 0; i < elementNodes * freedom_; i++)
        {
            index_i = element->vertexNode(i / freedom_) + (i % freedom_) * nodesCount;

            for (UInteger j = i; j < elementNodes * freedom_; j++)
            {
                index_j = element->vertexNode(j / freedom_) + (j % freedom_) * nodesCount;
                global(index_i, index_j) += surf(i, j);
                if (index_i != index_j) global(index_j, index_i) = global(index_i, index_j);
            } // for j
        } // for i
    } //for elNum

    // Учет сил
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::NODAL_FORCE)
        {
            // узловые нагрузки
            std::cout << "Nodal Forces...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    double f = (*condition)->value(point);
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(i) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(i + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(i + nodesCount + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                        force(i + nodesCount + nodesCount + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                        force(i + nodesCount + nodesCount + nodesCount + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                        force(i + nodesCount + nodesCount + nodesCount + nodesCount + nodesCount) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Edge Distributed Forces...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;

                if (mesh->isBorderElement(elNum))
                {
                    ElementPointer element = mesh->element(elNum);
                    for (unsigned int i = 0; i < elementNodes; i++)
                    {
                        if ((mesh->nodeType(element->vertexNode(i)) == BORDER || mesh->nodeType(element->vertexNode(i)) == CHARACTER) &&
                                (mesh->nodeType(element->vertexNode(i + 1)) == BORDER || mesh->nodeType(element->vertexNode(i + 1)) == CHARACTER))
                        {
                            PointPointer point0 = mesh->node(element->vertexNode(i));
                            PointPointer point1 = mesh->node(element->vertexNode(i + 1));
                            FemCondition::FemDirection dir = (*condition)->direction();
                            if ((*condition)->isApplied(point0) && (*condition)->isApplied(point1))
                            {
                                Point3D p0(point0->x(), point0->y(), point0->z());
                                Point3D p1(point1->x(), point1->y(), point1->z());
                                double l = p0.distanceTo(p1);
                                double jacobian = l / 2.0;
                                double f0 = 0.0;
                                double f1 = 0.0;
                                for (int ixi = 0; ixi < line_count; ixi++)
                                {

                                    double xi = line_points(ixi);
                                    double w = line_weights(ixi);
                                    double N0 = (1.0 - xi) / 2.0;
                                    double N1 = (1.0 + xi) / 2.0;
                                    Point3D p = Point3D(p0.x() * N0 + p1.x() * N1,
                                                        p0.y() * N0 + p1.y() * N1,
                                                        p1.z() * N0 + p1.z() * N1);
                                    f0 += N0 * jacobian * w * (*condition)->value(&p);
                                    f1 += N1 * jacobian * w * (*condition)->value(&p);
                                } // for ixi
                                if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                                {
                                    force(element->vertexNode(i)) += f0;
                                    force(element->vertexNode(i + 1)) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                                {
                                    force(element->vertexNode(i) + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                                {
                                    force(element->vertexNode(i) + 2UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 2UL * nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                                {
                                    force(element->vertexNode(i) + 3UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 3UL * nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                                {
                                    force(element->vertexNode(i) + 4UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 4UL * nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                                {
                                    force(element->vertexNode(i) + 5UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 5UL * nodesCount) += f1;
                                }
                            }
                        } // if
                    } // for i
                } // if
            } // for elNum
        }
        else if ((*condition)->type() == FemCondition::VOLUME_FORCE)
        {
            // объемные силы
            std::cout << "Surface Distributed Forces...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                ElementPointer element = mesh->element(elNum);
                Point3D A = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(0))));
                Point3D B = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(1))));
                Point3D C = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(2))));
                DoubleMatrix lambda = cosinuses(A, B, C);
                DoubleMatrix lambdaT = lambda.transpose();
                double x[elementNodes];
                double y[elementNodes];
                // извлечение координат узлов
                for (UInteger i = 0; i < elementNodes; i++)
                {
                    Point3D point = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(i))));
                    Point3D pp = point - A;
                    x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
                    y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
                }

                double element_force[elementNodes];
                for (UInteger i = 0; i < elementNodes; i++)
                    element_force[i] = 0.0;

                FemCondition::FemDirection dir = (*condition)->direction();

                for (int ig = 0; ig < gaussPoints; ig++)
                {
                    double xi = gxi(ig);
                    double eta = geta(ig);
                    double w = gweight(ig);
                    // значения функций формы
                    DoubleVector N(elementNodes);
                    // значения производных функций формы
                    DoubleVector dNdX(elementNodes);
                    DoubleVector dNdY(elementNodes);
                    // якобиан
                    double jacobian = 1.0;
                    if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
                    {
                        jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
                    {
                        jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    double xLocal = 0.0;
                    double yLocal = 0.0;
                    for (UInteger i = 0; i < elementNodes; i++)
                    {
                        xLocal += x[i] * N[i];
                        yLocal += y[i] * N[i];
                    }
                    Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
                               lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
                               lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
                    // вычисление объемных сил
                    Point3D pLocal = pl + A;
                    double fLocal = (*condition)->value(&pLocal);

                    for (unsigned int i = 0; i < elementNodes; i++)
                    {
                        element_force[i] += N[i] * jacobian * w * fLocal;
                    }
                } // ig
                // ансамбль объемных сил
                for (UInteger i = 0 ; i < elementNodes; i++)
                {
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                    {
                        force(element->vertexNode(i)) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                    {
                        force(element->vertexNode(i) + nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                    {
                        force(element->vertexNode(i) + 2UL * nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                    {
                        force(element->vertexNode(i) + 3UL * nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                    {
                        force(element->vertexNode(i) + 4UL * nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                    {
                        force(element->vertexNode(i) + 5UL * nodesCount) += element_force[i];
                    }
                }
            } //for elNum
        }
    } // iterator

    //учет условий закрепления
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::INITIAL_VALUE)
        {
            // учет граничных условий
            std::cout << "Boundary Conditions...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        setInitialNodalValue(global, force, i, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        setInitialNodalValue(global, force, i + nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        setInitialNodalValue(global, force, i + 2UL * nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                        setInitialNodalValue(global, force, i + 3UL * nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                        setInitialNodalValue(global, force, i + 4UL * nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                        setInitialNodalValue(global, force, i + 5UL * nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator
//    std::cout << "force: " << force.size() << std::endl;
//    force.print(' ');

    DoubleVector displacement = solve(global, force);
//    DoubleVector displacement = global.cholesky(force);

    std::vector<double> xxx(nodesCount);
    std::vector<double> yyy(nodesCount);
    std::vector<double> zzz(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    std::vector<double> theta_z(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        xxx[i] = displacement[i];
        yyy[i] = displacement[i + nodesCount];
        zzz[i] = displacement[i + 2UL * nodesCount];
        theta_x[i] = displacement[i + 3UL * nodesCount];
        theta_y[i] = displacement[i + 4UL * nodesCount];
        theta_z[i] = displacement[i + 5UL * nodesCount];
    }
    nodeValues_.push_back(NamedVector("U", xxx));
    nodeValues_.push_back(NamedVector("V", yyy));
    nodeValues_.push_back(NamedVector("W", zzz));
    nodeValues_.push_back(NamedVector("Theta X", theta_x));
    nodeValues_.push_back(NamedVector("Theta Y", theta_y));
    nodeValues_.push_back(NamedVector("Theta Z", theta_z));

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount);
    std::vector<double> SigmaY(nodesCount);
    std::vector<double> SigmaZ(nodesCount);
    std::vector<double> TauXY(nodesCount);
    std::vector<double> TauXZ(nodesCount);
    std::vector<double> TauYZ(nodesCount);
    std::vector<double> mises(nodesCount);

    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] = 0.0;
        SigmaY[i] = 0.0;
        SigmaZ[i] = 0.0;
        TauXY[i] = 0.0;
        TauXZ[i] = 0.0;
        TauYZ[i] = 0.0;
        mises[i] = 0.0;
    }

    double xi[elementNodes];
    double eta[elementNodes];
    if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
    {
        xi[0] = 0.0; eta[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) !=NULL)
    {
        xi[0] = -1.0; eta[0] = -1.0;
        xi[1] =  1.0; eta[1] = -1.0;
        xi[2] =  1.0; eta[2] =  1.0;
        xi[3] = -1.0; eta[3] =  1.0;
    }

    std::cout << "Stresses Recovery...";
    progressBar.restart(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

        ElementPointer element = mesh->element(elNum);
        Point3D A = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(0))));
        Point3D B = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(1))));
        Point3D C = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(2))));
        DoubleMatrix lambda = cosinuses(A, B, C);
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);

        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        DoubleMatrix lambdaT = lambda.transpose();

        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов
        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
        }

        for (UInteger inode = 0; inode < elementNodes; inode++)
        {
            DoubleVector N(elementNodes);
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            DoubleVector dis((size_type)(freedom_ * elementNodes), 0.0);
            DoubleVector sigma_membrane((size_type)3, 0.0);
            DoubleVector sigma_plate((size_type)3, 0.0);
            DoubleVector sigma((size_type)3, 0.0);
            DoubleVector tau((size_type)2, 0.0);

            // якобиан
            if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
            {
                isoTriangle3(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
            {
                isoQuad4(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bm(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);
            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bm(0, i * freedom_) = dNdX(i);
                                                    Bm(1, i * freedom_ + 1) = dNdY(i);
                Bm(2, i * freedom_) = dNdY(i);      Bm(2, i * freedom_ + 1) = dNdX(i);

                Bf(0, i * freedom_ + 3) = dNdX(i);
                                                    Bf(1, i * freedom_ + 4) = dNdY(i);
                Bf(2, i * freedom_ + 3) = dNdY(i);  Bf(2, i * freedom_ + 4) = dNdX(i);

                Bc(0, i * freedom_ + 2) = dNdX(i);  Bc(0, i * freedom_ + 3) = N(i);
                Bc(1, i * freedom_ + 2) = dNdY(i);  Bc(1, i * freedom_ + 4) = N(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(freedom_ * i) = displacement[element->vertexNode(i)];
                dis(freedom_ * i + 1) = displacement[element->vertexNode(i) + nodesCount];
                dis(freedom_ * i + 2) = displacement[element->vertexNode(i) + 2UL * nodesCount];
                dis(freedom_ * i + 3) = displacement[element->vertexNode(i) + 3UL * nodesCount];
                dis(freedom_ * i + 4) = displacement[element->vertexNode(i) + 4UL * nodesCount];
                dis(freedom_ * i + 5) = displacement[element->vertexNode(i) + 5UL * nodesCount];
            }

            DoubleVector dLocal = T * dis;

            sigma_membrane = (D * Bm) * dLocal;
            sigma_plate = thickness / 2.0 * ((D * Bf) * dLocal);
            sigma[0] = sigma_membrane[0] + sigma_plate[0];
            sigma[1] = sigma_membrane[1] + sigma_plate[1];
            sigma[2] = sigma_membrane[2] + sigma_plate[2];
            tau = (Dc * Bc) * dLocal;

            DoubleMatrix localSigma(3, 3, 0.0);
            localSigma(0, 0) = sigma(0); localSigma(0, 1) = sigma(2); localSigma(0, 2) = tau(0);
            localSigma(1, 0) = sigma(2); localSigma(1, 1) = sigma(1); localSigma(1, 2) = tau(1);
            localSigma(2, 0) = tau(0);   localSigma(2, 1) = tau(1);   localSigma(2, 2) = 0.0;

            DoubleMatrix SG = lambdaT * localSigma * lambda;

            double von = sqrt(0.5) * sqrt((SG(0, 0) - SG(1, 1)) * (SG(0, 0) - SG(1, 1)) +
                                          (SG(1, 1) - SG(2, 2)) * (SG(1, 1) - SG(2, 2)) +
                                          (SG(2, 2) - SG(0, 0)) * (SG(2, 2) - SG(0, 0)) +
                                          6.0 * (SG(0, 1) * SG(0, 1) + SG(1, 2) * SG(1, 2) + SG(0, 2) * SG(0, 2)));
            SigmaX[element->vertexNode(inode)] += SG(0, 0);
            SigmaY[element->vertexNode(inode)] += SG(1, 1);
            SigmaZ[element->vertexNode(inode)] += SG(2, 2);
            TauXY[element->vertexNode(inode)] += SG(0, 1);
            TauXZ[element->vertexNode(inode)] += SG(0, 2);
            TauYZ[element->vertexNode(inode)] += SG(1, 2);
            mises[element->vertexNode(inode)] += von;
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] /= (double)mesh->adjacentCount(i);
        SigmaY[i] /= (double)mesh->adjacentCount(i);
        SigmaZ[i] /= (double)mesh->adjacentCount(i);
        TauXY[i] /= (double)mesh->adjacentCount(i);
        TauXZ[i] /= (double)mesh->adjacentCount(i);
        TauYZ[i] /= (double)mesh->adjacentCount(i);
        mises[i] /= (double)mesh->adjacentCount(i);
    }

    nodeValues_.push_back(NamedVector("Sigma X", SigmaX));
    nodeValues_.push_back(NamedVector("Sigma Y", SigmaY));
    nodeValues_.push_back(NamedVector("Sigma Z", SigmaZ));
    nodeValues_.push_back(NamedVector("Tau XY", TauXY));
    nodeValues_.push_back(NamedVector("Tau XZ", TauXZ));
    nodeValues_.push_back(NamedVector("Tau YZ", TauYZ));
    nodeValues_.push_back(NamedVector("von Mises", mises));
}

MindlinShellBending::MindlinShellBending(Mesh3D *mesh, const std::vector<double> &thickness, const std::vector<ElasticMatrix> &elasticMatrix, std::list<FemCondition *> conditions) :
    Fem2D(mesh, 6)
{
    const double kappa = 5.0 / 6.0;

    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 5; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
    {
        gaussPoints = line_count * line_count;
        gxi.resize(gaussPoints);
        geta.resize(gaussPoints);
        gweight.resize(gaussPoints);
        elementNodes = 4;
        for (int i = 0; i < line_count; i++)
        {
            for (int j = 0; j < line_count; j++)
            {
                gxi[i * line_count + j] = line_points[i];
                geta[i * line_count + j] = line_points[j];
                gweight[i * line_count + j] = line_weights[i] * line_weights[j];
            }
        }
    }

    unsigned layers_count = thickness.size();
    std::vector<DoubleMatrix> D(layers_count); // матрица упругости
    std::vector<DoubleMatrix> Dc(layers_count);
    double H = 0.0; // толщина

    for (unsigned i = 0; i < layers_count; i++)
    {
        D[i] = elasticMatrix[i].D();
        Dc[i].resize(2, 2);
        Dc[i] (0, 1) = Dc[i] (1, 0) = 0.0;
        Dc[i] (0, 0) = D[i] (2, 2);     Dc[i] (1, 1) = D[i] (2, 2);
        std::cout << "D[" << i << "]:" << std::endl;
        D[i].print();
        std::cout << "Dc["<< i << "]:" << std::endl;
        Dc[i].print();

        H += thickness[i];
    }
    std::cout << "Thickness: " << H << std::endl;

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    MappedDoubleMatrix global (dimension_); // глобальная матрица жесткости
    DoubleVector force(dimension_, 0.0); // вектор сил

    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {
        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        ElementPointer element = mesh->element(elNum);
        Point3D A = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(0))));
        Point3D B = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(1))));
        Point3D C = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(2))));
        DoubleMatrix lambda = cosinuses(A, B, C);
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов

        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
        }

        for (int ig = 0; ig < gaussPoints; ig++)
        {
            double xi = gxi(ig);
            double eta = geta(ig);
            double w = gweight(ig);
            // значения функций формы
            DoubleVector N(elementNodes);
            // значения производных функций формы
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            // якобиан
            double jacobian = 1.0;
            if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
            {
                jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
            {
                jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bm(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);

            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bm(0, i * freedom_) = dNdX(i);
                                                    Bm(1, i * freedom_ + 1) = dNdY(i);
                Bm(2, i * freedom_) = dNdY(i);      Bm(2, i * freedom_ + 1) = dNdX(i);

                Bf(0, i * freedom_ + 3) = dNdX(i);
                                                    Bf(1, i * freedom_ + 4) = dNdY(i);
                Bf(2, i * freedom_ + 3) = dNdY(i);  Bf(2, i * freedom_ + 4) = dNdX(i);

                Bc(0, i * freedom_ + 2) = dNdX(i);  Bc(0, i * freedom_ + 3) = N(i);
                Bc(1, i * freedom_ + 2) = dNdY(i);  Bc(1, i * freedom_ + 4) = N(i);
            }

            double z0 = -H / 2.0;
            double z1 = 0.0;
            for (unsigned i = 0; i < layers_count; i++)
            {
                z1 = z0 + thickness[i];
                local += jacobian * w * (z1 - z0) * (Bm.transpose() * D[i] * Bm);
                local += jacobian * w * (z1*z1 - z0*z0) / 2.0 * (Bm.transpose() * D[i] * Bf);
                local += jacobian * w * (z1*z1 - z0*z0) / 2.0 * (Bf.transpose() * D[i] * Bm);
                local += jacobian * w * (z1*z1*z1 - z0*z0*z0) / 3.0 * (Bf.transpose() * D[i] * Bf);
                local += jacobian * w * kappa * (z1 - z0) * (Bc.transpose() * Dc[i] * Bc);
                z0 = z1;
            }
        } // ig

        DoubleMatrix surf = T.transpose() * local * T;
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (UInteger i = 0; i < elementNodes * freedom_; i++)
        {
            index_i = element->vertexNode(i / freedom_) + (i % freedom_) * nodesCount;

            for (UInteger j = i; j < elementNodes * freedom_; j++)
            {
                index_j = element->vertexNode(j / freedom_) + (j % freedom_) * nodesCount;
                global(index_i, index_j) += surf(i, j);
                if (index_i != index_j) global(index_j, index_i) = global(index_i, index_j);
            } // for j
        } // for i
    } //for elNum

    // Учет сил
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::NODAL_FORCE)
        {
            // узловые нагрузки
            std::cout << "Nodal Forces...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    double f = (*condition)->value(point);
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(i) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(i + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(i + nodesCount + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                        force(i + nodesCount + nodesCount + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                        force(i + nodesCount + nodesCount + nodesCount + nodesCount) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                        force(i + nodesCount + nodesCount + nodesCount + nodesCount + nodesCount) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Edge Distributed Forces...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;

                if (mesh->isBorderElement(elNum))
                {
                    ElementPointer element = mesh->element(elNum);
                    for (unsigned int i = 0; i < elementNodes; i++)
                    {
                        if ((mesh->nodeType(element->vertexNode(i)) == BORDER || mesh->nodeType(element->vertexNode(i)) == CHARACTER) &&
                                (mesh->nodeType(element->vertexNode(i + 1)) == BORDER || mesh->nodeType(element->vertexNode(i + 1)) == CHARACTER))
                        {
                            PointPointer point0 = mesh->node(element->vertexNode(i));
                            PointPointer point1 = mesh->node(element->vertexNode(i + 1));
                            FemCondition::FemDirection dir = (*condition)->direction();
                            if ((*condition)->isApplied(point0) && (*condition)->isApplied(point1))
                            {
                                Point3D p0(point0->x(), point0->y(), point0->z());
                                Point3D p1(point1->x(), point1->y(), point1->z());
                                double l = p0.distanceTo(p1);
                                double jacobian = l / 2.0;
                                double f0 = 0.0;
                                double f1 = 0.0;
                                for (int ixi = 0; ixi < line_count; ixi++)
                                {

                                    double xi = line_points(ixi);
                                    double w = line_weights(ixi);
                                    double N0 = (1.0 - xi) / 2.0;
                                    double N1 = (1.0 + xi) / 2.0;
                                    Point3D p = Point3D(p0.x() * N0 + p1.x() * N1,
                                                        p0.y() * N0 + p1.y() * N1,
                                                        p1.z() * N0 + p1.z() * N1);
                                    f0 += N0 * jacobian * w * (*condition)->value(&p);
                                    f1 += N1 * jacobian * w * (*condition)->value(&p);
                                } // for ixi
                                if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                                {
                                    force(element->vertexNode(i)) += f0;
                                    force(element->vertexNode(i + 1)) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                                {
                                    force(element->vertexNode(i) + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                                {
                                    force(element->vertexNode(i) + 2UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 2UL * nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                                {
                                    force(element->vertexNode(i) + 3UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 3UL * nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                                {
                                    force(element->vertexNode(i) + 4UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 4UL * nodesCount) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                                {
                                    force(element->vertexNode(i) + 5UL * nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + 5UL * nodesCount) += f1;
                                }
                            }
                        } // if
                    } // for i
                } // if
            } // for elNum
        }
        else if ((*condition)->type() == FemCondition::VOLUME_FORCE)
        {
            // объемные силы
            std::cout << "Surface Distributed Forces...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                ElementPointer element = mesh->element(elNum);
                Point3D A = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(0))));
                Point3D B = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(1))));
                Point3D C = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(2))));
                DoubleMatrix lambda = cosinuses(A, B, C);
                DoubleMatrix lambdaT = lambda.transpose();
                double x[elementNodes];
                double y[elementNodes];
                // извлечение координат узлов
                for (UInteger i = 0; i < elementNodes; i++)
                {
                    Point3D point = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(i))));
                    Point3D pp = point - A;
                    x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
                    y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
                }

                double element_force[elementNodes];
                for (UInteger i = 0; i < elementNodes; i++)
                    element_force[i] = 0.0;

                FemCondition::FemDirection dir = (*condition)->direction();

                for (int ig = 0; ig < gaussPoints; ig++)
                {
                    double xi = gxi(ig);
                    double eta = geta(ig);
                    double w = gweight(ig);
                    // значения функций формы
                    DoubleVector N(elementNodes);
                    // значения производных функций формы
                    DoubleVector dNdX(elementNodes);
                    DoubleVector dNdY(elementNodes);
                    // якобиан
                    double jacobian = 1.0;
                    if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
                    {
                        jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
                    {
                        jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    double xLocal = 0.0;
                    double yLocal = 0.0;
                    for (UInteger i = 0; i < elementNodes; i++)
                    {
                        xLocal += x[i] * N[i];
                        yLocal += y[i] * N[i];
                    }
                    Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
                               lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
                               lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
                    // вычисление объемных сил
                    Point3D pLocal = pl + A;
                    double fLocal = (*condition)->value(&pLocal);

                    for (unsigned int i = 0; i < elementNodes; i++)
                    {
                        element_force[i] += N[i] * jacobian * w * fLocal;
                    }
                } // ig
                // ансамбль объемных сил
                for (UInteger i = 0 ; i < elementNodes; i++)
                {
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                    {
                        force(element->vertexNode(i)) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                    {
                        force(element->vertexNode(i) + nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                    {
                        force(element->vertexNode(i) + 2UL * nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                    {
                        force(element->vertexNode(i) + 3UL * nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                    {
                        force(element->vertexNode(i) + 4UL * nodesCount) += element_force[i];
                    }
                    if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                    {
                        force(element->vertexNode(i) + 5UL * nodesCount) += element_force[i];
                    }
                }
            } //for elNum
        }
    } // iterator

    //учет условий закрепления
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::INITIAL_VALUE)
        {
            // учет граничных условий
            std::cout << "Boundary Conditions...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        setInitialNodalValue(global, force, i, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        setInitialNodalValue(global, force, i + nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        setInitialNodalValue(global, force, i + 2UL * nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::FOURTH)
                        setInitialNodalValue(global, force, i + 3UL * nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::FIFTH)
                        setInitialNodalValue(global, force, i + 4UL * nodesCount, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SIXTH)
                        setInitialNodalValue(global, force, i + 5UL * nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator
//    std::cout << "force: " << force.size() << std::endl;
//    force.print(' ');

    DoubleVector displacement = solve(global, force);
//    DoubleVector displacement = global.cholesky(force);

    std::vector<double> xxx(nodesCount);
    std::vector<double> yyy(nodesCount);
    std::vector<double> zzz(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    std::vector<double> theta_z(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        xxx[i] = displacement[i];
        yyy[i] = displacement[i + nodesCount];
        zzz[i] = displacement[i + 2UL * nodesCount];
        theta_x[i] = displacement[i + 3UL * nodesCount];
        theta_y[i] = displacement[i + 4UL * nodesCount];
        theta_z[i] = displacement[i + 5UL * nodesCount];
    }
    nodeValues_.push_back(NamedVector("U", xxx));
    nodeValues_.push_back(NamedVector("V", yyy));
    nodeValues_.push_back(NamedVector("W", zzz));
    nodeValues_.push_back(NamedVector("Theta X", theta_x));
    nodeValues_.push_back(NamedVector("Theta Y", theta_y));
    nodeValues_.push_back(NamedVector("Theta Z", theta_z));

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount);
    std::vector<double> SigmaY(nodesCount);
    std::vector<double> SigmaZ(nodesCount);
    std::vector<double> TauXY(nodesCount);
    std::vector<double> TauXZ(nodesCount);
    std::vector<double> TauYZ(nodesCount);
    std::vector<double> mises(nodesCount);

    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] = 0.0;
        SigmaY[i] = 0.0;
        SigmaZ[i] = 0.0;
        TauXY[i] = 0.0;
        TauXZ[i] = 0.0;
        TauYZ[i] = 0.0;
        mises[i] = 0.0;
    }

    double xi[elementNodes];
    double eta[elementNodes];
    if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
    {
        xi[0] = 0.0; eta[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) !=NULL)
    {
        xi[0] = -1.0; eta[0] = -1.0;
        xi[1] =  1.0; eta[1] = -1.0;
        xi[2] =  1.0; eta[2] =  1.0;
        xi[3] = -1.0; eta[3] =  1.0;
    }

    std::cout << "Stresses Recovery...";
    progressBar.restart(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

        ElementPointer element = mesh->element(elNum);
        Point3D A = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(0))));
        Point3D B = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(1))));
        Point3D C = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(2))));
        DoubleMatrix lambda = cosinuses(A, B, C);
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);

        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        DoubleMatrix lambdaT = lambda.transpose();

        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов
        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
        }

        for (UInteger inode = 0; inode < elementNodes; inode++)
        {
            DoubleVector N(elementNodes);
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            DoubleVector dis((size_type)(freedom_ * elementNodes), 0.0);
            DoubleVector sigma_membrane((size_type)3, 0.0);
            DoubleVector sigma_plate((size_type)3, 0.0);
            DoubleVector sigma((size_type)3, 0.0);
            DoubleVector tau((size_type)2, 0.0);

            // якобиан
            if (dynamic_cast<TriangleMesh3D*>(mesh) != NULL)
            {
                isoTriangle3(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh) != NULL)
            {
                isoQuad4(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bm(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);
            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bm(0, i * freedom_) = dNdX(i);
                                                    Bm(1, i * freedom_ + 1) = dNdY(i);
                Bm(2, i * freedom_) = dNdY(i);      Bm(2, i * freedom_ + 1) = dNdX(i);

                Bf(0, i * freedom_ + 3) = dNdX(i);
                                                    Bf(1, i * freedom_ + 4) = dNdY(i);
                Bf(2, i * freedom_ + 3) = dNdY(i);  Bf(2, i * freedom_ + 4) = dNdX(i);

                Bc(0, i * freedom_ + 2) = dNdX(i);  Bc(0, i * freedom_ + 3) = N(i);
                Bc(1, i * freedom_ + 2) = dNdY(i);  Bc(1, i * freedom_ + 4) = N(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(freedom_ * i) = displacement[element->vertexNode(i)];
                dis(freedom_ * i + 1) = displacement[element->vertexNode(i) + nodesCount];
                dis(freedom_ * i + 2) = displacement[element->vertexNode(i) + 2UL * nodesCount];
                dis(freedom_ * i + 3) = displacement[element->vertexNode(i) + 3UL * nodesCount];
                dis(freedom_ * i + 4) = displacement[element->vertexNode(i) + 4UL * nodesCount];
                dis(freedom_ * i + 5) = displacement[element->vertexNode(i) + 5UL * nodesCount];
            }

            DoubleVector dLocal = T * dis;

            sigma_membrane = (D[layers_count - 1] * Bm) * dLocal;
            sigma_plate = H / 2.0 * ((D[layers_count - 1] * Bf) * dLocal);
            sigma[0] = sigma_membrane[0] + sigma_plate[0];
            sigma[1] = sigma_membrane[1] + sigma_plate[1];
            sigma[2] = sigma_membrane[2] + sigma_plate[2];
            tau = (Dc[layers_count - 1] * Bc) * dLocal;

            DoubleMatrix localSigma(3, 3, 0.0);
            localSigma(0, 0) = sigma(0); localSigma(0, 1) = sigma(2); localSigma(0, 2) = tau(0);
            localSigma(1, 0) = sigma(2); localSigma(1, 1) = sigma(1); localSigma(1, 2) = tau(1);
            localSigma(2, 0) = tau(0);   localSigma(2, 1) = tau(1);   localSigma(2, 2) = 0.0;

            DoubleMatrix SG = lambdaT * localSigma * lambda;

            double von = sqrt(0.5) * sqrt((SG(0, 0) - SG(1, 1)) * (SG(0, 0) - SG(1, 1)) +
                                          (SG(1, 1) - SG(2, 2)) * (SG(1, 1) - SG(2, 2)) +
                                          (SG(2, 2) - SG(0, 0)) * (SG(2, 2) - SG(0, 0)) +
                                          6.0 * (SG(0, 1) * SG(0, 1) + SG(1, 2) * SG(1, 2) + SG(0, 2) * SG(0, 2)));
            SigmaX[element->vertexNode(inode)] += SG(0, 0);
            SigmaY[element->vertexNode(inode)] += SG(1, 1);
            SigmaZ[element->vertexNode(inode)] += SG(2, 2);
            TauXY[element->vertexNode(inode)] += SG(0, 1);
            TauXZ[element->vertexNode(inode)] += SG(0, 2);
            TauYZ[element->vertexNode(inode)] += SG(1, 2);
            mises[element->vertexNode(inode)] += von;
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] /= (double)mesh->adjacentCount(i);
        SigmaY[i] /= (double)mesh->adjacentCount(i);
        SigmaZ[i] /= (double)mesh->adjacentCount(i);
        TauXY[i] /= (double)mesh->adjacentCount(i);
        TauXZ[i] /= (double)mesh->adjacentCount(i);
        TauYZ[i] /= (double)mesh->adjacentCount(i);
        mises[i] /= (double)mesh->adjacentCount(i);
    }

    nodeValues_.push_back(NamedVector("Sigma X", SigmaX));
    nodeValues_.push_back(NamedVector("Sigma Y", SigmaY));
    nodeValues_.push_back(NamedVector("Sigma Z", SigmaZ));
    nodeValues_.push_back(NamedVector("Tau XY", TauXY));
    nodeValues_.push_back(NamedVector("Tau XZ", TauXZ));
    nodeValues_.push_back(NamedVector("Tau YZ", TauYZ));
    nodeValues_.push_back(NamedVector("von Mises", mises));
}

DoubleMatrix MindlinShellBending::cosinuses(const Point3D &A, const Point3D &B, const Point3D &C)
{
    DoubleMatrix lambda(3, 3, 0.0);
    Point3D AB = B - A;
    Point3D AC = C - A;
    Point3D N = AC.product(AB); // со scilab отличается порядок обхода узлов, значит меняется порядок произвденеия для нахождения нормали
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vx.product(Vz);
    lambda(0, 0) = Vx.x(); lambda(0, 1) = Vx.y(); lambda(0, 2) = Vx.z();
    lambda(1, 0) = Vy.x(); lambda(1, 1) = Vy.y(); lambda(1, 2) = Vy.z();
    lambda(2, 0) = Vz.x(); lambda(2, 1) = Vz.y(); lambda(2, 2) = Vz.z();
    return lambda;
}
