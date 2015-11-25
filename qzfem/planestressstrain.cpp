#include "planestressstrain.h"

#include <iostream>

#include "rowdoublematrix.h"
#include "mappeddoublematrix.h"

#include "consoleprogress.h"

PlaneStressStrain::PlaneStressStrain(Mesh2D *mesh,
                                     double thickness,
                                     const ElasticMatrix &elasticMatrix,
                                     std::list<FemCondition *> conditions) :
    Fem2D(mesh, 2)
{
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 2; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 3;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
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
        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов
        ElementPointer element = mesh->element(elNum);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
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
            if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
            {
                jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
            {
                jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix B(3, elementNodes * freedom_, 0.0);
            for (UInteger i = 0; i < elementNodes; i++)
            {
                B(0, i * freedom_) = dNdX(i);
                                                B(1, i * freedom_ + 1) = dNdY(i);
                B(2, i * freedom_) = dNdY(i);   B(2, i * freedom_ + 1) = dNdX(i);
            }

            local += jacobian * w * thickness * (B.transpose() * D * B);
        } // ig
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (UInteger i = 0; i < elementNodes * freedom_; i++)
        {
            index_i = element->vertexNode(i / freedom_) + (i % freedom_) * nodesCount;

            for (UInteger j = i; j < elementNodes * freedom_; j++)
            {
                index_j = element->vertexNode(j / freedom_) + (j % freedom_) * nodesCount;
                global(index_i, index_j) += local(i, j);
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
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Surface Forces...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;

                if (mesh->isBorderElement(elNum))
                {
                    ElementPointer element = mesh->element(elNum);
                    for (UInteger i = 0; i < elementNodes; i++)
                    {
                        if ((mesh->nodeType(element->vertexNode(i)) == BORDER || mesh->nodeType(element->vertexNode(i)) == CHARACTER) &&
                                (mesh->nodeType(element->vertexNode(i + 1)) == BORDER || mesh->nodeType(element->vertexNode(i + 1)) == CHARACTER))
                        {
                            PointPointer point0 = mesh->node(element->vertexNode(i));
                            PointPointer point1 = mesh->node(element->vertexNode(i + 1));
                            FemCondition::FemDirection dir = (*condition)->direction();
                            if ((*condition)->isApplied(point0) && (*condition)->isApplied(point1))
                            {
                                Point2D p0(point0->x(), point0->y());
                                Point2D p1(point1->x(), point1->y());
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
                                    Point2D p = Point2D(p0.x() * N0 + p1.x() * N1,
                                                             p0.y() * N0 + p1.y() * N1);
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
                            }
                        } // if
                    } // for i
                } // if
            } // for elNum
        }
        else if ((*condition)->type() == FemCondition::VOLUME_FORCE)
        {
            // объемные силы
            std::cout << "Volume Forces...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                double x[elementNodes];
                double y[elementNodes];
                double vForce[elementNodes]; // значения объемных сил в узлах
                // извлечение координат узлов
                ElementPointer element = mesh->element(elNum);
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    PointPointer point = mesh->node(element->vertexNode(i));
                    x[i] = point->x();
                    y[i] = point->y();
                    vForce[i] = 0.0;
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
                    double xLocal = 0.0;
                    double yLocal = 0.0;
                    // якобиан
                    double jacobian = 1.0;
                    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
                    {
                        jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
                    {
                        jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    for (UInteger i = 0; i < elementNodes; i++)
                    {
                        xLocal += x[i] * N[i];
                        yLocal += y[i] * N[i];
                    }
                    // вычисление объемных сил
                    Point2D pLocal(xLocal, yLocal);
                    double fLocal = (*condition)->value(&pLocal);
                    for (UInteger i = 0; i < elementNodes; i++)
                    {
                        vForce[i] = vForce[i] + (N[i] * jacobian * w) * fLocal;
                    }
                } // ig
                FemCondition::FemDirection dir = (*condition)->direction();
                // ансамбль объемных сил
                for (UInteger i = 0 ; i < elementNodes; i++)
                {
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(element->vertexNode(i)) += vForce[i];
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(element->vertexNode(i) + nodesCount) += vForce[i];
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
                }
                ++progressBar;
            } // for i
        }
    } // iterator

    // решение СЛАУ
    DoubleVector displacement = solve(global, force);
    std::vector<double> u(nodesCount);
    std::vector<double> v(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        u[i] = displacement[i];
        v[i] = displacement[i + nodesCount];
    }
    nodeValues_.push_back(NamedVector("U", u));
    nodeValues_.push_back(NamedVector("V", v));

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount);
    std::vector<double> SigmaY(nodesCount);
    std::vector<double> TauXY(nodesCount);

    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] = 0.0;
        SigmaY[i] = 0.0;
        TauXY[i] = 0.0;
    }

    double xi[elementNodes];
    double eta[elementNodes];
    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
    {
        xi[0] = 0.0; eta[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0;
    }
    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) !=NULL)
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

        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов
        ElementPointer element = mesh->element(elNum);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
        }

        for (UInteger inode = 0; inode < elementNodes; inode++)
        {
            DoubleVector N(elementNodes);
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            DoubleMatrix dis((size_type)(freedom_ * elementNodes), (size_type)1);
            DoubleMatrix sigma((size_type)3, (size_type)1);
            if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
            {
                isoTriangle3(x[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
            {
                isoQuad4(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix B(3, elementNodes * freedom_, 0.0);
            for (UInteger i = 0; i < elementNodes; i++)
            {
                B(0, i * freedom_) = dNdX(i);
                                                B(1, i * freedom_ + 1) = dNdY(i);
                B(2, i * freedom_) = dNdY(i);   B(2, i * freedom_ + 1) = dNdX(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(i * freedom_, 0) = displacement[element->vertexNode(i)];
                dis(i * freedom_ + 1, 0) = displacement[element->vertexNode(i) + nodesCount];
            }

            sigma = (D * B) * dis;
            SigmaX[element->vertexNode(inode)] += sigma(0, 0);
            SigmaY[element->vertexNode(inode)] += sigma(1, 0);
            TauXY[element->vertexNode(inode)] += sigma(2, 0);
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] /= (double)mesh->adjacentCount(i);
        SigmaY[i] /= (double)mesh->adjacentCount(i);
        TauXY[i] /= (double)mesh->adjacentCount(i);
    }
    nodeValues_.push_back(NamedVector("Sigma X", SigmaX));
    nodeValues_.push_back(NamedVector("Sigma Y", SigmaY));
    nodeValues_.push_back(NamedVector("Tau XY", TauXY));
}
