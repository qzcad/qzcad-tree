#include "mindlinplatebending.h"

#include "consoleprogress.h"
#include "rowdoublematrix.h"

#include <math.h>

MindlinPlateBending::MindlinPlateBending(Mesh2D *mesh,
                                         double thickness,
                                         const ElasticMatrix &elasticMatrix,
                                         std::list<FemCondition *> conditions) : Fem2D(mesh, 3)
{
    const double kappa = 5.0 / 6.0;

    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 3; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
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
    DoubleMatrix Dc(2, 0.0);
    Dc(0, 0) = D(2, 2); Dc(1, 1) = D(2, 2);

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
                jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
            {
                jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);

            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bf(0, i * freedom_ + 1) = dNdX(i);
                Bf(1, i * freedom_ + 2) = dNdY(i);
                Bf(2, i * freedom_ + 1) = dNdY(i);  Bf(2, i * freedom_ + 2) = dNdX(i);

                Bc(0, i * freedom_) = dNdX(i);  Bc(0, i * freedom_ + 1) = N(i);
                Bc(1, i * freedom_) = dNdY(i);  Bc(1, i * freedom_ + 2) = N(i);
            }

            local += jacobian * w * thickness*thickness*thickness / 12.0 * (Bf.transpose() * D * Bf);
            local += jacobian * w * kappa * thickness * (Bc.transpose() * Dc * Bc);
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(i + nodesCount + nodesCount) += f;
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
                                if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                                {
                                    force(element->vertexNode(i) + nodesCount + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount + nodesCount) += f1;
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
                    // якобиан
                    double jacobian = 1.0;
                    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
                    {
                        jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
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
                    // вычисление объемных сил
                    Point2D pLocal(xLocal, yLocal);
                    double fLocal = (*condition)->value(&pLocal);
                    for (unsigned int i = 0; i < elementNodes; i++)
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(element->vertexNode(i) + nodesCount + nodesCount) += vForce[i];
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
                        setInitialNodalValue(global, force, i + nodesCount + nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator

    DoubleVector displacement = solve(global, force);
    std::vector<double> w(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        w[i] = displacement[i];
        theta_x[i] = displacement[i + nodesCount];
        theta_y[i] = displacement[i + nodesCount + nodesCount];
    }
    nodeValues_.push_back(NamedVector("W", w));
    nodeValues_.push_back(NamedVector("Theta X", theta_x));
    nodeValues_.push_back(NamedVector("Theta Y", theta_y));

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount);
    std::vector<double> SigmaY(nodesCount);
    std::vector<double> TauXY(nodesCount);
    std::vector<double> TauXZ(nodesCount);
    std::vector<double> TauYZ(nodesCount);
    std::vector<double> mises(nodesCount);

    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] = 0.0;
        SigmaY[i] = 0.0;
        TauXY[i] = 0.0;
        TauXZ[i] = 0.0;
        TauYZ[i] = 0.0;
        mises[i] = 0.0;
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
            DoubleMatrix tau((size_type)2, (size_type)1);

            // якобиан
            if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
            {
                isoTriangle3(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
            {
                isoQuad4(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);

            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bf(0, i * freedom_ + 1) = dNdX(i);
                Bf(1, i * freedom_ + 2) = dNdY(i);
                Bf(2, i * freedom_ + 1) = dNdY(i);  Bf(2, i * freedom_ + 2) = dNdX(i);

                Bc(0, i * freedom_) = dNdX(i);  Bc(0, i * freedom_ + 1) = N(i);
                Bc(1, i * freedom_) = dNdY(i);  Bc(1, i * freedom_ + 2) = N(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(3 * i, 0) = displacement[element->vertexNode(i)];
                dis(3 * i + 1, 0) = displacement[element->vertexNode(i) + nodesCount];
                dis(3 * i + 2, 0) = displacement[element->vertexNode(i) + nodesCount + nodesCount];
            }

            sigma = thickness / 2.0 * ((D * Bf) * dis);
            tau = ((Dc * Bc) * dis);
//            double von = sqrt( 0.5 * ((sigma(0,0) - sigma(1,0))*(sigma(0,0) - sigma(1,0)) +
//                          sigma(1,0)*sigma(1,0) +
//                          sigma(0,0)*sigma(0,0) +
//                          6.0 * (sigma(2, 0)*sigma(2, 0) + tau(0, 0)*tau(0, 0) + tau(1, 0)*tau(1, 0))) );
            double von = sqrt(sigma(0,0)*sigma(0,0) - sigma(0,0)*sigma(1,0) + sigma(1,0)*sigma(1,0) + 3.0 * sigma(2,0)*sigma(2,0));

            SigmaX[element->vertexNode(inode)] += sigma(0, 0);
            SigmaY[element->vertexNode(inode)] += sigma(1, 0);
            TauXY[element->vertexNode(inode)] += sigma(2, 0);
            TauXZ[element->vertexNode(inode)] += tau(0, 0);
            TauYZ[element->vertexNode(inode)] += tau(1, 0);
            mises[element->vertexNode(inode)] += von;
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] /= (double)mesh->adjacentCount(i);
        SigmaY[i] /= (double)mesh->adjacentCount(i);
        TauXY[i] /= (double)mesh->adjacentCount(i);
        TauXZ[i] /= (double)mesh->adjacentCount(i);
        TauYZ[i] /= (double)mesh->adjacentCount(i);
        mises[i] /= (double)mesh->adjacentCount(i);
    }
    nodeValues_.push_back(NamedVector("Sigma X", SigmaX));
    nodeValues_.push_back(NamedVector("Sigma Y", SigmaY));
    nodeValues_.push_back(NamedVector("Tau XY", TauXY));
    nodeValues_.push_back(NamedVector("Tau XZ", TauXZ));
    nodeValues_.push_back(NamedVector("Tau YZ", TauYZ));
    nodeValues_.push_back(NamedVector("von Mises", mises));
}

MindlinPlateBending::MindlinPlateBending(Mesh2D *mesh, const std::vector<double> &thickness, const std::vector<ElasticMatrix> &elasticMatrix, std::list<FemCondition *> conditions) :
    Fem2D(mesh, 3)
{
    const double kappa = 5.0 / 6.0;

    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 3; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
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
    unsigned layers_count = thickness.size();
    std::vector<DoubleMatrix> D(layers_count); // матрица упругости
    std::vector<DoubleMatrix> Dc(layers_count);
    double H = 0.0; // толщина пластинки

    for (unsigned i = 0; i < layers_count; i++)
    {
        D[i] = elasticMatrix[i].D();
        Dc[i].resize(2, 2);
        Dc[i] (0, 1) = Dc[i] (1, 0) = 0.0;
        Dc[i] (0, 0) = D[i] (2, 2);     Dc[i] (1, 1) = D[i] (2, 2);

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
                jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
            {
                jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);

            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bf(0, i * freedom_ + 1) = dNdX(i);
                Bf(1, i * freedom_ + 2) = dNdY(i);
                Bf(2, i * freedom_ + 1) = dNdY(i);  Bf(2, i * freedom_ + 2) = dNdX(i);

                Bc(0, i * freedom_) = dNdX(i);  Bc(0, i * freedom_ + 1) = N(i);
                Bc(1, i * freedom_) = dNdY(i);  Bc(1, i * freedom_ + 2) = N(i);
            }
            double z0 = -H / 2.0;
            double z1 = 0.0;
            for (unsigned i = 0; i < layers_count; i++)
            {
                z1 = z0 + thickness[i];
                local += (jacobian * w * (z1*z1*z1 - z0*z0*z0) / 3.0) * (Bf.transpose() * D[i] * Bf);
                local += (jacobian * w * kappa * (z1 - z0)) * (Bc.transpose() * Dc[i] * Bc);
                z0 = z1;
            }
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(i + nodesCount + nodesCount) += f;
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
                                if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                                {
                                    force(element->vertexNode(i) + nodesCount + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount + nodesCount) += f1;
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
                    // якобиан
                    double jacobian = 1.0;
                    if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
                    {
                        jacobian = 0.5 * isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
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
                    // вычисление объемных сил
                    Point2D pLocal(xLocal, yLocal);
                    double fLocal = (*condition)->value(&pLocal);
                    for (unsigned int i = 0; i < elementNodes; i++)
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(element->vertexNode(i) + nodesCount + nodesCount) += vForce[i];
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
                        setInitialNodalValue(global, force, i + nodesCount + nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator

    DoubleVector displacement = solve(global, force);
    std::vector<double> w(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        w[i] = displacement[i];
        theta_x[i] = displacement[i + nodesCount];
        theta_y[i] = displacement[i + nodesCount + nodesCount];
    }
    nodeValues_.push_back(NamedVector("W", w));
    nodeValues_.push_back(NamedVector("Theta X", theta_x));
    nodeValues_.push_back(NamedVector("Theta Y", theta_y));

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount);
    std::vector<double> SigmaY(nodesCount);
    std::vector<double> TauXY(nodesCount);
    std::vector<double> TauXZ(nodesCount);
    std::vector<double> TauYZ(nodesCount);
    std::vector<double> mises(nodesCount);

    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] = 0.0;
        SigmaY[i] = 0.0;
        TauXY[i] = 0.0;
        TauXZ[i] = 0.0;
        TauYZ[i] = 0.0;
        mises[i] = 0.0;
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
            DoubleMatrix tau((size_type)2, (size_type)1);

            // якобиан
            if (dynamic_cast<TriangleMesh2D*>(mesh) != NULL)
            {
                isoTriangle3(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh) != NULL)
            {
                isoQuad4(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            //
            DoubleMatrix Bf(3, freedom_ * elementNodes, 0.0);
            DoubleMatrix Bc(2, freedom_ * elementNodes, 0.0);

            for (UInteger i = 0; i < elementNodes; i++)
            {
                Bf(0, i * freedom_ + 1) = dNdX(i);
                Bf(1, i * freedom_ + 2) = dNdY(i);
                Bf(2, i * freedom_ + 1) = dNdY(i);  Bf(2, i * freedom_ + 2) = dNdX(i);

                Bc(0, i * freedom_) = dNdX(i);  Bc(0, i * freedom_ + 1) = N(i);
                Bc(1, i * freedom_) = dNdY(i);  Bc(1, i * freedom_ + 2) = N(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(3 * i, 0) = displacement[element->vertexNode(i)];
                dis(3 * i + 1, 0) = displacement[element->vertexNode(i) + nodesCount];
                dis(3 * i + 2, 0) = displacement[element->vertexNode(i) + nodesCount + nodesCount];
            }

            sigma = H / 2.0 * ((D[layers_count - 1] * Bf) * dis);
            tau = ((Dc[layers_count - 1] * Bc) * dis);
//            double von = sqrt( 0.5 * ((sigma(0,0) - sigma(1,0))*(sigma(0,0) - sigma(1,0)) +
//                          sigma(1,0)*sigma(1,0) +
//                          sigma(0,0)*sigma(0,0) +
//                          6.0 * (sigma(2, 0)*sigma(2, 0) + tau(0, 0)*tau(0, 0) + tau(1, 0)*tau(1, 0))) );
            double von = sqrt(sigma(0,0)*sigma(0,0) - sigma(0,0)*sigma(1,0) + sigma(1,0)*sigma(1,0) + 3.0 * sigma(2,0)*sigma(2,0));

            SigmaX[element->vertexNode(inode)] += sigma(0, 0);
            SigmaY[element->vertexNode(inode)] += sigma(1, 0);
            TauXY[element->vertexNode(inode)] += sigma(2, 0);
            TauXZ[element->vertexNode(inode)] += tau(0, 0);
            TauYZ[element->vertexNode(inode)] += tau(1, 0);
            mises[element->vertexNode(inode)] += von;
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] /= (double)mesh->adjacentCount(i);
        SigmaY[i] /= (double)mesh->adjacentCount(i);
        TauXY[i] /= (double)mesh->adjacentCount(i);
        TauXZ[i] /= (double)mesh->adjacentCount(i);
        TauYZ[i] /= (double)mesh->adjacentCount(i);
        mises[i] /= (double)mesh->adjacentCount(i);
    }
    nodeValues_.push_back(NamedVector("Sigma X (z = h / 2)", SigmaX));
    nodeValues_.push_back(NamedVector("Sigma Y (z = h / 2)", SigmaY));
    nodeValues_.push_back(NamedVector("Tau XY (z = h / 2)", TauXY));
    nodeValues_.push_back(NamedVector("Tau XZ (z = h / 2)", TauXZ));
    nodeValues_.push_back(NamedVector("Tau YZ (z = h / 2)", TauYZ));
    nodeValues_.push_back(NamedVector("von Mises (z = h / 2)", mises));
}
