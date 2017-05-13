#include "mindlinplatebending.h"

#include "consoleprogress.h"
#include "rowdoublematrix.h"

#include <math.h>

MindlinPlateBending::MindlinPlateBending(Mesh2D *mesh,
                                         double thickness,
                                         const DoubleMatrix &planeStressMatrix,
                                         const std::list<FemCondition *> &conditions) :
    Fem2D(mesh, 3, conditions)
{
    thickness_ = thickness;
    D_ = planeStressMatrix;
}

void MindlinPlateBending::buildGlobalMatrix()
{
    const double kappa = 5.0 / 6.0;
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 3; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
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
    DoubleMatrix D = D_; // матрица упругости
    DoubleMatrix Dc(2, 0.0);
    Dc(0, 0) = D(2, 2); Dc(1, 1) = D(2, 2);
    std::cout << "D" << std::endl;
    D.print();
    std::cout << "Dc:" << std::endl;
    Dc.print();
    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        // извлечение координат узлов
        ElementPointer element = mesh_->element(elNum);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh_->node(element->vertexNode(i));
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
            if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
            {
                jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
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

            local += jacobian * w * thickness_*thickness_*thickness_ / 12.0 * (Bf.transpose() * D * Bf);
            local += jacobian * w * kappa * thickness_ * (Bc.transpose() * Dc * Bc);
        } // ig
        // Ансамблирование
        assembly(element, local);
    } //for elNum
}

void MindlinPlateBending::buildGlobalVector()
{
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 3; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
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
    // Учет сил
    for (std::list<FemCondition *>::iterator condition = conditions_.begin(); condition != conditions_.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::NODAL_FORCE)
        {
            // узловые нагрузки
            std::cout << "Nodal Forces...";
            ConsoleProgress progressBar(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh_->node(i);
                if ((*condition)->isApplied(point))
                {
                    double f = (*condition)->value(point);
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force_(freedom_ * i) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force_(freedom_ * i + 1) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force_(freedom_ * i + 2) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Edge Distributed Forces...";
            ConsoleProgress progressBar(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;
                ElementPointer element = mesh_->element(elNum);
                if (mesh_->isBorderElement(element))
                {
                    for (unsigned int i = 0; i < elementNodes; i++)
                    {
                        if ((mesh_->nodeType(element->vertexNode(i)) == BORDER || mesh_->nodeType(element->vertexNode(i)) == CHARACTER) &&
                                (mesh_->nodeType(element->vertexNode(i + 1)) == BORDER || mesh_->nodeType(element->vertexNode(i + 1)) == CHARACTER))
                        {
                            PointPointer point0 = mesh_->node(element->vertexNode(i));
                            PointPointer point1 = mesh_->node(element->vertexNode(i + 1));
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
                                    force_(freedom_ * element->vertexNode(i)) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1)) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 1) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 1) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 2) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 2) += f1;
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
            ConsoleProgress progressBar(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                DoubleVector x(elementNodes);
                DoubleVector y(elementNodes);
                double vForce[elementNodes]; // значения объемных сил в узлах
                // извлечение координат узлов
                ElementPointer element = mesh_->element(elNum);
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    PointPointer point = mesh_->node(element->vertexNode(i));
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
                    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
                    {
                        jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
                    {
                        jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    // вычисление объемных сил
                    Point2D pLocal(x * N, y * N);
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
                        force_(freedom_ * element->vertexNode(i)) += vForce[i];
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force_(freedom_ * element->vertexNode(i) + 1) += vForce[i];
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force_(freedom_ * element->vertexNode(i) + 2) += vForce[i];
                }
            } //for elNum
        }
    } // iterator
}

void MindlinPlateBending::processSolution(const DoubleVector &displacement)
{
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
    {
        elementNodes = 3;
    }
    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
    {
        elementNodes = 4;
    }
    DoubleMatrix D = D_; // матрица упругости
    DoubleMatrix Dc(2, 0.0);
    Dc(0, 0) = D(2, 2); Dc(1, 1) = D(2, 2);
    std::cout << "D" << std::endl;
    D.print();
    std::cout << "Dc:" << std::endl;
    Dc.print();
    std::vector<double> w(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        w[i] = displacement[freedom_ * i];
        theta_x[i] = displacement[freedom_ * i + 1];
        theta_y[i] = displacement[freedom_ * i + 2];
    }

    mesh_->addDataVector("W", w);
    mesh_->addDataVector("Theta X", theta_x);
    mesh_->addDataVector("Theta Y", theta_y);

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount, 0.0);
    std::vector<double> SigmaY(nodesCount, 0.0);
    std::vector<double> TauXY(nodesCount, 0.0);
    std::vector<double> TauXZ(nodesCount, 0.0);
    std::vector<double> TauYZ(nodesCount, 0.0);
    std::vector<double> mises(nodesCount, 0.0);

    double xi[elementNodes];
    double eta[elementNodes];
    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
    {
        xi[0] = 0.0; eta[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0;
    }
    else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) !=NULL)
    {
        xi[0] = -1.0; eta[0] = -1.0;
        xi[1] =  1.0; eta[1] = -1.0;
        xi[2] =  1.0; eta[2] =  1.0;
        xi[3] = -1.0; eta[3] =  1.0;
    }

    std::cout << "Stresses Recovery...";
    ConsoleProgress progressBar(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        // извлечение координат узлов
        ElementPointer element = mesh_->element(elNum);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh_->node(element->vertexNode(i));
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
            if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
            {
                isoTriangle3(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
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
                dis(freedom_ * i, 0) = displacement[freedom_ * element->vertexNode(i)];
                dis(freedom_ * i + 1, 0) = displacement[freedom_ * element->vertexNode(i) + 1];
                dis(freedom_ * i + 2, 0) = displacement[freedom_ * element->vertexNode(i) + 2];
            }

            sigma = (thickness_ / 2.0) * ((D * Bf) * dis);
            tau = ((Dc * Bc) * dis);
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
        SigmaX[i] /= (double)mesh_->adjacentCount(i);
        SigmaY[i] /= (double)mesh_->adjacentCount(i);
        TauXY[i] /= (double)mesh_->adjacentCount(i);
        TauXZ[i] /= (double)mesh_->adjacentCount(i);
        TauYZ[i] /= (double)mesh_->adjacentCount(i);
        mises[i] /= (double)mesh_->adjacentCount(i);
    }
    mesh_->addDataVector("Sigma X", SigmaX);
    mesh_->addDataVector("Sigma Y", SigmaY);
    mesh_->addDataVector("Tau XY", TauXY);
    mesh_->addDataVector("Tau XZ", TauXZ);
    mesh_->addDataVector("Tau YZ", TauYZ);
    mesh_->addDataVector("von Mises", mises);
}
