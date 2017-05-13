#include "planestressstrain.h"

#include <iostream>
#include <math.h>

#include "rowdoublematrix.h"
#include "mappeddoublematrix.h"

#include "consoleprogress.h"

PlaneStressStrain::PlaneStressStrain(Mesh2D *mesh,
                                     double thickness,
                                     const DoubleMatrix &elasticMatrix,
                                     const std::list<FemCondition *> &conditions, double alphaT) :
    Fem2D(mesh, 2, conditions)
{
    thickness_ = thickness;
    D_ = elasticMatrix;
    alpha_ = alphaT;
}

void PlaneStressStrain::buildGlobalMatrix()
{
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 2; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 3;
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

    // температурные деформации
    DoubleVector epsilon0(3, 0.0);
    epsilon0(0) = alpha_;
    epsilon0(1) = alpha_;

    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        DoubleVector epsilonForce(freedom_ * elementNodes, 0.0);
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
            DoubleMatrix B(3, elementNodes * freedom_, 0.0);
            for (UInteger i = 0; i < elementNodes; i++)
            {
                B(0, i * freedom_) = dNdX(i);
                                                B(1, i * freedom_ + 1) = dNdY(i);
                B(2, i * freedom_) = dNdY(i);   B(2, i * freedom_ + 1) = dNdX(i);
            }

            local += jacobian * w * thickness_ * (B.transpose() * D_ * B);

            epsilonForce += jacobian * w * thickness_ * ((B.transpose() * D_) * epsilon0);
        } // ig
        // Ансамблирование
        assembly(element, local);
        for (UInteger i = 0; i < freedom_ * elementNodes; i++)
        {
            UInteger index = freedom_ * element->vertexNode(i / freedom_) + (i % freedom_);
            force_(index) += epsilonForce(i);
        }
    } //for elNum
}

void PlaneStressStrain::buildGlobalVector()
{
    // Учет сил
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 2; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 3;
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
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Surface Forces...";
            ConsoleProgress progressBar(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;
                ElementPointer element = mesh_->element(elNum);
                if (mesh_->isBorderElement(element))
                {
                    for (UInteger i = 0; i < elementNodes; i++)
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
                        force_(freedom_ * element->vertexNode(i)) += vForce[i];
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force_(freedom_ * element->vertexNode(i) + 1) += vForce[i];
                }
            } //for elNum
        }
    } // iterator
}

void PlaneStressStrain::processSolution(const DoubleVector &displacement)
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

    std::vector<double> u(nodesCount);
    std::vector<double> v(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        u[i] = displacement[freedom_ * i];
        v[i] = displacement[freedom_ * i + 1];
    }

    mesh_->addDataVector("X", u);
    mesh_->addDataVector("Y", v);

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount, 0.0);
    std::vector<double> SigmaY(nodesCount, 0.0);
    std::vector<double> TauXY(nodesCount, 0.0);
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
            if (dynamic_cast<TriangleMesh2D*>(mesh_) != NULL)
            {
                isoTriangle3(x[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh2D*>(mesh_) != NULL)
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
                dis(i * freedom_, 0) = displacement[freedom_ * element->vertexNode(i)];
                dis(i * freedom_ + 1, 0) = displacement[freedom_ * element->vertexNode(i) + 1];
            }

            sigma = (D_ * B) * dis;
            SigmaX[element->vertexNode(inode)] += sigma(0, 0);
            SigmaY[element->vertexNode(inode)] += sigma(1, 0);
            TauXY[element->vertexNode(inode)] += sigma(2, 0);
            mises[element->vertexNode(inode)] += sqrt(sigma(0,0)*sigma(0,0) - sigma(0,0)*sigma(1,0) + sigma(1,0)*sigma(1,0) + 3.0 * sigma(2,0)*sigma(2,0)); // general plane stress
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        SigmaX[i] /= (double)mesh_->adjacentCount(i);
        SigmaY[i] /= (double)mesh_->adjacentCount(i);
        TauXY[i] /= (double)mesh_->adjacentCount(i);
        mises[i] /= (double)mesh_->adjacentCount(i);
    }

    mesh_->addDataVector("Sigma X", SigmaX);
    mesh_->addDataVector("Sigma Y", SigmaY);
    mesh_->addDataVector("Tau XY", TauXY);
    mesh_->addDataVector("von Mises", mises);
}
