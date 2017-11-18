#include "mindlinshellbending.h"

#include <iostream>
#include "consoleprogress.h"
#include "rowdoublematrix.h"
#include <math.h>

MindlinShellBending::MindlinShellBending(Mesh3D *mesh,
                                         double thickness,
                                         const DoubleMatrix &planeStressMatrix,
                                         const std::list<FemCondition *> &conditions,
                                         double alphaT) :
    Fem2D(mesh, 6, conditions)
{
    thickness_ = thickness;
    D_ = planeStressMatrix;
    Dc_.resize(2, 2);
    Dc_(0, 0) = D_(2, 2); Dc_(0, 1) = 0.0;
    Dc_(1, 0) = 0.0; Dc_(1, 1) = D_(2, 2);
    alpha_ = alphaT;
}

MindlinShellBending::MindlinShellBending(Mesh3D *mesh,
                                         double thickness,
                                         const DoubleMatrix &D,
                                         const DoubleMatrix &Dc,
                                         const std::list<FemCondition *> &conditions,
                                         double alphaT) :
    Fem2D(mesh, 6, conditions)
{
    thickness_ = thickness;
    D_ = D;
    Dc_ = Dc;
    alpha_ = alphaT;
}

/*
MindlinShellBending::MindlinShellBending(Mesh3D *mesh,
                                         const std::vector<double> &thickness,
                                         const std::vector<DoubleMatrix> &planeStressMatrix,
                                         const std::list<FemCondition *> &conditions) :
    Fem2D(mesh, 6, conditions)
{


    force_ = evalForces(conditions);

    processInitialValues();

    DoubleVector displacement = solveLinearSystem();

    std::vector<double> xxx(nodesCount);
    std::vector<double> yyy(nodesCount);
    std::vector<double> zzz(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    std::vector<double> theta_z(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        xxx[i] = displacement[freedom_ * i];
        yyy[i] = displacement[freedom_ * i + 1UL];
        zzz[i] = displacement[freedom_ * i + 2UL];
        theta_x[i] = displacement[freedom_ * i + 3UL];
        theta_y[i] = displacement[freedom_ * i + 4UL];
        theta_z[i] = displacement[freedom_ * i + 5UL];
    }

    mesh_->addDataVector("X", xxx);
    mesh_->addDataVector("Y", yyy);
    mesh_->addDataVector("Z", zzz);
    mesh_->addDataVector("Theta X", theta_x);
    mesh_->addDataVector("Theta Y", theta_y);
    mesh_->addDataVector("Theta Z", theta_z);

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount, 0.0);
    std::vector<double> SigmaY(nodesCount, 0.0);
    std::vector<double> SigmaZ(nodesCount, 0.0);
    std::vector<double> TauXY(nodesCount, 0.0);
    std::vector<double> TauXZ(nodesCount, 0.0);
    std::vector<double> TauYZ(nodesCount, 0.0);
    std::vector<double> mises(nodesCount, 0.0);

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
                dis(freedom_ * i) = displacement[freedom_ * element->vertexNode(i)];
                dis(freedom_ * i + 1) = displacement[freedom_ * element->vertexNode(i) + 1UL];
                dis(freedom_ * i + 2) = displacement[freedom_ * element->vertexNode(i) + 2UL];
                dis(freedom_ * i + 3) = displacement[freedom_ * element->vertexNode(i) + 3UL];
                dis(freedom_ * i + 4) = displacement[freedom_ * element->vertexNode(i) + 4UL];
                dis(freedom_ * i + 5) = displacement[freedom_ * element->vertexNode(i) + 5UL];
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

    mesh_->addDataVector("Sigma X", SigmaX);
    mesh_->addDataVector("Sigma Y", SigmaY);
    mesh_->addDataVector("Sigma Z", SigmaZ);
    mesh_->addDataVector("Tau XY", TauXY);
    mesh_->addDataVector("Tau XZ", TauXZ);
    mesh_->addDataVector("Tau YZ", TauYZ);
    mesh_->addDataVector("von Mises", mises);
}
*/

void MindlinShellBending::buildGlobalMatrix()
{
    const double kappa = 5.0 / 6.0;
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 5; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
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

    std::cout << "D:" << std::endl;
    D_.print();
    std::cout << "Dc:" << std::endl;
    Dc_.print();
    std::cout << "Thickness: " << thickness_ << std::endl;

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
        ElementPointer element = mesh_->element(elNum);
        Point3D A = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(0))));
        Point3D B = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(1))));
        Point3D C = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(2))));
        DoubleMatrix lambda = cosinuses(A, B, C);
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        DoubleVector epsilonForce(freedom_ * elementNodes, 0.0);
        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        // извлечение координат узлов

        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
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
            if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
            {
                jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
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

            local += jacobian * w * thickness_ * (Bm.transpose() * D_ * Bm);
            local += jacobian * w * thickness_*thickness_*thickness_ / 12.0 * (Bf.transpose() * D_ * Bf);
            local += jacobian * w * kappa * thickness_ * (Bc.transpose() * Dc_ * Bc);

            epsilonForce += jacobian * w * thickness_ * ((Bm.transpose() * D_) * epsilon0);
        } // ig

        DoubleMatrix surf = T.transpose() * local * T;
        // Ансамблирование
        assembly(element, surf);
        DoubleVector localForce = T.transpose() * epsilonForce;
        for (UInteger i = 0; i < freedom_ * elementNodes; i++)
        {
            UInteger index = freedom_ * element->vertexNode(i / freedom_) + (i % freedom_);
            force_(index) += localForce(i);
        }
    } //for elNum

    std::cout << "force norm: " << force_.norm_2() << std::endl;

}

void MindlinShellBending::buildGlobalVector()
{
    UInteger nodesCount = mesh_->nodesCount();
    UInteger elementsCount = mesh_->elementsCount();
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 5; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
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
                    int dir = (*condition)->direction();
                    if (dir & FemCondition::FIRST)
                        force_(freedom_ * i) += f;
                    if (dir & FemCondition::SECOND)
                        force_(freedom_ * i + 1UL) += f;
                    if (dir & FemCondition::THIRD)
                        force_(freedom_ * i + 2UL) += f;
                    if (dir & FemCondition::FOURTH)
                        force_(freedom_ * i + 3UL) += f;
                    if (dir & FemCondition::FIFTH)
                        force_(freedom_ * i + 4UL) += f;
                    if (dir & FemCondition::SIXTH)
                        force_(freedom_ * i + 5UL) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки (в данном случае это нагрузки на ребра)
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
                            int dir = (*condition)->direction();
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
                                if (dir & FemCondition::FIRST)
                                {
                                    force_(freedom_ * element->vertexNode(i)) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1)) += f1;
                                }
                                if (dir & FemCondition::SECOND)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 1UL) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 1UL) += f1;
                                }
                                if (dir & FemCondition::THIRD)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 2UL) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 2UL) += f1;
                                }
                                if (dir & FemCondition::FOURTH)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 3UL) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 3UL) += f1;
                                }
                                if (dir & FemCondition::FIFTH)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 4UL) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 4UL) += f1;
                                }
                                if (dir & FemCondition::SIXTH)
                                {
                                    force_(freedom_ * element->vertexNode(i) + 5UL) += f0;
                                    force_(freedom_ * element->vertexNode(i + 1) + 5UL) += f1;
                                }
                            }
                        } // if
                    } // for i
                } // if
            } // for elNum
        }
        else if ((*condition)->type() == FemCondition::VOLUME_FORCE)
        {
            // объемные силы (в данном случае - распрделенные по поверхности оболочки)
            std::cout << "Surface Distributed Forces...";
            ConsoleProgress progressBar(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                ElementPointer element = mesh_->element(elNum);
                Point3D A = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(0))));
                Point3D B = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(1))));
                Point3D C = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(2))));
                DoubleMatrix lambda = cosinuses(A, B, C);
                DoubleMatrix lambdaT = lambda.transpose();
                DoubleVector x(elementNodes);
                DoubleVector y(elementNodes);
                // извлечение координат узлов
                for (UInteger i = 0; i < elementNodes; i++)
                {
                    Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
                    Point3D pp = point - A;
                    x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
                    y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
                }

                double element_force[elementNodes];
                for (UInteger i = 0; i < elementNodes; i++)
                    element_force[i] = 0.0;

                int dir = (*condition)->direction();

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
                    if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
                    {
                        jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
                    {
                        jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                    }
                    double xLocal = x * N;
                    double yLocal = y * N;
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
                    if (dir & FemCondition::FIRST)
                    {
                        force_(freedom_ * element->vertexNode(i)) += element_force[i];
                    }
                    if (dir & FemCondition::SECOND)
                    {
                        force_(freedom_ * element->vertexNode(i) + 1UL) += element_force[i];
                    }
                    if (dir & FemCondition::THIRD)
                    {
                        force_(freedom_ * element->vertexNode(i) + 2UL) += element_force[i];
                    }
                    if (dir & FemCondition::FOURTH)
                    {
                        force_(freedom_ * element->vertexNode(i) + 3UL) += element_force[i];
                    }
                    if (dir & FemCondition::FIFTH)
                    {
                        force_(freedom_ * element->vertexNode(i) + 4UL) += element_force[i];
                    }
                    if (dir & FemCondition::SIXTH)
                    {
                        force_(freedom_ * element->vertexNode(i) + 5UL) += element_force[i];
                    }
                }
            } // for elNum
        }
    } // iterator
}

void MindlinShellBending::processSolution(const DoubleVector &displacement)
{
    UInteger nodesCount = mesh_->nodesCount();
    UInteger elementsCount = mesh_->elementsCount();
    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
    {
        elementNodes = 3;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
    {
        elementNodes = 4;
    }
    std::vector<double> xxx(nodesCount);
    std::vector<double> yyy(nodesCount);
    std::vector<double> zzz(nodesCount);
    std::vector<double> theta_x(nodesCount);
    std::vector<double> theta_y(nodesCount);
    std::vector<double> theta_z(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        xxx[i] = displacement[freedom_ * i];
        yyy[i] = displacement[freedom_ * i + 1UL];
        zzz[i] = displacement[freedom_ * i + 2UL];
        theta_x[i] = displacement[freedom_ * i + 3UL];
        theta_y[i] = displacement[freedom_ * i + 4UL];
        theta_z[i] = displacement[freedom_ * i + 5UL];
    }

    mesh_->addDataVector("X", xxx);
    mesh_->addDataVector("Y", yyy);
    mesh_->addDataVector("Z", zzz);
    mesh_->addDataVector("Theta X", theta_x);
    mesh_->addDataVector("Theta Y", theta_y);
    mesh_->addDataVector("Theta Z", theta_z);

    // вычисление напряжений
    std::vector<double> SigmaX(nodesCount, 0.0);
    std::vector<double> SigmaY(nodesCount, 0.0);
    std::vector<double> SigmaZ(nodesCount, 0.0);
    std::vector<double> TauXY(nodesCount, 0.0);
    std::vector<double> TauXZ(nodesCount, 0.0);
    std::vector<double> TauYZ(nodesCount, 0.0);
    std::vector<double> mises(nodesCount, 0.0);

    double xi[elementNodes];
    double eta[elementNodes];
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
    {
        xi[0] = 0.0; eta[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) !=NULL)
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

        ElementPointer element = mesh_->element(elNum);
        Point3D A = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(0))));
        Point3D B = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(1))));
        Point3D C = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(2))));
        DoubleMatrix lambda = cosinuses(A, B, C);
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);

        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        DoubleMatrix lambdaT = lambda.transpose();

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        // извлечение координат узлов
        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
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
            if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
            {
                isoTriangle3(xi[inode], eta[inode], x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
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
                dis(freedom_ * i) = displacement[freedom_ * element->vertexNode(i)];
                dis(freedom_ * i + 1) = displacement[freedom_ * element->vertexNode(i) + 1UL];
                dis(freedom_ * i + 2) = displacement[freedom_ * element->vertexNode(i) + 2UL];
                dis(freedom_ * i + 3) = displacement[freedom_ * element->vertexNode(i) + 3UL];
                dis(freedom_ * i + 4) = displacement[freedom_ * element->vertexNode(i) + 4UL];
                dis(freedom_ * i + 5) = displacement[freedom_ * element->vertexNode(i) + 5UL];
            }

            DoubleVector dLocal = T * dis;

            sigma_membrane = (D_ * Bm) * dLocal;
            sigma_plate = thickness_ / 2.0 * ((D_ * Bf) * dLocal);
            sigma[0] = sigma_membrane[0] + sigma_plate[0];
            sigma[1] = sigma_membrane[1] + sigma_plate[1];
            sigma[2] = sigma_membrane[2] + sigma_plate[2];
            tau = (Dc_ * Bc) * dLocal;

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
        SigmaX[i] /= (double)mesh_->adjacentCount(i);
        SigmaY[i] /= (double)mesh_->adjacentCount(i);
        SigmaZ[i] /= (double)mesh_->adjacentCount(i);
        TauXY[i] /= (double)mesh_->adjacentCount(i);
        TauXZ[i] /= (double)mesh_->adjacentCount(i);
        TauYZ[i] /= (double)mesh_->adjacentCount(i);
        mises[i] /= (double)mesh_->adjacentCount(i);
    }

    mesh_->addDataVector("Sigma X", SigmaX);
    mesh_->addDataVector("Sigma Y", SigmaY);
    mesh_->addDataVector("Sigma Z", SigmaZ);
    mesh_->addDataVector("Tau XY", TauXY);
    mesh_->addDataVector("Tau XZ", TauXZ);
    mesh_->addDataVector("Tau YZ", TauYZ);
    mesh_->addDataVector("von Mises", mises);
}
/*
MindlinShellBending::MindlinShellBending(Mesh3D *mesh,
                                         double thickness,
                                         const std::vector<double> &strain,
                                         const std::vector<double> &stress,
                                         double nu,
                                         const std::list<FemCondition *> &conditions) :
    Fem2D(mesh, 6, conditions)
{
    const double kappa = 5.0 / 6.0;
    // координаты квадратур Гаусса
    DoubleVector gxi;
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

    std::vector<double>::size_type zones_count = stress.size();
    std::vector<DoubleMatrix> D(zones_count);
    std::vector<DoubleMatrix> Dc(zones_count);

    for (std::vector<double>::size_type i = 0; i < zones_count; i++)
    {
        double E = stress[i] / strain[i];
        DoubleMatrix d = evalPlaneStressMatrix(E, nu);
        DoubleMatrix dc(2, 0.0);
        dc(0, 0) = d(2, 2); dc(1, 1) = d(2, 2);
        D[i] = d;
        Dc[i] = dc;
        std::cout << i << ": " << std::endl;
        std::cout << "D:" << std::endl;
        d.print();
        std::cout << "Dc:" << std::endl;
        dc.print();
    }

    std::cout << "Thickness: " << thickness << std::endl;

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    force_ = evalForces(conditions); // Сила

    bool isBreak = false;

    mesh->clearLayers();
    for (msh::UInteger i = 0; i < elementsCount; i++) mesh->pushLayer(0);

    std::vector<double> xxx(nodesCount, 0.0);
    std::vector<double> yyy(nodesCount, 0.0);
    std::vector<double> zzz(nodesCount, 0.0);
    std::vector<double> theta_x(nodesCount, 0.0);
    std::vector<double> theta_y(nodesCount, 0.0);
    std::vector<double> theta_z(nodesCount, 0.0);
    std::vector<double> SigmaX(nodesCount, 0.0);
    std::vector<double> SigmaY(nodesCount, 0.0);
    std::vector<double> SigmaZ(nodesCount, 0.0);
    std::vector<double> TauXY(nodesCount, 0.0);
    std::vector<double> TauXZ(nodesCount, 0.0);
    std::vector<double> TauYZ(nodesCount, 0.0);
    std::vector<double> mises(nodesCount, 0.0);

    int iteration = 0;
    double force_factor = 0.0;

    while (!isBreak)
    {
        // построение глобальной матрицы жесткости
        std::cout << "Stiffness Matrix...";
        ConsoleProgress progressBar(elementsCount);

        global_.resize(dimension_);

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
            DoubleMatrix eD = D[mesh->layer(elNum)];
            DoubleMatrix eDc = Dc[mesh->layer(elNum)];

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
                    jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
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

                local += jacobian * w * thickness * (Bm.transpose() * eD * Bm);
                local += jacobian * w * thickness*thickness*thickness / 12.0 * (Bf.transpose() * eD * Bf);
                local += jacobian * w * kappa * thickness * (Bc.transpose() * eDc * Bc);
            } // ig

            DoubleMatrix surf = T.transpose() * local * T;
            // Ансамблирование
            assembly(element, surf);
        } //for elNum

        processInitialValues();

        DoubleVector displacement = solveLinearSystem();

        // перемещения
        for (UInteger i = 0; i < nodesCount; i++)
        {
            xxx[i] += displacement[freedom_ * i];
            yyy[i] += displacement[freedom_ * i + 1UL];
            zzz[i] += displacement[freedom_ * i + 2UL];
            theta_x[i] += displacement[freedom_ * i + 3UL];
            theta_y[i] += displacement[freedom_ * i + 4UL];
            theta_z[i] += displacement[freedom_ * i + 5UL];
        }

        std::vector<double> localSigmaX(nodesCount, 0.0);
        std::vector<double> localSigmaY(nodesCount, 0.0);
        std::vector<double> localSigmaZ(nodesCount, 0.0);
        std::vector<double> localTauXY(nodesCount, 0.0);
        std::vector<double> localTauXZ(nodesCount, 0.0);
        std::vector<double> localTauYZ(nodesCount, 0.0);
        std::vector<double> localMises(nodesCount, 0.0);

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
            DoubleMatrix eD = D[mesh->layer(elNum)];
            DoubleMatrix eDc = Dc[mesh->layer(elNum)];

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
                    dis(freedom_ * i) = displacement[freedom_ * element->vertexNode(i)];
                    dis(freedom_ * i + 1) = displacement[freedom_ * element->vertexNode(i) + 1UL];
                    dis(freedom_ * i + 2) = displacement[freedom_ * element->vertexNode(i) + 2UL];
                    dis(freedom_ * i + 3) = displacement[freedom_ * element->vertexNode(i) + 3UL];
                    dis(freedom_ * i + 4) = displacement[freedom_ * element->vertexNode(i) + 4UL];
                    dis(freedom_ * i + 5) = displacement[freedom_ * element->vertexNode(i) + 5UL];
                }

                DoubleVector dLocal = T * dis;

                sigma_membrane = (eD * Bm) * dLocal;
                sigma_plate = thickness / 2.0 * ((eD * Bf) * dLocal);
                sigma[0] = sigma_membrane[0] + sigma_plate[0];
                sigma[1] = sigma_membrane[1] + sigma_plate[1];
                sigma[2] = sigma_membrane[2] + sigma_plate[2];
                tau = (eDc * Bc) * dLocal;

                DoubleMatrix localSigma(3, 3, 0.0);
                localSigma(0, 0) = sigma(0); localSigma(0, 1) = sigma(2); localSigma(0, 2) = tau(0);
                localSigma(1, 0) = sigma(2); localSigma(1, 1) = sigma(1); localSigma(1, 2) = tau(1);
                localSigma(2, 0) = tau(0);   localSigma(2, 1) = tau(1);   localSigma(2, 2) = 0.0;

                DoubleMatrix SG = lambdaT * localSigma * lambda;

                double von = sqrt(0.5) * sqrt((SG(0, 0) - SG(1, 1)) * (SG(0, 0) - SG(1, 1)) +
                                              (SG(1, 1) - SG(2, 2)) * (SG(1, 1) - SG(2, 2)) +
                                              (SG(2, 2) - SG(0, 0)) * (SG(2, 2) - SG(0, 0)) +
                                              6.0 * (SG(0, 1) * SG(0, 1) + SG(1, 2) * SG(1, 2) + SG(0, 2) * SG(0, 2)));
                localSigmaX[element->vertexNode(inode)] += SG(0, 0);
                localSigmaY[element->vertexNode(inode)] += SG(1, 1);
                localSigmaZ[element->vertexNode(inode)] += SG(2, 2);
                localTauXY[element->vertexNode(inode)] += SG(0, 1);
                localTauXZ[element->vertexNode(inode)] += SG(0, 2);
                localTauYZ[element->vertexNode(inode)] += SG(1, 2);
                localMises[element->vertexNode(inode)] += von;
            }
        } //for elNum

        for (UInteger i = 0; i < nodesCount; i++)
        {
            localSigmaX[i] /= (double)mesh->adjacentCount(i);
            localSigmaY[i] /= (double)mesh->adjacentCount(i);
            localSigmaZ[i] /= (double)mesh->adjacentCount(i);
            localTauXY[i] /= (double)mesh->adjacentCount(i);
            localTauXZ[i] /= (double)mesh->adjacentCount(i);
            localTauYZ[i] /= (double)mesh->adjacentCount(i);
            localMises[i] /= (double)mesh->adjacentCount(i);

            SigmaX[i] += localSigmaX[i];
            SigmaY[i] += localSigmaY[i];
            SigmaZ[i] += localSigmaZ[i];
            TauXY[i] += localTauXY[i];
            TauXZ[i] += localTauXZ[i];
            TauYZ[i] += localTauYZ[i];
            mises[i] += localMises[i];
        }

        double maxMises = 0.0;
        for (UInteger i = 0; i < nodesCount; i++)
        {
            if (maxMises < mises[i]) maxMises = mises[i];

        }
        if (stress[zones_count - 1] <= maxMises)
        {
            isBreak = true;
        }
        std::cout << "MAX MISES: " << maxMises << std::endl;

        if (iteration == 0)
        {
            double factor = stress[0] / maxMises;
            std::cout << "Initial factor: " << factor << std::endl;
            for (UInteger i = 0; i < nodesCount; i++)
            {
                xxx[i] *= factor;
                yyy[i] *= factor;
                zzz[i] *= factor;
                theta_x[i] *= factor;
                theta_y[i] *= factor;
                theta_z[i] *= factor;
                SigmaX[i] *= factor;
                SigmaY[i] *= factor;
                SigmaZ[i] *= factor;
                TauXY[i] *= factor;
                TauXZ[i] *= factor;
                TauYZ[i] *= factor;
                mises[i] *= factor;
            }
            force_factor += factor;
        }
        else force_factor += 1.0;

        for (UInteger elNum = 0; elNum < elementsCount; elNum++)
        {
            ElementPointer element = mesh->element(elNum);
            double elementMises = 0.0;
            for (UInteger i = 0; i < elementNodes; i++)
            {
                if (elementMises < mises[element->vertexNode(i)]) elementMises = mises[element->vertexNode(i)];
            }
            for (std::vector<double>::size_type j = 0; j < zones_count; j++)
            {
                if (stress[j] < elementMises) mesh->setLayer(elNum, j + 1);
            }
        }
        iteration++;
    }
    std::cout << "<=== FROCE FACTOR: " << force_factor << " ===>" << std::endl;
    mesh_->addDataVector("X", xxx);
    mesh_->addDataVector("Y", yyy);
    mesh_->addDataVector("Z", zzz);
    mesh_->addDataVector("Theta X", theta_x);
    mesh_->addDataVector("Theta Y", theta_y);
    mesh_->addDataVector("Theta Z", theta_z);
    mesh_->addDataVector("Sigma X", SigmaX);
    mesh_->addDataVector("Sigma Y", SigmaY);
    mesh_->addDataVector("Sigma Z", SigmaZ);
    mesh_->addDataVector("Tau XY", TauXY);
    mesh_->addDataVector("Tau XZ", TauXZ);
    mesh_->addDataVector("Tau YZ", TauYZ);
    mesh_->addDataVector("von Mises", mises);
}*/

DoubleMatrix MindlinShellBending::cosinuses(const Point3D &A, const Point3D &B, const Point3D &C)
{
    DoubleMatrix lambda(3, 3, 0.0);
    Point3D AB = B - A;
    Point3D AC = C - A;
    Point3D N = AB.product(AC);
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vz.product(Vx);
    lambda(0, 0) = Vx.x(); lambda(0, 1) = Vx.y(); lambda(0, 2) = Vx.z();
    lambda(1, 0) = Vy.x(); lambda(1, 1) = Vy.y(); lambda(1, 2) = Vy.z();
    lambda(2, 0) = Vz.x(); lambda(2, 1) = Vz.y(); lambda(2, 2) = Vz.z();
    return lambda;
}
