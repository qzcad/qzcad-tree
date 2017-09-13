#include "mindlinshelllaminated.h"
#include <iostream>
#include <math.h>
#include "consoleprogress.h"

MindlinShellLaminated::MindlinShellLaminated(Mesh3D *mesh,
                                             const std::vector<double> &thickness,
                                             const std::vector<DoubleMatrix> &planeStressMatrix,
                                             const std::list<FemCondition *> &conditions,
                                             double alphaT) :
    MindlinShellBending(mesh, thickness[0], planeStressMatrix[0], conditions, alphaT)
{
    thickness_ = thickness;
    D_ = planeStressMatrix;
    thickness_func_ = NULL;
}

MindlinShellLaminated::MindlinShellLaminated(Mesh3D *mesh,
                                             std::function<std::vector<double> (double, double, double)> thickness,
                                             const std::vector<DoubleMatrix> &planeStressMatrix,
                                             const std::list<FemCondition *> &conditions,
                                             double alphaT) :
    MindlinShellBending(mesh, 1.0, planeStressMatrix[0], conditions, alphaT)
{
    D_ = planeStressMatrix;
    thickness_func_ = thickness;
}

void MindlinShellLaminated::buildGlobalMatrix()
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

    unsigned layers_count = D_.size();
    std::vector<DoubleMatrix> D(layers_count); // матрица упругости
    std::vector<DoubleMatrix> Dc(layers_count);
    std::vector<double> thickness = thickness_;
    double H = 0.0; // толщина

    for (unsigned i = 0; i < layers_count; i++)
    {
        D[i] = D_[i];
        Dc[i].resize(2, 2);
        Dc[i] (0, 1) = Dc[i] (1, 0) = 0.0;
        Dc[i] (0, 0) = D[i] (2, 2);     Dc[i] (1, 1) = D[i] (2, 2);
        std::cout << "D[" << i << "]:" << std::endl;
        D[i].print();
        std::cout << "Dc["<< i << "]:" << std::endl;
        Dc[i].print();

        if (!thickness_.empty()) H += thickness_[i];
    }
    if (!thickness_.empty())
        std::cout << "Thickness: " << H << std::endl;
    else
        std::cout << "The laminated shell with a variable thickness." << std::endl;

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
        DoubleMatrix lambdaT = lambda.transpose();
        DoubleMatrix T(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
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
            if (thickness_func_ != NULL)
            {
                double xLocal = x * N;
                double yLocal = y * N;

                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
                // вычисление объемных сил
                Point3D pLocal = pl + A;
                thickness = thickness_func_(pLocal.x(), pLocal.y(), pLocal.z());
                H = 0.0;
                for (double h: thickness)
                    H += h;
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
        assembly(element, surf);
    } //for elNum
}

void MindlinShellLaminated::processSolution(const DoubleVector &displacement)
{
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != NULL)
    {
        elementNodes = 3;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != NULL)
    {
        elementNodes = 4;
    }

    unsigned layers_count = D_.size();
    std::vector<DoubleMatrix> D(layers_count); // матрица упругости
    std::vector<DoubleMatrix> Dc(layers_count);
    double H = 0.0; // толщина
    for (unsigned i = 0; i < layers_count; i++)
    {
        D[i] = D_[i];
        Dc[i].resize(2, 2);
        Dc[i] (0, 1) = Dc[i] (1, 0) = 0.0;
        Dc[i] (0, 0) = D[i] (2, 2);     Dc[i] (1, 1) = D[i] (2, 2);
        std::cout << "D[" << i << "]:" << std::endl;
        D[i].print();
        std::cout << "Dc["<< i << "]:" << std::endl;
        Dc[i].print();

        if (!thickness_.empty()) H += thickness_[i];
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

            sigma_membrane = (D[layers_count - 1] * Bm) * dLocal;
            if (thickness_func_ != NULL)
            {
                double xLocal = x * N;
                double yLocal = y * N;

                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
                // вычисление объемных сил
                Point3D pLocal = pl + A;
                std::vector<double> thickness = thickness_func_(pLocal.x(), pLocal.y(), pLocal.z());
                H = 0.0;
                for (double h: thickness)
                    H += h;
            }
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
