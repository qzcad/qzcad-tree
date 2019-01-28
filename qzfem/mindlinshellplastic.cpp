#include "mindlinshellplastic.h"

#include <iostream>
#include <math.h>

#include "consoleprogress.h"

MindlinShellPlastic::MindlinShellPlastic(Mesh3D *mesh, double thickness, std::function<double(double)> E, const double  &sigma_max, const double &nu, const std::list<FemCondition *> &conditions, double alphaT) :
    MindlinShellBending(mesh, thickness, evalPlaneStressMatrix(E(0), nu), conditions, alphaT)
{
    E_ = E;
    nu_ = nu;
    sigma_max_ = sigma_max;
}

MindlinShellPlastic::MindlinShellPlastic(Mesh3D *mesh, std::function<double (double, double, double)> thickness_func, std::function<double (double)> E, const double  &sigma_max, const double &nu, const std::list<FemCondition *> &conditions, double alphaT) :
    MindlinShellBending(mesh, thickness_func, evalPlaneStressMatrix(E(0), nu), conditions, alphaT)
{
    E_ = E;
    nu_ = nu;
    sigma_max_ = sigma_max;
}

void MindlinShellPlastic::solve()
{
    UInteger nodesCount = mesh_->nodesCount();
    mises_.resize(nodesCount);
    mises_.set(0.0);
    buildGlobalMatrix();
    buildGlobalVector();
    processInitialValues();
    DoubleVector solution = solveLinearSystem();
    processSolution(solution);
    double f = 1.0;
    double m = mises_.max();
    double E0 = E_(m);
    std::cout << "Plastic analysis" << std::endl;
    std::cout << "Initial Elastic Modulus: " << E0 << std::endl;
    std::cout << "Initial von Mises Stress: " << m << "(linear fracture load factor is " << sigma_max_ / m << ")" << std::endl;
    while (fabs(E0 - E_(f * m)) < 1.0E-7)
    {
        f += 1.0;
    }
    mises_.scale(f);
    std::cout << "Yield Elastic Modulus: " << E_(f * m) << std::endl;
    std::cout << f << ": Current Sigma: " << f * m << std::endl;
    do
    {
        alpha_ = 0; // !!!!
        buildGlobalMatrix();
        force_.set(0.0);
        buildGlobalVector();
        processInitialValues();
        DoubleVector solution = solveLinearSystem();
        processSolution(solution);
        f += 1.0;
        m = mises_.max();
        std::cout << f << ": Current Sigma: " << m << ". Elastic Modulus: " << E_(m) << std::endl;
    } while (m < sigma_max_);
    std::vector<double> mises(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        mises[i] = mises_[i];
    }
    mesh_->addDataVector("mises", mises);
    std::cout << "Mises stress: " << m << ". Force factor: " << f << std::endl;
}

void MindlinShellPlastic::buildGlobalMatrix()
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
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != nullptr)
    {
        elementNodes = 3;
        gaussPoints = 4;
        quadrature(gaussPoints, gxi, geta, gweight);
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != nullptr)
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
    if (thickness_func_ == nullptr)
        std::cout << "Thickness: " << thickness_ << std::endl;
    else
        std::cout << "Thickness is variable." << std::endl;
    // температурные деформации
    DoubleVector epsilon0(3, 0.0);
    epsilon0(0) = alpha_;
    epsilon0(1) = alpha_;

    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix...";
    ConsoleProgress progressBar(elementsCount);
    UInteger plastic=0;
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
        DoubleVector epsilonForce(freedom_ * elementNodes, 0.0);
        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        double sigma=0.0;
        // извлечение координат узлов

        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
            if (sigma < mises_[element->vertexNode(i)])
                sigma = mises_[element->vertexNode(i)];
        }
        D_ = evalPlaneStressMatrix(E_(sigma), nu_);
        Dc_(0, 0) = D_(2, 2); Dc_(0, 1) = 0.0;
        Dc_(1, 0) = 0.0; Dc_(1, 1) = D_(2, 2);

        if (sigma >= 90.0) plastic++;

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
            if (dynamic_cast<TriangleMesh3D*>(mesh_) != nullptr)
            {
                jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            }
            else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != nullptr)
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
            if (thickness_func_ != nullptr)
            {
                double xLocal = x * N;
                double yLocal = y * N;
                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
                Point3D pLocal = pl + A;
                thickness_ = thickness_func_(pLocal.x(), pLocal.y(), pLocal.z());
            }



            local += jacobian * w * thickness_ * (Bm.transpose() * D_ * Bm);
            local += jacobian * w * thickness_*thickness_*thickness_ / 12.0 * (Bf.transpose() * D_ * Bf);
            local += jacobian * w * kappa * thickness_ * (Bc.transpose() * Dc_ * Bc);

            epsilonForce += jacobian * w * thickness_ * ((Bm.transpose() * D_) * epsilon0);
        } // ig
        double max_local = -1.0E20;
        for (UInteger i = 0; i < freedom_ * elementNodes; i++)
        {
            for (UInteger j = 0; j < freedom_ * elementNodes; j++)
            {
                if (max_local < local(i, j)) max_local = local(i, j);
            }
        }
        for (UInteger i = 0; i < elementNodes; i++)
        {
            local(i * freedom_ + 5, i * freedom_ + 5) = max_local * 1.0E-3;
        }

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
    std::cout << "Number of plastic element: " << plastic << std::endl;
}

void MindlinShellPlastic::processSolution(const DoubleVector &displacement)
{
    UInteger nodesCount = mesh_->nodesCount();
    UInteger elementsCount = mesh_->elementsCount();
    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != nullptr)
    {
        elementNodes = 3;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) != nullptr)
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

    double xi[elementNodes];
    double eta[elementNodes];
    if (dynamic_cast<TriangleMesh3D*>(mesh_) != nullptr)
    {
        xi[0] = 0.0; eta[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0;
    }
    else if (dynamic_cast<QuadrilateralMesh3D*>(mesh_) !=nullptr)
    {
        xi[0] = -1.0; eta[0] = -1.0;
        xi[1] =  1.0; eta[1] = -1.0;
        xi[2] =  1.0; eta[2] =  1.0;
        xi[3] = -1.0; eta[3] =  1.0;
    }

    std::cout << "Stresses Recovery...";
    DoubleVector mises(nodesCount, 0.0);
    UInteger plastic = 0;
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
        double s=0.0;
        // извлечение координат узлов
        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
            if (s < mises_[element->vertexNode(i)])
                s = mises_[element->vertexNode(i)];
        }
        D_ = evalPlaneStressMatrix(E_(s), nu_);
        Dc_(0, 0) = D_(2, 2); Dc_(0, 1) = 0.0;
        Dc_(1, 0) = 0.0; Dc_(1, 1) = D_(2, 2);

        if (s >= 90.0) plastic++;

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

            if (thickness_func_ != nullptr)
            {
                double xLocal = x * N;
                double yLocal = y * N;

                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
                // вычисление объемных сил
                Point3D pLocal = pl + A;
                thickness_ = thickness_func_(pLocal.x(), pLocal.y(), pLocal.z());
            }


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
            mises[element->vertexNode(inode)] += von;
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
        mises[i] /= (double)mesh_->adjacentCount(i);
    }
    mises_ += mises;
    std::cout << "Number of plastic element: " << plastic << std::endl;
}
