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

    x_.resize(nodesCount, 0.0);
    y_.resize(nodesCount, 0.0);
    z_.resize(nodesCount, 0.0);
    theta_x_.resize(nodesCount, 0.0);
    theta_y_.resize(nodesCount, 0.0);
    theta_z_.resize(nodesCount, 0.0);
    sigma_x_.resize(nodesCount, 0.0);
    sigma_y_.resize(nodesCount, 0.0);
    sigma_z_.resize(nodesCount, 0.0);
    tau_xy_.resize(nodesCount, 0.0);
    tau_xz_.resize(nodesCount, 0.0);
    tau_yz_.resize(nodesCount, 0.0);
    mises_.resize(nodesCount, 0.0);

    buildGlobalMatrix();
    buildGlobalVector();
    processInitialValues();
    DoubleVector solution = solveLinearSystem();
    processSolution(solution);
    f_ = 1.0;
    double f_max[] = {6.26,	32.34,	64.28,	96.22,	130.82,	165.6,	196.24};
    double m = mises_.max();
    double E0 = E_(m);
    int idx = 0;
    std::cout << "Plastic analysis" << std::endl;
    std::cout << "Initial Elastic Modulus: " << E0 << std::endl;
    std::cout << "Initial von Mises Stress: " << m << "(linear fracture load factor is " << sigma_max_ / m << ")" << std::endl;
    while (fabs(E0 - E_(f_ * m)) < 1.0E-7 && f_ < f_max[0])
    {
        f_ += 1.0;
    }
    x_.scale(f_);
    y_.scale(f_);
    z_.scale(f_);
    theta_x_.scale(f_);
    theta_y_.scale(f_);
    theta_z_.scale(f_);
    sigma_x_.scale(f_);
    sigma_y_.scale(f_);
    sigma_z_.scale(f_);
    tau_xy_.scale(f_);
    tau_xz_.scale(f_);
    tau_yz_.scale(f_);
    mises_.scale(f_);
    std::cout << "Initial Force Factor: " << f_ << std::endl;
    std::cout << "Yield Elastic Modulus: " << E_(f_ * m) << std::endl;
    std::cout << "Current Sigma: " << f_ * m << std::endl;
    do
    {
        if (idx < 7 && (f_ > f_max[idx]))
        {
            double factor = f_max[idx] / f_;
            std::cout << "f_max: " << f_max[idx] << std::endl;
            std::cout << "Factor: " << factor << std::endl;
            double polus = 1000.0, shpan = 1000.0;
            UInteger ip=0, is=0;
            for (UInteger i = 0; i < mesh_->nodesCount(); i++)
            {
                PointPointer p = mesh_->node(i);
                if (fabs(0.11 + 0.0435 + 1.037 - p->y()) < polus)
                {
                    polus = fabs(0.11 + 0.0435 + 1.037 - p->y());
                    ip = i;
                }
                if (fabs(0.11 + 0.0435 - p->y()) < shpan)
                {
                    shpan = fabs(0.11 + 0.0435 - p->y());
                    is = i;
                }
            }
            std::cout << "=================================================================================" << endl;
            std::cout << "Polus distance: " << polus << " Polus displacement: " << y_[ip] * factor << endl;
            std::cout << "Shpangout distance: " << shpan << " Shpangout displacement: " << y_[is] * factor << endl;
            std::cout << "=================================================================================" << endl;
            idx++;
        }
        if (idx >= 7)
        {
            break;
        }
        alpha_ = 0; // !!!!
        buildGlobalMatrix();
        force_.set(0.0);
        buildGlobalVector();
        processInitialValues();
        DoubleVector solution = solveLinearSystem();
        processSolution(solution);
        f_ += 1.0;
        m = mises_.max();
        std::cout << "Current Force Factor: " << f_ << std::endl;
        std::cout << "Current Sigma: " << m << ". Elastic Modulus: " << E_(m) << std::endl;
        if(m >= sigma_max_)
        {
            std::cout << "DAMAGE FORCE: " << f_ << std::endl;
        }
    } while(1);// while (m < sigma_max_);

    mesh_->clearDataVectors();
    mesh_->addDataVector("X", x_.to_std());
    mesh_->addDataVector("Y", y_.to_std());
    mesh_->addDataVector("Z", z_.to_std());
    mesh_->addDataVector("Theta X", theta_x_.to_std());
    mesh_->addDataVector("Theta Y", theta_y_.to_std());
    mesh_->addDataVector("Theta Z", theta_z_.to_std());
    mesh_->addDataVector("Sigma X", sigma_x_.to_std());
    mesh_->addDataVector("Sigma Y", sigma_y_.to_std());
    mesh_->addDataVector("Sigma Z", sigma_z_.to_std());
    mesh_->addDataVector("Tau XY", tau_xy_.to_std());
    mesh_->addDataVector("Tau XZ", tau_xz_.to_std());
    mesh_->addDataVector("Tau YZ", tau_yz_.to_std());
    mesh_->addDataVector("mises", mises_.to_std());
    std::cout << "Mises stress: " << m << ". Force factor: " << f_ << std::endl;
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
    double min_e = E_(0.0), max_e = min_e;

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
        DoubleMatrix W(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        DoubleVector epsilonForce(freedom_ * elementNodes, 0.0);
        for (UInteger i = 0; i <= (freedom_ * elementNodes - 3); i += 3)
        {
            for (UInteger ii = 0; ii < 3; ii++)
                for (UInteger jj = 0; jj < 3; jj++)
                    T(ii + i, jj + i) = lambda(ii, jj);
        }

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        DoubleVector z(elementNodes);
        DoubleVector nodal_thickness(elementNodes);
        double thickness = 0.0;
        double sigma=0.0;
        double elastic_modulus;
        // извлечение координат узлов

        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
            z[i] = lambda(2, 0) * pp.x() + lambda(2, 1) * pp.y() + lambda(2, 2) * pp.z();

            if (thickness_func_ != nullptr)
                nodal_thickness[i] = thickness_func_(point.x(), point.y(), point.z());
            else
                nodal_thickness[i] = thickness_;

            if (sigma < mises_[element->vertexNode(i)])
                sigma = mises_[element->vertexNode(i)];
        }

        for (UInteger i = 0; i < elementNodes; i++)
        {
            UInteger offset = i * freedom_;
            for (UInteger ii = 0; ii < freedom_; ii++)
                    W(ii + offset, ii + offset) = 1.0;

            W(3 + offset, 1 + offset) = z[i];
            W(4 + offset, 0 + offset) = -z[i];
        }

        elastic_modulus = E_(sigma);
        D_ = evalPlaneStressMatrix(elastic_modulus, nu_);
        Dc_(0, 0) = D_(2, 2); Dc_(0, 1) = 0.0;
        Dc_(1, 0) = 0.0; Dc_(1, 1) = D_(2, 2);

        if (elastic_modulus > max_e)
            max_e = elastic_modulus;

        if (elastic_modulus < min_e)
            min_e = elastic_modulus;

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
//            if (thickness_func_ != nullptr)
//            {
//                double xLocal = x * N;
//                double yLocal = y * N;
//                double zLocal = z * N;
//                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal + lambdaT(0, 2) * zLocal,
//                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal + lambdaT(1, 2) * zLocal,
//                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal + lambdaT(2, 2) * zLocal);
//                Point3D pLocal = pl + A;
//                thickness_ = thickness_func_(pLocal.x(), pLocal.y(), pLocal.z());
//            }

            thickness = nodal_thickness * N;

            local += jacobian * w * thickness * (Bm.transpose() * D_ * Bm);
            local += jacobian * w * thickness*thickness*thickness / 12.0 * (Bf.transpose() * D_ * Bf);
            local += jacobian * w * kappa * thickness * (Bc.transpose() * Dc_ * Bc);

            epsilonForce += jacobian * w * thickness * ((Bm.transpose() * D_) * epsilon0);
        } // ig
//        It is noted that when an isoparametric membrane element is used as a part of
//        the flat shell element, there is no stiffness associated with the rotation degree
//        of freedom θz . This lack of stiffness lead to singularity in the global stiffness
//        matrix when all the elements are coplanar. A simple method for remedying this
//        singularity is to insert a small fictitious stiffness to each drilling degree of freedom
//        (Zienkiewicz and Taylor, 2000). This is done by simply replacing the null values of
//        the stiffness corresponding to the drilling degree of freedom by a value of 1/1000
//        of the largest diagonal term of the element stiffness matrix.
//        My observation: max_local * 1.0E+3 is the best
        double max_local = local(0, 0);
        for (UInteger i = 0; i < freedom_ * elementNodes; i++)
        {
            if (max_local < local(i, i)) max_local = local(i, i);
        }
        for (UInteger i = 0; i < elementNodes; i++)
        {
            if (fabs(local(i * freedom_ + 5, i * freedom_ + 5)) < 1.0E-3)
                local(i * freedom_ + 5, i * freedom_ + 5) = max_local * 1.0;
        }

        DoubleMatrix surf = T.transpose() * (W * local * W.transpose()) * T;
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
    std::cout << min_e << " <= Elastic Modulus <= " << max_e << std::endl;
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
    DoubleVector xxx(nodesCount);
    DoubleVector yyy(nodesCount);
    DoubleVector zzz(nodesCount);
    DoubleVector theta_x(nodesCount);
    DoubleVector theta_y(nodesCount);
    DoubleVector theta_z(nodesCount);
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
    DoubleVector SigmaX(nodesCount, 0.0);
    DoubleVector SigmaY(nodesCount, 0.0);
    DoubleVector SigmaZ(nodesCount, 0.0);
    DoubleVector TauXY(nodesCount, 0.0);
    DoubleVector TauXZ(nodesCount, 0.0);
    DoubleVector TauYZ(nodesCount, 0.0);
    DoubleVector mises(nodesCount, 0.0);

    ConsoleProgress progressBar(elementsCount);
    double min_e = E_(0.0), max_e = min_e;

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
        DoubleVector z(elementNodes);
        DoubleVector nodal_thickness(elementNodes);
        double thickness = 0.0;
        double s = 0.0;
        double elastic_modulus;
        // извлечение координат узлов
        for (UInteger i = 0; i < elementNodes; i++)
        {
            Point3D point = *(dynamic_cast<const Point3D *>(mesh_->node(element->vertexNode(i))));
            Point3D pp = point - A;
            x[i] = lambda(0, 0) * pp.x() + lambda(0, 1) * pp.y() + lambda(0, 2) * pp.z();
            y[i] = lambda(1, 0) * pp.x() + lambda(1, 1) * pp.y() + lambda(1, 2) * pp.z();
            z[i] = lambda(2, 0) * pp.x() + lambda(2, 1) * pp.y() + lambda(2, 2) * pp.z();

            if (thickness_func_ != nullptr)
                nodal_thickness[i] = thickness_func_(point.x(), point.y(), point.z());
            else
                nodal_thickness[i] = thickness_;

            if (s < mises_[element->vertexNode(i)])
                s = mises_[element->vertexNode(i)];
        }

        elastic_modulus = E_(s);
        D_ = evalPlaneStressMatrix(elastic_modulus, nu_);
        Dc_(0, 0) = D_(2, 2); Dc_(0, 1) = 0.0;
        Dc_(1, 0) = 0.0; Dc_(1, 1) = D_(2, 2);

        if (elastic_modulus > max_e)
            max_e = elastic_modulus;

        if (elastic_modulus < min_e)
            min_e = elastic_modulus;

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

//            if (thickness_func_ != nullptr)
//            {
//                double xLocal = x * N;
//                double yLocal = y * N;
//                double zLocal = z * N;
//                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal + lambdaT(0, 2) * zLocal,
//                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal + lambdaT(1, 2) * zLocal,
//                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal + lambdaT(2, 2) * zLocal);

////                Point3D pl(lambdaT(0, 0) * xLocal + lambdaT(0, 1) * yLocal,
////                           lambdaT(1, 0) * xLocal + lambdaT(1, 1) * yLocal,
////                           lambdaT(2, 0) * xLocal + lambdaT(2, 1) * yLocal);
//                // вычисление объемных сил
//                Point3D pLocal = pl + A;
//                thickness_ = thickness_func_(pLocal.x(), pLocal.y(), pLocal.z());
//            }
            thickness = nodal_thickness * N;

            sigma_membrane = (D_ * Bm) * dLocal;
            sigma_plate = thickness / 2.0 * ((D_ * Bf) * dLocal);
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
    x_ += xxx;
    y_ += yyy;
    z_ += zzz;
    theta_x_ += theta_x;
    theta_y_ += theta_y;
    theta_z_ += theta_z;
    sigma_x_ += SigmaX;
    sigma_y_ += SigmaY;
    sigma_z_ += SigmaZ;
    tau_xy_ += TauXY;
    tau_xz_ += TauXZ;
    tau_yz_ += TauYZ;
    mises_ += mises;
    std::cout << min_e << " <= Elastic Modulus <= " << max_e << std::endl;
}
