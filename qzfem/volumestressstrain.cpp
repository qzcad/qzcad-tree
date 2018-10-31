#include "volumestressstrain.h"

#include <iostream>
#include <math.h>
#include "consoleprogress.h"

#include "hexahedralmesh3d.h"
#include "tetrahedralmesh3d.h"

VolumeStressStrain::VolumeStressStrain(Mesh3D *mesh, const DoubleMatrix &D, const std::list<FemCondition *> &conditions) : Fem3D(mesh, 3, conditions)
{
    D_ = D;
}

void VolumeStressStrain::solve(std::function<double (double, double, double)> func, double delta, int maxiter)
{
    mesh_->clearDataVectors();
    for (int i = 0; i < maxiter; i++)
    {
        DoubleVector char_vec;
        dimension_ = mesh_->nodesCount() * freedom_;
        global_.resize(dimension_);
        force_.resize(dimension_);
        force_.set(0.0);
        std::cout << "Iteration " << i + 1 << " from " << maxiter << std::endl;
        buildGlobalMatrix();
        buildGlobalVector();
        processInitialValues();
        DoubleVector solution = solveLinearSystem();
        char_vec = adaptationVector(solution);
        double d = char_vec.max() - char_vec.min();
        std::list<UInteger> elements;
        AdjacentSet elset;
        for (UInteger elnum = 0; elnum < mesh_->elementsCount(); elnum++)
        {
            ElementPointer element = mesh_->element(elnum);
            bool needSubdivision = false;
            for (int j = 0; j < element->verticesCount(); j++)
            {
                for (int k = j + 1; k < element->verticesCount(); k++)
                    if ((fabs(char_vec[element->vertexNode(j)] - char_vec[element->vertexNode(k)]) / fabs(d)) >= delta)
                        needSubdivision = true;
            }
//            if (needSubdivision) elements.push_back(elnum);
            if (needSubdivision)
                for (int j = 0; j < element->verticesCount(); j++)
                {
                    AdjacentSet a = mesh_->adjacent(element->vertexNode(j));
                    elset.insert(a.begin(), a.end());
                }
        }
        std::copy(elset.begin(), elset.end(), std::back_inserter(elements));
        if (!elements.empty() && i < maxiter - 1)
        {
            std::cout << elements.size() << " elements must be subdivided!" << std::endl;
            if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
                dynamic_cast<TetrahedralMesh3D *>(mesh_)->subdivide(elements, func);
            else if (dynamic_cast<HexahedralMesh3D *>(mesh_) != nullptr)
                dynamic_cast<HexahedralMesh3D*>(mesh_)->subdivide(elements, func);
        }
        else
        {
            processSolution(solution);
            break;
        }
    }
}

void VolumeStressStrain::buildGlobalMatrix()
{
    const UInteger elementsCount = mesh_->elementsCount();
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gmu;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 2; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    UInteger elementNodes = 0; // количество узлов в элементе

    D_.print();

    quadrature(line_count, line_points, line_weights);

    if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
    {
        elementNodes = 4;
        gaussPoints = 5;
        quadrature(gaussPoints, gxi, geta, gmu, gweight);
        std::cout << "Terahedral Mesh" << std::endl;
    }
    if (dynamic_cast<HexahedralMesh3D*>(mesh_) != nullptr)
    {
        gaussPoints = line_count * line_count * line_count;
        gxi.resize(gaussPoints);
        geta.resize(gaussPoints);
        gmu.resize(gaussPoints);
        gweight.resize(gaussPoints);
        elementNodes = 8;
        for (int i = 0; i < line_count; i++)
        {
            for (int j = 0; j < line_count; j++)
            {
                for (int k = 0; k < line_count; k++)
                {
                    gxi[i * line_count * line_count + j * line_count + k] = line_points[i];
                    geta[i * line_count * line_count + j * line_count + k] = line_points[j];
                    gmu[i * line_count * line_count + j * line_count + k] = line_points[k];
                    gweight[i * line_count * line_count + j * line_count + k] = line_weights[i] * line_weights[j] * line_weights[k];
                }
            }
        }
        std::cout << "Hexahedral Mesh" << std::endl;
    }

    std::cout << "Stiffness Matrix...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        DoubleVector z(elementNodes);
        // извлечение координат узлов
        ElementPointer element = mesh_->element(elementNumber);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh_->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }
        // интегрирование квадратурами
        for (int ig = 0; ig < gaussPoints; ig++)
        {
            double xi = gxi(ig);
            double eta = geta(ig);
            double mu = gmu(ig);
            double w = gweight(ig);
            // значения функций формы
            DoubleVector N(elementNodes);
            // значения производных функций формы в глобальных координатах
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            DoubleVector dNdZ(elementNodes);
            double jacobian = 1.0;
            if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
            {
                jacobian = isoTet4(xi, eta, mu, x, y, z, N, dNdX, dNdY, dNdZ);
            }
            if (dynamic_cast<HexahedralMesh3D *>(mesh_) != nullptr)
            {
                jacobian = isoHex8(xi, eta, mu, x, y, z, N, dNdX, dNdY, dNdZ);
            }
                    // матрицы вариационной постановки
            DoubleMatrix B(6, elementNodes * freedom_, 0.0);
            for (int i = 0; i < elementNodes; i++)
            {
                B(0, i * freedom_) = dNdX(i);
                                              B(1, i * freedom_ + 1) = dNdY(i);
                                                                                B(2, i * freedom_ + 2) = dNdZ(i);
                B(3, i * freedom_) = dNdY(i); B(3, i * freedom_ + 1) = dNdX(i);
                                              B(4, i * freedom_ + 1) = dNdZ(i); B(4, i * freedom_ + 2) = dNdY(i);
                B(5, i * freedom_) = dNdZ(i);                                   B(5, i * freedom_ + 2) = dNdX(i);
            }

            local += jacobian * w * (B.transpose() * D_ * B);
        } // ig
        // Ансамблирование
        // Ансамблирование
        assembly(element, local);
    } // for elementNumber
}

void VolumeStressStrain::buildGlobalVector()
{
    // Учет сил
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    DoubleVector gxi; // координаты квадратур Гаусса
    DoubleVector geta;
    DoubleVector gmu;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 0; // количество точек квадратур
    int line_count = 2; // количество точек квадратур при интегрировании вдоль линии
    DoubleVector line_points;
    DoubleVector line_weights;
    quadrature(line_count, line_points, line_weights);

    UInteger elementNodes = 0; // количество узлов в элементе

    if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
    {
        elementNodes = 4;
        gaussPoints = 5;
        quadrature(gaussPoints, gxi, geta, gmu, gweight);
    }
    if (dynamic_cast<HexahedralMesh3D*>(mesh_) != nullptr)
    {
        gaussPoints = line_count * line_count * line_count;
        gxi.resize(gaussPoints);
        geta.resize(gaussPoints);
        gmu.resize(gaussPoints);
        gweight.resize(gaussPoints);
        elementNodes = 8;
        for (int i = 0; i < line_count; i++)
        {
            for (int j = 0; j < line_count; j++)
            {
                for (int k = 0; k < line_count; k++)
                {
                    gxi[i * line_count * line_count + j * line_count + k] = line_points[i];
                    geta[i * line_count * line_count + j * line_count + k] = line_points[j];
                    gmu[i * line_count * line_count + j * line_count + k] = line_points[k];
                    gweight[i * line_count * line_count + j * line_count + k] = line_weights[i] * line_weights[j] * line_weights[k];
                }
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
                    int dir = (*condition)->direction();
                    if (dir & FemCondition::FIRST)
                        force_(freedom_ * i) += f;
                    if (dir & FemCondition::SECOND)
                        force_(freedom_ * i + 1) += f;
                    if (dir & FemCondition::THIRD)
                        force_(freedom_ * i + 2) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Surface Forces...";
            ConsoleProgress progressBar(elementsCount);
            int dir = (*condition)->direction();
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;
                ElementPointer element = mesh_->element(elNum);
                for (int j = 0; j < element->facesCount(); j++)
                {
                    UIntegerVector face = element->face(j);\
                    bool is_border_face = mesh_->isBorderFace(face);
                    for (UInteger ee = 0; ee < elementsCount && is_border_face; ee++)
                        if (ee != elNum)
                        {
                            ElementPointer inner = mesh_->element(ee);
                            unsigned int k = 0;
                            for (unsigned int ii = 0; ii < face.size(); ii++)
                                if (inner->in(face[ii])) k++;
                            if (k == face.size()) is_border_face = false;
                        }
                    if (is_border_face)
                    {
                        double faceArea = mesh_->faceArea(face); // площадь грани
                        Point3D center (0.0, 0.0, 0.0); // центр грани
                        bool is_applied = true;
                        double val = 0.0;
                        for (int k = 0; k < face.size(); k++)
                        {
                            center = center + dynamic_cast<Mesh3D *>(mesh_)->point3d(face[k]);
                            if (!(*condition)->isApplied(mesh_->node(face[k])))
                                is_applied = false;
                            else
                                val += (*condition)->value(mesh_->node(face[k]));
                        }
                        center.scale(1.0 / static_cast<double>(face.size()));
                        val /= static_cast<double>(face.size());
                        if (/*(*condition)->isApplied(&center)*/is_applied)
                        {
                            double f = /*(*condition)->value(&center)*/val * faceArea / static_cast<double>(face.size());
                            for (int k = 0; k < face.size(); k++)
                            {
                                if (dir & FemCondition::FIRST)
                                {
                                    force_(freedom_ * face[k]) += f;
                                }
                                if (dir & FemCondition::SECOND)
                                {
                                    force_(freedom_ * face[k] + 1) += f;
                                }
                                if (dir & FemCondition::THIRD)
                                {
                                    force_(freedom_ * face[k] + 2) += f;
                                }
                            }
                        }
                    }// if
                } // for j
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
                DoubleVector z(elementNodes);
                double vForce[elementNodes]; // значения объемных сил в узлах
                // извлечение координат узлов
                ElementPointer element = mesh_->element(elNum);
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    PointPointer point = mesh_->node(element->vertexNode(i));
                    x[i] = point->x();
                    y[i] = point->y();
                    z[i] = point->z();
                    vForce[i] = 0.0;
                }
                for (int ig = 0; ig < gaussPoints; ig++)
                {
                    double xi = gxi(ig);
                    double eta = geta(ig);
                    double mu = gmu(ig);
                    double w = gweight(ig);
                    // значения функций формы
                    DoubleVector N(elementNodes);
                    // значения производных функций формы
                    DoubleVector dNdX(elementNodes);
                    DoubleVector dNdY(elementNodes);
                    DoubleVector dNdZ(elementNodes);
                    // якобиан
                    double jacobian = 1.0;
                    if (dynamic_cast<TetrahedralMesh3D*>(mesh_) != nullptr)
                    {
                        jacobian = isoTet4(xi, eta, mu, x, y, z, N, dNdX, dNdY, dNdZ);
                    }
                    if (dynamic_cast<HexahedralMesh3D*>(mesh_) != nullptr)
                    {
                        jacobian = isoHex8(xi, eta, mu, x, y, z, N, dNdX, dNdY, dNdZ);
                    }
                    // вычисление объемных сил
                    Point3D pLocal(x * N, y * N, z * N);
                    double fLocal = (*condition)->value(&pLocal);
                    for (UInteger i = 0; i < elementNodes; i++)
                    {
                        vForce[i] = vForce[i] + (N[i] * jacobian * w) * fLocal;
                    }
                } // ig
                int dir = (*condition)->direction();
                // ансамбль объемных сил
                for (UInteger i = 0 ; i < elementNodes; i++)
                {
                    if (dir & FemCondition::FIRST)
                        force_(freedom_ * element->vertexNode(i)) += vForce[i];
                    if (dir & FemCondition::SECOND)
                        force_(freedom_ * element->vertexNode(i) + 1) += vForce[i];
                    if (dir & FemCondition::THIRD)
                        force_(freedom_ * element->vertexNode(i) + 2) += vForce[i];
                }
            } //for elNum
        }
    } // iterator
}

void VolumeStressStrain::processSolution(const DoubleVector &displacement)
{
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TetrahedralMesh3D*>(mesh_) != nullptr)
    {
        elementNodes = 4;
    }
    if (dynamic_cast<HexahedralMesh3D *>(mesh_) != nullptr)
    {
        elementNodes = 8;
    }

    std::vector<double> u(nodesCount);
    std::vector<double> v(nodesCount);
    std::vector<double> w(nodesCount);
//    std::vector<double> fx(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        u[i] = displacement[freedom_ * i];
        v[i] = displacement[freedom_ * i + 1];
        w[i] = displacement[freedom_ * i + 2];
//        fx[i] = force_[freedom_ * i];
    }

    mesh_->addDataVector("X", u);
    mesh_->addDataVector("Y", v);
    mesh_->addDataVector("Z", w);

    std::vector<double> SigmaX(nodesCount, 0.0);
    std::vector<double> SigmaY(nodesCount, 0.0);
    std::vector<double> SigmaZ(nodesCount, 0.0);
    std::vector<double> TauXY(nodesCount, 0.0);
    std::vector<double> TauXZ(nodesCount, 0.0);
    std::vector<double> TauYZ(nodesCount, 0.0);
    std::vector<double> mises(nodesCount, 0.0);

    double xi[elementNodes];
    double eta[elementNodes];
    double mu[elementNodes];
    if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
    {
        xi[0] = 0.0; eta[0] = 0.0; mu[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0; mu[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0; mu[2] = 0.0;
        xi[3] = 0.0; eta[3] = 0.0; mu[3] = 1.0;
    }
    else if (dynamic_cast<HexahedralMesh3D*>(mesh_) !=nullptr)
    {
        xi[0] = -1.0; eta[0] = -1.0; mu[0] = -1.0;
        xi[1] = -1.0; eta[1] = -1.0; mu[1] =  1.0;
        xi[2] =  1.0; eta[2] = -1.0; mu[2] =  1.0;
        xi[3] =  1.0; eta[3] = -1.0; mu[3] = -1.0;
        xi[4] = -1.0; eta[4] =  1.0; mu[4] = -1.0;
        xi[5] = -1.0; eta[5] =  1.0; mu[5] =  1.0;
        xi[6] =  1.0; eta[6] =  1.0; mu[6] =  1.0;
        xi[7] =  1.0; eta[7] =  1.0; mu[7] = -1.0;
    }

    std::cout << "Stresses Recovery...";
    ConsoleProgress progressBar(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        DoubleVector z(elementNodes);
        // извлечение координат узлов
        ElementPointer element = mesh_->element(elNum);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh_->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }

        for (UInteger inode = 0; inode < elementNodes; inode++)
        {
            DoubleVector N(elementNodes);
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            DoubleVector dNdZ(elementNodes);
            DoubleMatrix dis((size_type)(freedom_ * elementNodes), (size_type)1);
            DoubleMatrix sigma((size_type)6, (size_type)1);
            if (dynamic_cast<TetrahedralMesh3D*>(mesh_) != nullptr)
            {
                isoTet4(xi[inode], eta[inode], mu[inode], x, y, z, N, dNdX, dNdY, dNdZ);
            }
            else if (dynamic_cast<HexahedralMesh3D*>(mesh_) != nullptr)
            {
                isoHex8(xi[inode], eta[inode], mu[inode], x, y, z, N, dNdX, dNdY, dNdZ);
            }
            //
            DoubleMatrix B(6, elementNodes * freedom_, 0.0);
            for (int i = 0; i < elementNodes; i++)
            {
                B(0, i * freedom_) = dNdX(i);
                                              B(1, i * freedom_ + 1) = dNdY(i);
                                                                                B(2, i * freedom_ + 2) = dNdZ(i);
                B(3, i * freedom_) = dNdY(i); B(3, i * freedom_ + 1) = dNdX(i);
                                              B(4, i * freedom_ + 1) = dNdZ(i); B(4, i * freedom_ + 2) = dNdY(i);
                B(5, i * freedom_) = dNdZ(i);                                   B(5, i * freedom_ + 2) = dNdX(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(i * freedom_, 0) = displacement[freedom_ * element->vertexNode(i)];
                dis(i * freedom_ + 1, 0) = displacement[freedom_ * element->vertexNode(i) + 1];
                dis(i * freedom_ + 2, 0) = displacement[freedom_ * element->vertexNode(i) + 2];
            }

            sigma = (D_ * B) * dis;
            SigmaX[element->vertexNode(inode)] += sigma(0, 0);
            SigmaY[element->vertexNode(inode)] += sigma(1, 0);
            SigmaZ[element->vertexNode(inode)] += sigma(2, 0);
            TauXY[element->vertexNode(inode)] += sigma(3, 0);
            TauYZ[element->vertexNode(inode)] += sigma(4, 0);
            TauXZ[element->vertexNode(inode)] += sigma(5, 0);
//            mises[element->vertexNode(inode)] += sqrt(sigma(0,0)*sigma(0,0) - sigma(0,0)*sigma(1,0) + sigma(1,0)*sigma(1,0) + 3.0 * sigma(2,0)*sigma(2,0)); // general plane stress
            mises[element->vertexNode(inode)] += sqrt(0.5 * ((sigma(0, 0) - sigma(1, 0)) * (sigma(0, 0) - sigma(1, 0)) + (sigma(1, 0) - sigma(2, 0)) * (sigma(1, 0) - sigma(2, 0)) + (sigma(2, 0) - sigma(0, 0)) * (sigma(2, 0) - sigma(0, 0)))
                                                      + 3.0 * (sigma(3, 0) * sigma(3, 0) + sigma(4, 0) * sigma(4, 0) + sigma(5, 0) * sigma(5, 0)));
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
    //    mesh_->addDataVector("Fx", fx);
}

DoubleVector VolumeStressStrain::adaptationVector(const DoubleVector &displacement)
{
    UInteger nodesCount = mesh_->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh_->elementsCount(); // количество элементов
    UInteger elementNodes = 0; // количество узлов в элементе
    if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
    {
        elementNodes = 4;
    }
    if (dynamic_cast<HexahedralMesh3D *>(mesh_) != nullptr)
    {
        elementNodes = 8;
    }

    DoubleVector u(nodesCount);
    DoubleVector v(nodesCount);
    DoubleVector w(nodesCount);
    DoubleVector vec(nodesCount);
//    std::vector<double> fx(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        u[i] = displacement[freedom_ * i];
        v[i] = displacement[freedom_ * i + 1];
        w[i] = displacement[freedom_ * i + 2];
        vec[i] = sqrt(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);
//        fx[i] = force_[freedom_ * i];
    }

    DoubleVector SigmaX(nodesCount, 0.0);
    DoubleVector SigmaY(nodesCount, 0.0);
    DoubleVector SigmaZ(nodesCount, 0.0);
    DoubleVector TauXY(nodesCount, 0.0);
    DoubleVector TauXZ(nodesCount, 0.0);
    DoubleVector TauYZ(nodesCount, 0.0);
    DoubleVector mises(nodesCount, 0.0);

    double xi[elementNodes];
    double eta[elementNodes];
    double mu[elementNodes];
    if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
    {
        xi[0] = 0.0; eta[0] = 0.0; mu[0] = 0.0;
        xi[1] = 1.0; eta[1] = 0.0; mu[1] = 0.0;
        xi[2] = 0.0; eta[2] = 1.0; mu[2] = 0.0;
        xi[3] = 0.0; eta[3] = 0.0; mu[3] = 1.0;
    }
    else if (dynamic_cast<HexahedralMesh3D*>(mesh_) !=nullptr)
    {
        xi[0] = -1.0; eta[0] = -1.0; mu[0] = -1.0;
        xi[1] = -1.0; eta[1] = -1.0; mu[1] =  1.0;
        xi[2] =  1.0; eta[2] = -1.0; mu[2] =  1.0;
        xi[3] =  1.0; eta[3] = -1.0; mu[3] = -1.0;
        xi[4] = -1.0; eta[4] =  1.0; mu[4] = -1.0;
        xi[5] = -1.0; eta[5] =  1.0; mu[5] =  1.0;
        xi[6] =  1.0; eta[6] =  1.0; mu[6] =  1.0;
        xi[7] =  1.0; eta[7] =  1.0; mu[7] = -1.0;
    }

    std::cout << "Stresses Recovery...";
    ConsoleProgress progressBar(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

        DoubleVector x(elementNodes);
        DoubleVector y(elementNodes);
        DoubleVector z(elementNodes);
        // извлечение координат узлов
        ElementPointer element = mesh_->element(elNum);
        for (UInteger i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh_->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }

        for (UInteger inode = 0; inode < elementNodes; inode++)
        {
            DoubleVector N(elementNodes);
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            DoubleVector dNdZ(elementNodes);
            DoubleMatrix dis((size_type)(freedom_ * elementNodes), (size_type)1);
            DoubleMatrix sigma((size_type)6, (size_type)1);
            if (dynamic_cast<TetrahedralMesh3D *>(mesh_) != nullptr)
            {
                isoTet4(xi[inode], eta[inode], mu[inode], x, y, z, N, dNdX, dNdY, dNdZ);
            }
            else if (dynamic_cast<HexahedralMesh3D *>(mesh_) != nullptr)
            {
                isoHex8(xi[inode], eta[inode], mu[inode], x, y, z, N, dNdX, dNdY, dNdZ);
            }
            //
            DoubleMatrix B(6, elementNodes * freedom_, 0.0);
            for (int i = 0; i < elementNodes; i++)
            {
                B(0, i * freedom_) = dNdX(i);
                                              B(1, i * freedom_ + 1) = dNdY(i);
                                                                                B(2, i * freedom_ + 2) = dNdZ(i);
                B(3, i * freedom_) = dNdY(i); B(3, i * freedom_ + 1) = dNdX(i);
                                              B(4, i * freedom_ + 1) = dNdZ(i); B(4, i * freedom_ + 2) = dNdY(i);
                B(5, i * freedom_) = dNdZ(i);                                   B(5, i * freedom_ + 2) = dNdX(i);
            }

            for (UInteger i = 0; i < elementNodes; i++)
            {
                dis(i * freedom_, 0) = displacement[freedom_ * element->vertexNode(i)];
                dis(i * freedom_ + 1, 0) = displacement[freedom_ * element->vertexNode(i) + 1];
                dis(i * freedom_ + 2, 0) = displacement[freedom_ * element->vertexNode(i) + 2];
            }

            sigma = (D_ * B) * dis;
//            SigmaX[element->vertexNode(inode)] += sigma(0, 0);
//            SigmaY[element->vertexNode(inode)] += sigma(1, 0);
//            SigmaZ[element->vertexNode(inode)] += sigma(2, 0);
//            TauXY[element->vertexNode(inode)] += sigma(3, 0);
//            TauYZ[element->vertexNode(inode)] += sigma(4, 0);
//            TauXZ[element->vertexNode(inode)] += sigma(5, 0);
//            mises[element->vertexNode(inode)] += sqrt(sigma(0,0)*sigma(0,0) - sigma(0,0)*sigma(1,0) + sigma(1,0)*sigma(1,0) + 3.0 * sigma(2,0)*sigma(2,0)); // general plane stress
            mises[element->vertexNode(inode)] += sqrt(0.5 * ((sigma(0, 0) - sigma(1, 0)) * (sigma(0, 0) - sigma(1, 0)) + (sigma(1, 0) - sigma(2, 0)) * (sigma(1, 0) - sigma(2, 0)) + (sigma(2, 0) - sigma(0, 0)) * (sigma(2, 0) - sigma(0, 0)))
                                                      + 3.0 * (sigma(3, 0) * sigma(3, 0) + sigma(4, 0) * sigma(4, 0) + sigma(5, 0) * sigma(5, 0)));
        }
    } //for elNum
    for (UInteger i = 0; i < nodesCount; i++)
    {
//        SigmaX[i] /= (double)mesh_->adjacentCount(i);
//        SigmaY[i] /= (double)mesh_->adjacentCount(i);
//        SigmaZ[i] /= (double)mesh_->adjacentCount(i);
//        TauXY[i] /= (double)mesh_->adjacentCount(i);
//        TauXZ[i] /= (double)mesh_->adjacentCount(i);
//        TauYZ[i] /= (double)mesh_->adjacentCount(i);
        mises[i] /= static_cast<double>(mesh_->adjacentCount(i));
    }

    return mises;
}

DoubleMatrix VolumeStressStrain::evalStressMatrix(const double &E, const double &nu)
{
    DoubleMatrix S(6, 6, 0.0);

    S(0, 0) = 1.0 / E;    S(0, 1) = - nu / E;    S(0, 2) = - nu / E;
    S(1, 0) = - nu / E;    S(1, 1) = 1.0 / E;    S(1, 2) = - nu / E;
    S(2, 0) = - nu / E;    S(2, 1) = - nu / E;    S(2, 2) = 1.0 / E;
    double G = E / (2.0 * (1.0 + nu));
    S(3, 3) = 1.0 / G;
    S(4, 4) = 1.0 / G;
    S(5, 5) = 1.0 / G;

    return S.inverted();
}
