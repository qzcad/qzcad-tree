#include "quadrilateralfem.h"

#include <iostream>

#include "gaussquadrature.h"
#include "doublevector.h"
#include "doublematrix.h"
#include "mappeddoublematrix.h"
#include "rowdoublematrix.h"

using namespace mtx;

QuadrilateralFEM::QuadrilateralFEM()
{
}

void QuadrilateralFEM::planeStressStrain(QuadrilateralMesh2D *mesh,
                                         double thickness,
                                         const ElasticMatrix &elasticMatrix,
                                         std::function<int (double, double)> fixFunc,
                                         std::function<Point2D(double, double)> boundaryValue,
                                         std::function<Point2D(double, double)> nodalForce,
                                         std::function<Point2D(double, double)> surfaceForce,
                                         std::function<Point2D(double, double)> volumeForce)
{
    GaussQuadrature quadrature;
    DoubleVector gpoint;
    DoubleVector gweight;
    int gaussPoints = 2; // количество точек квадратуры

    unsigned int elementNodes = 4;

    DoubleMatrix D = elasticMatrix.D();

    UInteger nodesCount = mesh->nodesCount();
    UInteger elementsCount = mesh->elementsCount();

    UInteger freedom = 2; // количество степеней свободы
    UInteger dimension = freedom * nodesCount; // размер системы

    MappedDoubleMatrix global (dimension);
    DoubleVector force(dimension, 0.0);
    DoubleVector displacement(dimension);

    quadrature.quadrature(gaussPoints, gpoint, gweight);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {
        DoubleMatrix local(freedom * elementNodes, freedom * elementNodes, 0.0);
        double x[elementNodes];
        double y[elementNodes];
        // извлечение координат узлов
        ElementPointer element = mesh->element(elNum);
        for (unsigned int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
        }
        for (int ixi = 0; ixi < gaussPoints; ixi++)
        {
            double xi = gpoint(ixi);
            double wi = gweight(ixi);
            for (int jeta = 0; jeta < gaussPoints; jeta++)
            {
                double eta = gpoint(jeta);
                double wj = gweight(jeta);

                DoubleVector N(elementNodes);
                DoubleVector dNdXi(elementNodes);
                DoubleVector dNdEta(elementNodes);
                DoubleMatrix Jacobi(2, 0.0); // матрица Якоби
                DoubleVector dNdX(elementNodes);
                DoubleVector dNdY(elementNodes);
                // билинейные функции формы
                N(0) = (1.0 - xi) * (1.0 - eta) / 4.0;
                N(1) = (1.0 + xi) * (1.0 - eta) / 4.0;
                N(2) = (1.0 + xi) * (1.0 + eta) / 4.0;
                N(3) = (1.0 - xi) * (1.0 + eta) / 4.0;
                // производные функций формы
                dNdXi(0) = -(1.0 - eta) / 4.0;
                dNdXi(1) = (1.0 - eta) / 4.0;
                dNdXi(2) = (1.0 + eta) / 4.0;
                dNdXi(3) = -(1.0 + eta) / 4.0;

                dNdEta(0) = -(1.0 - xi) / 4.0;
                dNdEta(1) = -(1.0 + xi) / 4.0;
                dNdEta(2) = (1.0 + xi) / 4.0;
                dNdEta(3) = (1.0 - xi) / 4.0;
                // матрица Якоби
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    Jacobi(0, 0) += dNdXi(i) * x[i];  Jacobi(0, 1) += dNdXi(i) * y[i];
                    Jacobi(1, 0) += dNdEta(i) * x[i]; Jacobi(1, 1) += dNdEta(i) * y[i];
                }
                // якобиан
                double jacobian = Jacobi.det2x2();
                // обратный якобиан
                DoubleMatrix invJacobi = Jacobi.inverted2x2();
                // производные в глобальных координатах
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    dNdX(i) = invJacobi(0, 0) * dNdXi(i) + invJacobi(0, 1) * dNdEta(i);
                    dNdY(i) = invJacobi(1, 0) * dNdXi(i) + invJacobi(1, 1) * dNdEta(i);
                }
                //
                DoubleMatrix B(3, 8, 0.0);
                B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2); B(0, 3) = dNdX(3);
                B(1, 4) = dNdY(0); B(1, 5) = dNdY(1); B(1, 6) = dNdY(2); B(1, 7) = dNdY(3);
                B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2); B(2, 3) = dNdY(3);
                B(2, 4) = dNdX(0); B(2, 5) = dNdX(1); B(2, 6) = dNdX(2); B(2, 7) = dNdX(3);
                local += jacobian * wi * wj * thickness * (B.transpose() * D * B);
            } // jeta
        } // ieta
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (unsigned int i = 0; i < elementNodes * freedom; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else //if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }

            for (unsigned int j = i; j < elementNodes * freedom; j++)
            {
                if (j < elementNodes)
                {
                    index_j = element->vertexNode(j);
                }
                else //if (j < 2 * elementNodes)
                {
                    index_j = element->vertexNode(j - elementNodes) + nodesCount;
                }

                global(index_i, index_j) += local(i, j);
                if (index_i != index_j) global(index_j, index_i) = global(index_i, index_j);
            } // for j
        } // for i
    } //for elNum
    // узловые нагрузки
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh->node(i);
        Point2D f = nodalForce(point->x(), point->y());
        force(i) += f.x();
        force(i + nodesCount) = f.y();
    } // for i
    // поверхностные нагрузки
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {
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
                    Point2D p0(point0->x(), point0->y());
                    Point2D p1(point1->x(), point1->y());
                    double l = p0.distanceTo(p1);
                    double jacobian = l / 2.0;
                    Point2D f0(0.0, 0.0);
                    Point2D f1(0.0, 0.0);
                    for (int ixi = 0; ixi < gaussPoints; ixi++)
                    {

                        double xi = gpoint(ixi);
                        double w = gweight(ixi);
                        double N0 = (1.0 - xi) / 2.0;
                        double N1 = (1.0 + xi) / 2.0;
                        Point2D p = surfaceForce(p0.x() * N0 + p1.x() * N1,
                                                 p0.y() * N0 + p1.y() * N1);
                        f0 = f0 + (N0 * jacobian * w * p);
                        f1 = f1 + (N1 * jacobian * w * p);
                    } // for ixi
                    force(element->vertexNode(i)) += f0.x();
                    force(element->vertexNode(i) + nodesCount) += f0.y();
                    force(element->vertexNode(i + 1)) += f1.x();
                    force(element->vertexNode(i + 1) + nodesCount) += f1.y();
                } // if
            } // for i
        } // if
    } // for elNum
    // учет граничных условий
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh->node(i);
        int fixType = fixFunc(point->x(), point->y());
        Point2D fixVal = boundaryValue(point->x(), point->y());
        if (fixType == 0 || fixType == 1)
        {
            for (UInteger j = 0; j < dimension; j++)
            {
                if (i != j && global.data(i, j) != 0)
                { // см. Зенкевич, стр. 485
                    force(j) = force(j) - global.data(i, j) * fixVal.x();
                }
            } // for j
            global.zeroSym(i);
            force(i) = fixVal.x();
            global(i, i) = 1.0;
        }
        if (fixType == 0 || fixType == 2)
        {
            UInteger rowNumber = i + nodesCount;
            for (UInteger j = 0; j < dimension; j++)
            {
                if (rowNumber != j && global.data(rowNumber, j) != 0)
                { // см. Зенкевич, стр. 485
                    force(j) = force(j) - global.data(rowNumber, j) * fixVal.y();
                }
            } // for j
            global.zeroSym(rowNumber);
            force(rowNumber) = fixVal.y();
            global(rowNumber, rowNumber) = 1.0;
        }
    } // for i

    RowDoubleMatrix rdm(global);
    displacement = rdm.conjugateGradient(force);
    u_.clear();
    v_.clear();
    for (UInteger i = 0; i < nodesCount; i++)
    {
        u_.push_back(displacement[i]);
        v_.push_back(displacement[i + nodesCount]);
    }
    double maxU = u_[0];
    double maxV = v_[0];
    double minU = u_[0];
    double minV = v_[0];
    for (UInteger i = 1; i < u_.size(); i++)
    {
        double u = u_[i], v = v_[i];
        // max
        if (maxU < u)
        {
            maxU = u;
        }
        if (maxV < v)
        {
            maxV = v;
        }
        // min
        if (minU > u)
        {
            minU = u;
        }
        if (minV > v)
        {
            minV = v;
        }
    } // for i
    std::cout << "Первое (x) направление:\t" << minU << " <= U <= " << maxU << std::endl;
    std::cout << "Второе (y) направление:\t" << minV << " <= V <= " << maxV << std::endl;
}

std::vector<double> QuadrilateralFEM::u() const
{
    return u_;
}

void QuadrilateralFEM::setU(const std::vector<double> &u)
{
    u_ = u;
}

std::vector<double> QuadrilateralFEM::v() const
{
    return v_;
}

void QuadrilateralFEM::setV(const std::vector<double> &v)
{
    v_ = v;
}

std::vector<double> QuadrilateralFEM::sigmaX() const
{
    return sigmaX_;
}

void QuadrilateralFEM::setSigmaX(const std::vector<double> &sigmaX)
{
    sigmaX_ = sigmaX;
}
std::vector<double> QuadrilateralFEM::sigmaY() const
{
    return sigmaY_;
}

void QuadrilateralFEM::setSigmaY(const std::vector<double> &sigmaY)
{
    sigmaY_ = sigmaY;
}
std::vector<double> QuadrilateralFEM::tauXY() const
{
    return tauXY_;
}

void QuadrilateralFEM::setTauXY(const std::vector<double> &tauXY)
{
    tauXY_ = tauXY;
}





