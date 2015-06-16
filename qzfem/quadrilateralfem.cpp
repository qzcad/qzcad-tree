#include "quadrilateralfem.h"
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

    int elementNodes = 4;

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
        for (int i = 0; i < elementNodes; i++)
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
                for (int i = 0; i < elementNodes; i++)
                {
                    Jacobi(0, 0) += dNdXi(i) * x[i];  Jacobi(0, 1) += dNdXi(i) * y[i];
                    Jacobi(1, 0) += dNdEta(i) * x[i]; Jacobi(1, 1) += dNdEta(i) * y[i];
                }
                // якобиан
                double jacobian = Jacobi.det2x2();
                // обратный якобиан
                DoubleMatrix invJacobi = Jacobi.inverted2x2();
                // производные в глобальных координатах
                for (int i = 0; i < elementNodes; i++)
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
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else //if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }

            for (int j = i; j < elementNodes * freedom; j++)
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
}
