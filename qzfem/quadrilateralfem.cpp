#include "quadrilateralfem.h"
#include "gaussquadrature.h"

QuadrilateralFEM::QuadrilateralFEM()
{
}

void QuadrilateralFEM::planeStressStrain(QuadrilateralMesh2D *mesh,
                                         double thickness,
                                         const ElasticMatrix &elasticMatrix,
                                         std::function<int (double, double)> fixFunc,
                                         std::function<double (double, double)> nodalForce,
                                         std::function<double (double, double)> surfaceForce,
                                         std::function<double (double, double)> volumeForce)
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
        for (int i = 0; i < gaussPoints; i++)
        {
            double xi = gpoint(i);
            double wi = gweight(i);
            for (int j = 0; j < gaussPoints; j++)
            {
                double eta = gpoint(j);
                double wj = gweight(j);

                DoubleVector N(elementNodes);
                DoubleVector dNdXi(elementNodes);
                DoubleVector dNdEta(elementNodes);
                DoubleMatrix Jacobi(2, 2); // матрица Якоби
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
                    // координаты
                    x = nodes(quads(elnum, :), 1)
                    y = nodes(quads(elnum, :), 2)
                    // матрица Якоби
                    Jacobi(1, 1) = sum(dNdXi .* x); Jacobi(1, 2) = sum(dNdXi .* y);
                    Jacobi(2, 1) = sum(dNdEta .* x); Jacobi(2, 2) = sum(dNdEta .* y);
                    // якобиан
                    jacobian = det(Jacobi)
                    // обратный якобиан
                    invJacobi = inv(Jacobi)
                    // производные в глобальных координатах
                    dNdX = invJacobi(1, 1) * dNdXi + invJacobi(1, 2) * dNdEta
                    dNdY = invJacobi(2, 1) * dNdXi + invJacobi(2, 2) * dNdEta
            }
        }
    }
}
