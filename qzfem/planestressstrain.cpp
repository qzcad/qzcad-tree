#include "planestressstrain.h"

#include <iostream>

#include "rowdoublematrix.h"
#include "mappeddoublematrix.h"

#include "consoleprogress.h"

PlaneStressStrain::PlaneStressStrain(QuadrilateralMesh2D *mesh,
                                     double thickness,
                                     const ElasticMatrix &elasticMatrix,
                                     BoundaryConditionFunction fixFunc,
                                     VectorFunction2D boundaryValue,
                                     VectorFunction2D nodalForce,
                                     VectorFunction2D surfaceForce,
                                     VectorFunction2D volumeForce) : Fem2D(mesh)
{
    freedom_ = 2; // количество степеней свободы
    DoubleVector gpoint; // координаты квадратур Гаусса
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 2; // количество точек квадратуры

    unsigned int elementNodes = 4; // количество узлов в элементе

    DoubleMatrix D = elasticMatrix.D(); // матрица упругости

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    UInteger dimension = freedom_ * nodesCount; // размер системы

    MappedDoubleMatrix global (dimension); // глобальная матрица жесткости
    DoubleVector force(dimension, 0.0); // вектор сил

    quadrature(gaussPoints, gpoint, gweight); // генерация квадратур
    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix & Volume Forces (матрица жесткости и объемные силы)...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        double x[elementNodes];
        double y[elementNodes];
        Point2D vForce[elementNodes]; // значения объемных сил в узлах
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
                // значения функций формы
                DoubleVector N(elementNodes);
                // значения производных функций формы
                DoubleVector dNdX(elementNodes);
                DoubleVector dNdY(elementNodes);
                // якобиан
                double jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                //
                DoubleMatrix B(3, 8, 0.0);
                B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2); B(0, 3) = dNdX(3);
                B(1, 4) = dNdY(0); B(1, 5) = dNdY(1); B(1, 6) = dNdY(2); B(1, 7) = dNdY(3);
                B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2); B(2, 3) = dNdY(3);
                B(2, 4) = dNdX(0); B(2, 5) = dNdX(1); B(2, 6) = dNdX(2); B(2, 7) = dNdX(3);
                local += jacobian * wi * wj * thickness * (B.transpose() * D * B);
                // вычисление объемных сил
                double xLocal = x[0] * N[0] + x[1] * N[1] + x[2] * N[2] + x[3] * N[3];
                double yLocal = y[0] * N[0] + y[1] * N[1] + y[2] * N[2] + y[3] * N[3];
                Point2D fLocal = volumeForce(xLocal, yLocal);
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    vForce[i] = vForce[i] + (N[i] * jacobian * wi * wj) * fLocal;
                }
            } // jeta
        } // ieta
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (unsigned int i = 0; i < elementNodes * freedom_; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else //if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }

            for (unsigned int j = i; j < elementNodes * freedom_; j++)
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
        // ансамбль объемных сил
        for (UInteger i = 0 ; i < elementNodes; i++)
        {
            force(element->vertexNode(i)) += vForce[i].x();
            force(element->vertexNode(i) + nodesCount) += vForce[i].y();
        }
    } //for elNum

    // узловые нагрузки
    std::cout << "Nodal Forces (узловые нагрузки)...";
    progressBar.restart(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh->node(i);
        Point2D f = nodalForce(point->x(), point->y());
        force(i) += f.x();
        force(i + nodesCount) += f.y();

        ++progressBar;
    } // for i

    // поверхностные нагрузки
    std::cout << "Surface Forces (поверхностные нагрузки)...";
    progressBar.restart(elementsCount);
    double length = 0.0;
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {
        ++progressBar;

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
                    length += l;
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
    std::cout << "length = " << length << std::endl;
    // учет граничных условий
    std::cout << "Boundary Conditions (граничные условия)...";
    progressBar.restart(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh->node(i);
        int fixType = fixFunc(point->x(), point->y());
        Point2D fixVal = boundaryValue(point->x(), point->y());
        if (fixType == 0 || fixType == 1)
        {
            setInitialNodalValue(global, force, i, fixVal.x());
        }
        if (fixType == 0 || fixType == 2)
        {
            setInitialNodalValue(global, force, i + nodesCount, fixVal.y());
        }

        ++progressBar;
    } // for i
    // решение СЛАУ
    std::cout << "Linear Equations (СЛАУ)..." << std::endl;
    RowDoubleMatrix rdm(global);
    nodeValues_ = rdm.conjugateGradient(force);

    // вычисление напряжений
    elementVectorsCount_ = 3;
    elementValues_.resize(elementVectorsCount_ * elementsCount);
    std::cout << "Stresses (напряжения)...";
    progressBar.restart(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

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

        double xi = 0.0;
        double eta = 0.0;

        DoubleVector N(elementNodes);
        DoubleVector dNdX(elementNodes);
        DoubleVector dNdY(elementNodes);
        DoubleMatrix displacement((size_type)(freedom_ * elementNodes), (size_type)1);
        DoubleMatrix sigma((size_type)3, (size_type)1);
        // функции формы
        isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
        //
        DoubleMatrix B(3, 8, 0.0);
        B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2); B(0, 3) = dNdX(3);
        B(1, 4) = dNdY(0); B(1, 5) = dNdY(1); B(1, 6) = dNdY(2); B(1, 7) = dNdY(3);
        B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2); B(2, 3) = dNdY(3);
        B(2, 4) = dNdX(0); B(2, 5) = dNdX(1); B(2, 6) = dNdX(2); B(2, 7) = dNdX(3);

        for (unsigned int i = 0; i < elementNodes; i++)
        {
            displacement(i, 0) = nodeValues_[element->vertexNode(i)];
            displacement(i + elementNodes, 0) = nodeValues_[element->vertexNode(i) + nodesCount];
        }

        sigma = (D * B) * displacement;

        elementValues_[elNum] = sigma(0, 0);
        elementValues_[elNum + elementsCount] = sigma(1, 0);
        elementValues_[elNum + 2UL * elementsCount] = sigma(2, 0);
    } //for elNum
}

PlaneStressStrain::PlaneStressStrain(QuadrilateralMesh2D *mesh,
                                     double thickness,
                                     const ElasticMatrix &elasticMatrix,
                                     std::list<FemCondition *> conditions) :
    Fem2D(mesh)
{
    freedom_ = 2; // количество степеней свободы
    DoubleVector gpoint; // координаты квадратур Гаусса
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 2; // количество точек квадратуры

    unsigned int elementNodes = 4; // количество узлов в элементе

    DoubleMatrix D = elasticMatrix.D(); // матрица упругости

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    UInteger dimension = freedom_ * nodesCount; // размер системы

    MappedDoubleMatrix global (dimension); // глобальная матрица жесткости
    DoubleVector force(dimension, 0.0); // вектор сил

    quadrature(gaussPoints, gpoint, gweight); // генерация квадратур
    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix & Volume Forces (матрица жесткости и объемные силы)...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
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
                // значения функций формы
                DoubleVector N(elementNodes);
                // значения производных функций формы
                DoubleVector dNdX(elementNodes);
                DoubleVector dNdY(elementNodes);
                // якобиан
                double jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
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
        for (unsigned int i = 0; i < elementNodes * freedom_; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else //if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }

            for (unsigned int j = i; j < elementNodes * freedom_; j++)
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

    // Учет сил
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::NODAL_FORCE)
        {
            // узловые нагрузки
            std::cout << "Nodal Forces (узловые нагрузки)...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    double f = (*condition)->value(point);
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(i) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(i + nodesCount) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Surface Forces (поверхностные нагрузки)...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;

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
                            FemCondition::FemDirection dir = (*condition)->direction();
                            if ((*condition)->isApplied(point0) && (*condition)->isApplied(point1))
                            {
                                Point2D p0(point0->x(), point0->y());
                                Point2D p1(point1->x(), point1->y());
                                double l = p0.distanceTo(p1);
                                double jacobian = l / 2.0;
                                double f0 = 0.0;
                                double f1 = 0.0;
                                for (int ixi = 0; ixi < gaussPoints; ixi++)
                                {

                                    double xi = gpoint(ixi);
                                    double w = gweight(ixi);
                                    double N0 = (1.0 - xi) / 2.0;
                                    double N1 = (1.0 + xi) / 2.0;
                                    Point2D p = Point2D(p0.x() * N0 + p1.x() * N1,
                                                             p0.y() * N0 + p1.y() * N1);
                                    f0 += N0 * jacobian * w * (*condition)->value(&p);
                                    f1 += N1 * jacobian * w * (*condition)->value(&p);
                                } // for ixi
                                if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                                {
                                    force(element->vertexNode(i)) += f0;
                                    force(element->vertexNode(i + 1)) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                                {
                                    force(element->vertexNode(i) + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount) += f1;
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
            std::cout << "Volume Forces (объемные силы)...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                double x[elementNodes];
                double y[elementNodes];
                double vForce[elementNodes]; // значения объемных сил в узлах
                // извлечение координат узлов
                ElementPointer element = mesh->element(elNum);
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    PointPointer point = mesh->node(element->vertexNode(i));
                    x[i] = point->x();
                    y[i] = point->y();
                    vForce[i] = 0.0;
                }
                for (int ixi = 0; ixi < gaussPoints; ixi++)
                {
                    double xi = gpoint(ixi);
                    double wi = gweight(ixi);
                    for (int jeta = 0; jeta < gaussPoints; jeta++)
                    {
                        double eta = gpoint(jeta);
                        double wj = gweight(jeta);
                        // значения функций формы
                        DoubleVector N(elementNodes);
                        // значения производных функций формы
                        DoubleVector dNdX(elementNodes);
                        DoubleVector dNdY(elementNodes);
                        // якобиан
                        double jacobian = isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
                        // вычисление объемных сил
                        double xLocal = x[0] * N[0] + x[1] * N[1] + x[2] * N[2] + x[3] * N[3];
                        double yLocal = y[0] * N[0] + y[1] * N[1] + y[2] * N[2] + y[3] * N[3];
                        Point2D pLocal(xLocal, yLocal);
                        double fLocal = (*condition)->value(&pLocal);
                        for (unsigned int i = 0; i < elementNodes; i++)
                        {
                            vForce[i] = vForce[i] + (N[i] * jacobian * wi * wj) * fLocal;
                        }
                    } // jeta
                } // ieta
                FemCondition::FemDirection dir = (*condition)->direction();
                // ансамбль объемных сил
                for (UInteger i = 0 ; i < elementNodes; i++)
                {
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(element->vertexNode(i)) += vForce[i];
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(element->vertexNode(i) + nodesCount) += vForce[i];
                }
            } //for elNum
        }
    } // iterator

    //учет условий закрепления
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::INITIAL_VALUE)
        {
            // учет граничных условий
            std::cout << "Boundary Conditions (граничные условия)...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        setInitialNodalValue(global, force, i, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        setInitialNodalValue(global, force, i + nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator

    // решение СЛАУ
    std::cout << "Linear Equations (СЛАУ)..." << std::endl;
    RowDoubleMatrix rdm(global);
    nodeValues_ = rdm.conjugateGradient(force);

    // вычисление напряжений
    elementVectorsCount_ = 3;
    elementValues_.resize(elementVectorsCount_ * elementsCount);
    std::cout << "Stresses (напряжения)...";
    progressBar.restart(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

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

        double xi = 0.0;
        double eta = 0.0;

        DoubleVector N(elementNodes);
        DoubleVector dNdX(elementNodes);
        DoubleVector dNdY(elementNodes);
        DoubleMatrix displacement((size_type)(freedom_ * elementNodes), (size_type)1);
        DoubleMatrix sigma((size_type)3, (size_type)1);
        // функции формы
        isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
        //
        DoubleMatrix B(3, 8, 0.0);
        B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2); B(0, 3) = dNdX(3);
        B(1, 4) = dNdY(0); B(1, 5) = dNdY(1); B(1, 6) = dNdY(2); B(1, 7) = dNdY(3);
        B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2); B(2, 3) = dNdY(3);
        B(2, 4) = dNdX(0); B(2, 5) = dNdX(1); B(2, 6) = dNdX(2); B(2, 7) = dNdX(3);

        for (unsigned int i = 0; i < elementNodes; i++)
        {
            displacement(i, 0) = nodeValues_[element->vertexNode(i)];
            displacement(i + elementNodes, 0) = nodeValues_[element->vertexNode(i) + nodesCount];
        }

        sigma = (D * B) * displacement;

        elementValues_[elNum] = sigma(0, 0);
        elementValues_[elNum + elementsCount] = sigma(1, 0);
        elementValues_[elNum + 2UL * elementsCount] = sigma(2, 0);
    } //for elNum
}

PlaneStressStrain::PlaneStressStrain(TriangleMesh2D *mesh,
                                     double thickness,
                                     const ElasticMatrix &elasticMatrix,
                                     BoundaryConditionFunction fixFunc,
                                     VectorFunction2D boundaryValue,
                                     VectorFunction2D nodalForce,
                                     VectorFunction2D surfaceForce,
                                     VectorFunction2D volumeForce) : Fem2D(mesh)
{
    freedom_ = 2; // количество степеней свободы
    // координаты квадратур Гаусса
    DoubleVector gxi;
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 3; // количество точек квадратуры

    unsigned int elementNodes = 3; // количество узлов в элементе

    DoubleMatrix D = elasticMatrix.D(); // матрица упругости

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    UInteger dimension = freedom_ * nodesCount; // размер системы

    MappedDoubleMatrix global (dimension); // глобальная матрица жесткости
    DoubleVector force(dimension, 0.0); // вектор сил

    quadrature(gaussPoints, gxi, geta, gweight); // генерация квадратур
    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix & Volume Forces (матрица жесткости и объемные силы)...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
        double x[elementNodes];
        double y[elementNodes];
        Point2D vForce[elementNodes]; // значения объемных сил в узлах
        // извлечение координат узлов
        ElementPointer element = mesh->element(elNum);
        for (unsigned int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh->node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
        }
        for (int ig = 0; ig < gaussPoints; ig++)
        {
            double xi = gxi(ig);
            double eta = geta(ig);
            double weight = gweight(ig);
            // значения функций формы
            DoubleVector N(elementNodes);
            // значения производных функций формы
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            // якобиан
            double jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            //
            DoubleMatrix B(3, 6, 0.0);
            B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2);
            B(1, 3) = dNdY(0); B(1, 4) = dNdY(1); B(1, 5) = dNdY(2);
            B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2);
            B(2, 3) = dNdX(0); B(2, 4) = dNdX(1); B(2, 5) = dNdX(2);
            local += 0.5 * jacobian * weight * thickness * (B.transpose() * D * B);
            // вычисление объемных сил
            double xLocal = x[0] * N[0] + x[1] * N[1] + x[2] * N[2];
            double yLocal = y[0] * N[0] + y[1] * N[1] + y[2] * N[2];
            Point2D fLocal = volumeForce(xLocal, yLocal);
            for (unsigned int i = 0; i < elementNodes; i++)
            {
                vForce[i] = vForce[i] + (N[i] * 0.5 * jacobian * weight) * fLocal;
            }
        } // ieta
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (unsigned int i = 0; i < elementNodes * freedom_; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else //if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }

            for (unsigned int j = i; j < elementNodes * freedom_; j++)
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
        // ансамбль объемных сил
        for (UInteger i = 0 ; i < elementNodes; i++)
        {
            force(element->vertexNode(i)) += vForce[i].x();
            force(element->vertexNode(i) + nodesCount) += vForce[i].y();
        }
    } //for elNum

    // узловые нагрузки
    std::cout << "Nodal Forces (узловые нагрузки)...";
    progressBar.restart(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh->node(i);
        Point2D f = nodalForce(point->x(), point->y());
        force(i) += f.x();
        force(i + nodesCount) += f.y();

        ++progressBar;
    } // for i

    // поверхностные нагрузки
    std::cout << "Surface Forces (поверхностные нагрузки)...";
    progressBar.restart(elementsCount);
    double length = 0.0;
    int line_points = 2;
    DoubleVector line_xi;
    DoubleVector line_weight;
    quadrature(line_points, line_xi, line_weight); // генерация квадратур
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {
        ++progressBar;

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
                    length += l;
                    Point2D f0(0.0, 0.0);
                    Point2D f1(0.0, 0.0);
                    for (int ixi = 0; ixi < line_points; ixi++)
                    {

                        double xi = line_xi(ixi);
                        double w = line_weight(ixi);
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
    std::cout << "length = " << length << std::endl;
    // учет граничных условий
    std::cout << "Boundary Conditions (граничные условия)...";
    progressBar.restart(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh->node(i);
        int fixType = fixFunc(point->x(), point->y());
        Point2D fixVal = boundaryValue(point->x(), point->y());
        if (fixType == 0 || fixType == 1)
        {
            setInitialNodalValue(global, force, i, fixVal.x());
        }
        if (fixType == 0 || fixType == 2)
        {
            setInitialNodalValue(global, force, i + nodesCount, fixVal.y());
        }

        ++progressBar;
    } // for i
    // решение СЛАУ
    std::cout << "Linear Equations (СЛАУ)..." << std::endl;
    RowDoubleMatrix rdm(global);
    nodeValues_ = rdm.conjugateGradient(force);

    // вычисление напряжений
    elementVectorsCount_ = 3;
    elementValues_.resize(elementVectorsCount_ * elementsCount);
    std::cout << "Stresses (напряжения)...";
    progressBar.restart(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

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
        // центр элемента
        double xi = 1.0 / 3.0;
        double eta = 1.0 / 3.0;

        DoubleVector N(elementNodes);
        DoubleVector dNdX(elementNodes);
        DoubleVector dNdY(elementNodes);
        DoubleMatrix displacement((size_type)(freedom_ * elementNodes), (size_type)1);
        DoubleMatrix sigma((size_type)3, (size_type)1);
        // функции формы
        isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
        //
        DoubleMatrix B(3, 6, 0.0);
        B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2);
        B(1, 3) = dNdY(0); B(1, 4) = dNdY(1); B(1, 5) = dNdY(2);
        B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2);
        B(2, 3) = dNdX(0); B(2, 4) = dNdX(1); B(2, 5) = dNdX(2);

        for (unsigned int i = 0; i < elementNodes; i++)
        {
            displacement(i, 0) = nodeValues_[element->vertexNode(i)];
            displacement(i + elementNodes, 0) = nodeValues_[element->vertexNode(i) + nodesCount];
        }

        sigma = (D * B) * displacement;

        elementValues_[elNum] = sigma(0, 0);
        elementValues_[elNum + elementsCount] = sigma(1, 0);
        elementValues_[elNum + 2UL * elementsCount] = sigma(2, 0);
    } //for elNum
}

PlaneStressStrain::PlaneStressStrain(TriangleMesh2D *mesh,
                                     double thickness,
                                     const ElasticMatrix &elasticMatrix,
                                     std::list<FemCondition *> conditions) :
    Fem2D(mesh)
{
    freedom_ = 2; // количество степеней свободы
    // координаты квадратур Гаусса
    DoubleVector gxi;
    DoubleVector geta;
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 3; // количество точек квадратуры

    unsigned int elementNodes = 3; // количество узлов в элементе

    DoubleMatrix D = elasticMatrix.D(); // матрица упругости

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    UInteger dimension = freedom_ * nodesCount; // размер системы

    MappedDoubleMatrix global (dimension); // глобальная матрица жесткости
    DoubleVector force(dimension, 0.0); // вектор сил

    quadrature(gaussPoints, gxi, geta, gweight); // генерация квадратур
    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix & Volume Forces (матрица жесткости и объемные силы)...";
    ConsoleProgress progressBar(elementsCount);

    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;
        DoubleMatrix local(freedom_ * elementNodes, freedom_ * elementNodes, 0.0);
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
        for (int ig = 0; ig < gaussPoints; ig++)
        {
            double xi = gxi(ig);
            double eta = geta(ig);
            double weight = gweight(ig);
            // значения функций формы
            DoubleVector N(elementNodes);
            // значения производных функций формы
            DoubleVector dNdX(elementNodes);
            DoubleVector dNdY(elementNodes);
            // якобиан
            double jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
            //
            DoubleMatrix B(3, 6, 0.0);
            B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2);
            B(1, 3) = dNdY(0); B(1, 4) = dNdY(1); B(1, 5) = dNdY(2);
            B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2);
            B(2, 3) = dNdX(0); B(2, 4) = dNdX(1); B(2, 5) = dNdX(2);
            local += 0.5 * jacobian * weight * thickness * (B.transpose() * D * B);
        } // ieta
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (unsigned int i = 0; i < elementNodes * freedom_; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else //if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }

            for (unsigned int j = i; j < elementNodes * freedom_; j++)
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

    // Учет сил
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::INITIAL_VALUE)
        {
            // учет граничных условий
            std::cout << "Boundary Conditions (граничные условия)...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        setInitialNodalValue(global, force, i, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        setInitialNodalValue(global, force, i + nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
        else if ((*condition)->type() == FemCondition::NODAL_FORCE)
        {
            // узловые нагрузки
            std::cout << "Nodal Forces (узловые нагрузки)...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    double f = (*condition)->value(point);
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(i) += f;
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(i + nodesCount) += f;
                }
                ++progressBar;
            } // for i

        }
        else if ((*condition)->type() == FemCondition::SURFACE_FORCE)
        {
            // поверхностные нагрузки
            std::cout << "Surface Forces (поверхностные нагрузки)...";
            int line_points = 2;
            DoubleVector  line_point;
            DoubleVector line_weight;
            quadrature(line_points,line_point, line_weight);
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {
                ++progressBar;

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
                            FemCondition::FemDirection dir = (*condition)->direction();
                            if ((*condition)->isApplied(point0) && (*condition)->isApplied(point1))
                            {
                                Point2D p0(point0->x(), point0->y());
                                Point2D p1(point1->x(), point1->y());
                                double l = p0.distanceTo(p1);
                                double jacobian = l / 2.0;
                                double f0 = 0.0;
                                double f1 = 0.0;
                                for (int ixi = 0; ixi < line_points; ixi++)
                                {

                                    double xi = line_point(ixi);
                                    double w = line_weight(ixi);
                                    double N0 = (1.0 - xi) / 2.0;
                                    double N1 = (1.0 + xi) / 2.0;
                                    Point2D p = Point2D(p0.x() * N0 + p1.x() * N1,
                                                             p0.y() * N0 + p1.y() * N1);
                                    f0 += N0 * jacobian * w * (*condition)->value(&p);
                                    f1 += N1 * jacobian * w * (*condition)->value(&p);
                                } // for ixi
                                if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                                {
                                    force(element->vertexNode(i)) += f0;
                                    force(element->vertexNode(i + 1)) += f1;
                                }
                                if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                                {
                                    force(element->vertexNode(i) + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount) += f1;
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
            std::cout << "Volume Forces (объемные силы)...";
            progressBar.restart(elementsCount);
            for (UInteger elNum = 0; elNum < elementsCount; elNum++)
            {

                ++progressBar;
                double x[elementNodes];
                double y[elementNodes];
                double vForce[elementNodes]; // значения объемных сил в узлах
                // извлечение координат узлов
                ElementPointer element = mesh->element(elNum);
                for (unsigned int i = 0; i < elementNodes; i++)
                {
                    PointPointer point = mesh->node(element->vertexNode(i));
                    x[i] = point->x();
                    y[i] = point->y();
                    vForce[i] = 0.0;
                }
                for (int ig = 0; ig < gaussPoints; ig++)
                {
                    double xi = gxi(ig);
                    double eta = geta(ig);
                    double weight = gweight(ig);
                    // значения функций формы
                    DoubleVector N(elementNodes);
                    // значения производных функций формы
                    DoubleVector dNdX(elementNodes);
                    DoubleVector dNdY(elementNodes);
                    // якобиан
                    double jacobian = isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
                    // вычисление объемных сил
                    double xLocal = x[0] * N[0] + x[1] * N[1] + x[2] * N[2];
                    double yLocal = y[0] * N[0] + y[1] * N[1] + y[2] * N[2];
                    Point2D pLocal(xLocal, yLocal);
                    double fLocal = (*condition)->value(&pLocal);
                    for (unsigned int i = 0; i < elementNodes; i++)
                    {
                        vForce[i] = vForce[i] + (N[i] * jacobian * 0.5 * weight) * fLocal;
                    }
                } // ig
                FemCondition::FemDirection dir = (*condition)->direction();
                // ансамбль объемных сил
                for (UInteger i = 0 ; i < elementNodes; i++)
                {
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        force(element->vertexNode(i)) += vForce[i];
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        force(element->vertexNode(i) + nodesCount) += vForce[i];
                }
            } //for elNum
        }
    } // iterator

    // Учет условий закрепления
    for (std::list<FemCondition *>::iterator condition = conditions.begin(); condition != conditions.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::INITIAL_VALUE)
        {
            // учет граничных условий
            std::cout << "Boundary Conditions (граничные условия)...";
            progressBar.restart(nodesCount);
            for (UInteger i = 0; i < nodesCount; i++)
            {
                PointPointer point = mesh->node(i);
                if ((*condition)->isApplied(point))
                {
                    FemCondition::FemDirection dir = (*condition)->direction();
                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        setInitialNodalValue(global, force, i, (*condition)->value(point));
                    if (dir == FemCondition::ALL || dir == FemCondition::SECOND)
                        setInitialNodalValue(global, force, i + nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator

    // решение СЛАУ
    std::cout << "Linear Equations (СЛАУ)..." << std::endl;
    RowDoubleMatrix rdm(global);
    nodeValues_ = rdm.conjugateGradient(force);

    // вычисление напряжений
    elementVectorsCount_ = 3;
    elementValues_.resize(elementVectorsCount_ * elementsCount);
    std::cout << "Stresses (напряжения)...";
    progressBar.restart(elementsCount);
    for (UInteger elNum = 0; elNum < elementsCount; elNum++)
    {

        ++progressBar;

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
        // центр элемента
        double xi = 1.0 / 3.0;
        double eta = 1.0 / 3.0;

        DoubleVector N(elementNodes);
        DoubleVector dNdX(elementNodes);
        DoubleVector dNdY(elementNodes);
        DoubleMatrix displacement((size_type)(freedom_ * elementNodes), (size_type)1);
        DoubleMatrix sigma((size_type)3, (size_type)1);
        // функции формы
        isoTriangle3(xi, eta, x, y, N, dNdX, dNdY);
        //
        DoubleMatrix B(3, 6, 0.0);
        B(0, 0) = dNdX(0); B(0, 1) = dNdX(1); B(0, 2) = dNdX(2);
        B(1, 3) = dNdY(0); B(1, 4) = dNdY(1); B(1, 5) = dNdY(2);
        B(2, 0) = dNdY(0); B(2, 1) = dNdY(1); B(2, 2) = dNdY(2);
        B(2, 3) = dNdX(0); B(2, 4) = dNdX(1); B(2, 5) = dNdX(2);

        for (unsigned int i = 0; i < elementNodes; i++)
        {
            displacement(i, 0) = nodeValues_[element->vertexNode(i)];
            displacement(i + elementNodes, 0) = nodeValues_[element->vertexNode(i) + nodesCount];
        }

        sigma = (D * B) * displacement;

        elementValues_[elNum] = sigma(0, 0);
        elementValues_[elNum + elementsCount] = sigma(1, 0);
        elementValues_[elNum + 2UL * elementsCount] = sigma(2, 0);
    } //for elNum
}

std::string PlaneStressStrain::nodeVectorName(UInteger num) const
{
    switch (num) {
    case 0:
        return "X";
    case 1:
        return "Y";
    default:
        return "undefined";
    }
}

std::string PlaneStressStrain::elementVectorName(UInteger num) const
{
    switch (num) {
    case 0:
        return "Sigma X";
    case 1:
        return "Sigma Y";
    case 2:
        return "Tau XY";
    default:
        return "undefined";
    }
}

