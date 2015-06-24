#include "mindlinplatebending.h"

#include "consoleprogress.h"
#include "rowdoublematrix.h"

MindlinPlateBending::MindlinPlateBending(QuadrilateralMesh2D *mesh,
                                         double thickness,
                                         const ElasticMatrix &elasticMatrix,
                                         std::list<FemCondition *> conditions) : Fem2D(mesh)
{
    freedom_ = 3; // количество степеней свободы
    DoubleVector gpoint; // координаты квадратур Гаусса
    DoubleVector gweight; // весовые коэффициенты квадратур
    int gaussPoints = 3; // количество точек квадратуры
    const double kappa = 5.0 / 6.0;

    unsigned int elementNodes = 4; // количество узлов в элементе

    DoubleMatrix D = elasticMatrix.D(); // матрица упругости
    DoubleMatrix Dc(2, 0.0);
    Dc(0, 0) = D(2, 2); Dc(1, 1) = D(2, 2);

    UInteger nodesCount = mesh->nodesCount(); // количество узлов сетки
    UInteger elementsCount = mesh->elementsCount(); // количество элементов

    UInteger dimension = freedom_ * nodesCount; // размер системы

    MappedDoubleMatrix global (dimension); // глобальная матрица жесткости
    DoubleVector force(dimension, 0.0); // вектор сил

    quadrature(gaussPoints, gpoint, gweight); // генерация квадратур
    // построение глобальной матрицы жесткости
    std::cout << "Stiffness Matrix (матрица жесткости)...";
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
                DoubleMatrix Bf(3, 12, 0.0);
                DoubleMatrix Bc(2, 12, 0.0);

                Bf(0, 0) = 0.0; Bf(0, 1) = dNdX(0); Bf(0, 2) = 0.0;     Bf(0, 3) = 0.0; Bf(0, 4) = dNdX(1); Bf(0, 5) = 0.0;     Bf(0, 6) = 0.0; Bf(0, 7) = dNdX(2); Bf(0, 8) = 0.0;     Bf(0, 9) = 0.0; Bf(0, 10) = dNdX(3);    Bf(0, 11) = 0.0;
                Bf(1, 0) = 0.0; Bf(1, 1) = 0.0;     Bf(1, 2) = dNdY(0); Bf(1, 3) = 0.0; Bf(1, 4) = 0.0;     Bf(1, 5) = dNdY(1); Bf(1, 6) = 0.0; Bf(1, 7) = 0.0;     Bf(1, 8) = dNdY(2); Bf(1, 9) = 0.0; Bf(1, 10) = 0.0;        Bf(1, 11) = dNdY(3);
                Bf(2, 0) = 0.0; Bf(2, 1) = dNdY(0); Bf(2, 2) = dNdX(0); Bf(2, 3) = 0.0; Bf(2, 4) = dNdY(1); Bf(2, 5) = dNdX(1); Bf(2, 6) = 0.0; Bf(2, 7) = dNdY(2); Bf(2, 8) = dNdX(2); Bf(2, 9) = 0.0; Bf(2, 10) = dNdY(3);    Bf(2, 11) = dNdX(3);

                Bc(0, 0) = dNdX(0); Bc(0, 1) = N(0);    Bc(0, 2) = 0.0;     Bc(0, 3) = dNdX(1); Bc(0, 4) = N(1);    Bc(0, 5) = 0.0;     Bc(0, 6) = dNdX(2); Bc(0, 7) = N(2);    Bc(0, 8) = 0.0;     Bc(0, 9) = dNdX(3); Bc(0, 10) = N(3);   Bc(0, 11) = 0.0;
                Bc(1, 0) = dNdY(0); Bc(1, 1) = 0.0;     Bc(1, 2) = N(0);    Bc(1, 3) = dNdY(1); Bc(1, 4) = 0.0;     Bc(1, 5) = N(1);    Bc(1, 6) = dNdY(2); Bc(1, 7) = 0.0;     Bc(1, 8) = N(2);    Bc(1, 9) = dNdY(3); Bc(1, 10) = 0.0;    Bc(1, 11) = N(3);
                local += jacobian * wi * wj * thickness*thickness*thickness / 12.0 * (Bf.transpose() * D * Bf);
                local += jacobian * wi * wj * kappa * thickness * (Bc.transpose() * Dc * Bc);
            } // jeta
        } // ieta
        // Ансамблирование
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (unsigned int i = 0; i < elementNodes * freedom_; i++)
        {
            if (i == 0 || i == 3 || i == 6 || i == 9)
            {
                index_i = element->vertexNode(i / 3);
            }
            else if (i == 1 || i == 4 || i == 7 || i == 10)
            {
                index_i = element->vertexNode(i / 3) + nodesCount;
            }
            else //if (i == 2 || i == 5 || i == 8 || i == 11)
            {
                index_i = element->vertexNode(i / 3) + nodesCount + nodesCount;
            }

            for (unsigned int j = i; j < elementNodes * freedom_; j++)
            {
                if (j == 0 || j == 3 || j == 6 || j == 9)
                {
                    index_j = element->vertexNode(j / 3);
                }
                else if (j == 1 || j == 4 || j == 7 || j == 10)
                {
                    index_j = element->vertexNode(j / 3) + nodesCount;
                }
                else //if (j == 2 || j == 5 || ij == 8 || j == 11)
                {
                    index_j = element->vertexNode(j / 3) + nodesCount + nodesCount;
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(i + nodesCount + nodesCount) += f;
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
                                if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                                {
                                    force(element->vertexNode(i) + nodesCount + nodesCount) += f0;
                                    force(element->vertexNode(i + 1) + nodesCount + nodesCount) += f1;
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        force(element->vertexNode(i) + nodesCount + nodesCount) += vForce[i];
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
                    if (dir == FemCondition::ALL || dir == FemCondition::THIRD)
                        setInitialNodalValue(global, force, i + nodesCount + nodesCount, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        }
    } // iterator

    solve(global, force);

    // вычисление напряжений
    elementVectorsCount_ = 5UL;
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
        DoubleMatrix tau((size_type)2, (size_type)1);
        // функции формы
        isoQuad4(xi, eta, x, y, N, dNdX, dNdY);
        //
        DoubleMatrix Bf(3, 12, 0.0);
        DoubleMatrix Bc(2, 12, 0.0);

        Bf(0, 0) = 0.0; Bf(0, 1) = dNdX(0); Bf(0, 2) = 0.0;     Bf(0, 3) = 0.0; Bf(0, 4) = dNdX(1); Bf(0, 5) = 0.0;     Bf(0, 6) = 0.0; Bf(0, 7) = dNdX(2); Bf(0, 8) = 0.0;     Bf(0, 9) = 0.0; Bf(0, 10) = dNdX(3);    Bf(0, 11) = 0.0;
        Bf(1, 0) = 0.0; Bf(1, 1) = 0.0;     Bf(1, 2) = dNdY(0); Bf(1, 3) = 0.0; Bf(1, 4) = 0.0;     Bf(1, 5) = dNdY(1); Bf(1, 6) = 0.0; Bf(1, 7) = 0.0;     Bf(1, 8) = dNdY(2); Bf(1, 9) = 0.0; Bf(1, 10) = 0.0;        Bf(1, 11) = dNdY(3);
        Bf(2, 0) = 0.0; Bf(2, 1) = dNdY(0); Bf(2, 2) = dNdX(0); Bf(2, 3) = 0.0; Bf(2, 4) = dNdY(1); Bf(2, 5) = dNdX(1); Bf(2, 6) = 0.0; Bf(2, 7) = dNdY(2); Bf(2, 8) = dNdX(2); Bf(2, 9) = 0.0; Bf(2, 10) = dNdY(3);    Bf(2, 11) = dNdX(3);

        Bc(0, 0) = dNdX(0); Bc(0, 1) = N(0);    Bc(0, 2) = 0.0;     Bc(0, 3) = dNdX(1); Bc(0, 4) = N(1);    Bc(0, 5) = 0.0;     Bc(0, 6) = dNdX(2); Bc(0, 7) = N(2);    Bc(0, 8) = 0.0;     Bc(0, 9) = dNdX(3); Bc(0, 10) = N(3);   Bc(0, 11) = 0.0;
        Bc(1, 0) = dNdY(0); Bc(1, 1) = 0.0;     Bc(1, 2) = N(0);    Bc(1, 3) = dNdY(1); Bc(1, 4) = 0.0;     Bc(1, 5) = N(1);    Bc(1, 6) = dNdY(2); Bc(1, 7) = 0.0;     Bc(1, 8) = N(2);    Bc(1, 9) = dNdY(3); Bc(1, 10) = 0.0;    Bc(1, 11) = N(3);

        for (unsigned int i = 0; i < elementNodes; i++)
        {
            displacement(3 * i, 0) = nodeValues_[element->vertexNode(i)];
            displacement(3 * i + 1, 0) = nodeValues_[element->vertexNode(i) + nodesCount];
            displacement(3 * i + 2, 0) = nodeValues_[element->vertexNode(i) + nodesCount + nodesCount];
        }

        sigma = (D * Bf) * displacement;
        tau = (Dc * Bc) * displacement;

        elementValues_[elNum] = sigma(0, 0);
        elementValues_[elNum + elementsCount] = sigma(1, 0);
        elementValues_[elNum + 2UL * elementsCount] = sigma(2, 0);
        elementValues_[elNum + 3UL * elementsCount] = tau(0, 0);
        elementValues_[elNum + 4UL * elementsCount] = tau(1, 0);
    } //for elNum
}

std::string MindlinPlateBending::nodeVectorName(UInteger num) const
{
    switch (num) {
    case 0:
        return "W";
    case 1:
        return "Theta X";
    case 2:
        return "Theta Y";
    default:
        return "undefined";
    }
}

std::string MindlinPlateBending::elementVectorName(UInteger num) const
{
    switch (num) {
    case 0:
        return "Sigma X";
    case 1:
        return "Sigma Y";
    case 2:
        return "Tau XY";
    case 3:
        return "Tau XZ";
    case 4:
        return "Tau YZ";
    default:
        return "undefined";
    }
}
