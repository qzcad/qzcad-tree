#include "quadrilateralunion2d.h"
#include <iostream>

namespace msh
{
QuadrilateralUnion2D::QuadrilateralUnion2D()
{
    isFirstRegion_ = true;
    layerNumber_ = 0;
    std::cout << "Создана сетка - объединение сеток четырехугольных элементов..." << std::endl;
}

void QuadrilateralUnion2D::addQuadRegion(const UInteger &xCount, const UInteger &yCount, const Point2D &v0, const Point2D &v1, const Point2D &v2, const Point2D &v3)
{
    double hx = 2.0 / (double)(xCount - 1);
    double hy = 2.0 / (double)(yCount - 1);
    UInteger nodeNumber[xCount * yCount];
    for (UInteger i = 0; i < xCount; i++)
    {
        double xi = -1.0 + (double) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double eta = -1.0 + (double) j * hy;
            Point2D point = isoFunc(0, xi, eta) * v0  + isoFunc(1, xi, eta) * v1 + isoFunc(2, xi, eta) * v2 + isoFunc(3, xi, eta) * v3;

            nodeNumber[i * yCount + j] = addNode(point, INNER);
        }
    }
    // формирование массива элементов
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            addElement(nodeNumber[i * yCount + j], nodeNumber[(i + 1) * yCount + j], nodeNumber[(i + 1) * yCount + j + 1], nodeNumber[i * yCount + j + 1]);

            pushLayer(layerNumber_);
        }
    }
    if (isFirstRegion_)
    {
        xMin_ = std::min(std::min(v0.x(), v1.x()), std::min(v2.x(), v3.x()));
        xMax_ = std::max(std::max(v0.x(), v1.x()), std::max(v2.x(), v3.x()));
        yMin_ = std::min(std::min(v0.y(), v1.y()), std::min(v2.y(), v3.y()));
        yMax_ = std::max(std::max(v0.y(), v1.y()), std::max(v2.y(), v3.y()));
        isFirstRegion_ = false;
    }
    else
    {
        xMin_ = std::min(xMin_, std::min(std::min(v0.x(), v1.x()), std::min(v2.x(), v3.x())));
        xMax_ = std::max(xMax_, std::max(std::max(v0.x(), v1.x()), std::max(v2.x(), v3.x())));
        yMin_ = std::min(yMin_, std::min(std::min(v0.y(), v1.y()), std::min(v2.y(), v3.y())));
        yMax_ = std::max(yMax_, std::max(std::max(v0.y(), v1.y()), std::max(v2.y(), v3.y())));
    }

    layerNumber_++;
}

void QuadrilateralUnion2D::addTriangleRegion(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2)
{
    Point2D center = (v0 + v1 + v2) / 3.0; // центр треугольника
    Point2D c01 = (v0 + v1) / 2.0; // центр стороны, соединяющей вершину 0 и 1
    Point2D c12 = (v1 + v2) / 2.0; // центр стороны, соединяющей вершину 1 и 2
    Point2D c20 = (v2 + v0) / 2.0; // центр стороны, соединяющей вершину 2 и 0
    UInteger sideCount = count / 2 + 1;
    double h = 2.0 / (double)(sideCount - 1); // шаг изо-сетки
    Point2D quads [][4] = {
        {c20, v0, c01, center},
        {c01, v1, c12, center},
        {c12, v2, c20, center}
    };
    for (int q = 0; q < 3; q++)
    {
        UInteger nodeNumber[sideCount * sideCount];
        for (UInteger i = 0; i < sideCount; i++)
        {
            double xi = -1.0 + (double) i * h;
            for (UInteger j = 0; j < sideCount; j++)
            {
                double eta = -1.0 + (double) j * h;
                Point2D point = isoFunc(0, xi, eta) * quads[q][0]  +
                        isoFunc(1, xi, eta) * quads[q][1] +
                        isoFunc(2, xi, eta) * quads[q][2] +
                        isoFunc(3, xi, eta) * quads[q][3];
                //                if (i == 0 || j == sideCount - 1)
                //                    nodeNumber[i * sideCount + j] = addNode(point, BORDER);
                //                else
                nodeNumber[i * sideCount + j] = addNode(point, INNER);
            }
        }
        // формирование массива элементов
        for (UInteger i = 0; i < sideCount - 1; i++)
        {
            for (UInteger j = 0; j < sideCount - 1; j++)
            {
                addElement(nodeNumber[i * sideCount + j],
                        nodeNumber[(i + 1) * sideCount + j],
                        nodeNumber[(i + 1) * sideCount + j + 1],
                        nodeNumber[i * sideCount + j + 1]);

                pushLayer(layerNumber_);
            }
        }
    }
    if (isFirstRegion_)
    {
        xMin_ = std::min(std::min(v0.x(), v1.x()), v2.x());
        xMax_ = std::max(std::max(v0.x(), v1.x()), v2.x());
        yMin_ = std::min(std::min(v0.y(), v1.y()), v2.y());
        yMax_ = std::max(std::max(v0.y(), v1.y()), v2.y());
        isFirstRegion_ = false;
    }
    else
    {
        xMin_ = std::min(xMin_, std::min(std::min(v0.x(), v1.x()), v2.x()));
        xMax_ = std::max(xMax_, std::max(std::max(v0.x(), v1.x()), v2.x()));
        yMin_ = std::min(yMin_, std::min(std::min(v0.y(), v1.y()), v2.y()));
        yMax_ = std::max(yMax_, std::max(std::max(v0.y(), v1.y()), v2.y()));
    }

    layerNumber_++;
}

void QuadrilateralUnion2D::addMesh(QuadrilateralMesh2D *mesh)
{
    UInteger nodeNumber[mesh->nodesCount()]; // номера узлов в сетке после добавления
    // все узлы добавляютя как внутрение, для востановления типов необходимо выполнить отдельную процедуру
    for (UInteger i = 0; i < mesh->nodesCount(); i++)
    {
        // внутренние узлы добавляем без перебеора (считаем, что наложения областей нет)
        PointPointer point = mesh->node(i);
        UInteger currSize = node_.size(); // текущее количество узлов в сетке
        if (mesh->nodeType(i) == INNER)
            nodeNumber[i] = pushNode(Point2D(point->x(), point->y()), INNER);
        else if (mesh->nodeType(i) == CHARACTER)
            nodeNumber[i] = addNode(Point2D(point->x(), point->y()), CHARACTER);
        else
        {
            nodeNumber[i] = addNode(Point2D(point->x(), point->y()), INNER);
            if (nodeNumber[i] == currSize)
            {
                // граничный узел не был найден в массиве узлов => сохраняем его статус
                node_[currSize].type = BORDER;
            }
        }
    }

    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        ElementPointer element = mesh->element(i);
        UInteger node0 = nodeNumber[element->vertexNode(0)];
        UInteger node1 = nodeNumber[element->vertexNode(1)];
        UInteger node2 = nodeNumber[element->vertexNode(2)];
        UInteger node3 = nodeNumber[element->vertexNode(3)];
        addElement(node0, node1, node2, node3);

        pushLayer(layerNumber_);
    }
    if (isFirstRegion_)
    {
        xMin_ = mesh->xMin();
        xMax_ = mesh->xMax();
        yMin_ = mesh->yMin();
        yMax_ = mesh->yMax();
        isFirstRegion_ = false;
    }
    else
    {
        xMin_ = std::min(xMin_, mesh->xMin());
        xMax_ = std::max(xMax_, mesh->xMax());
        yMin_ = std::min(yMin_, mesh->yMin());
        yMax_ = std::max(yMax_, mesh->yMax());
    }

    std::cout << "Добавлена сетка четырехугольных элементов: узлов - " << mesh->nodesCount() << ", элементов - " << mesh->elementsCount() << "; номер слоя - " << layerNumber_ << "." << std::endl;
    std::cout << "Объединенная дискретная модель: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;

    layerNumber_++;
}

}
