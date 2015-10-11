#include "segmentmesh2d.h"
#include "segment.h"

#include <iostream>
#include <math.h>

namespace msh
{
SegmentMesh2D::SegmentMesh2D()
{
    xMax_ = yMax_ = 1.0;
    xMin_ = yMin_ = -1.0;
}

SegmentMesh2D::SegmentMesh2D(const SegmentMesh2D &mesh)
{
    xMin_ = mesh.xMin_;
    xMax_ = mesh.xMax_;
    yMin_ = mesh.yMin_;
    yMax_ = mesh.yMax_;
    node_ = mesh.node_;
    element_ = mesh.element_;
}

SegmentMesh2D::SegmentMesh2D(const SegmentMesh2D *mesh)
{
    element_ = mesh->element_;
    node_ = mesh->node_;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
}

SegmentMesh2D::SegmentMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint)
{
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    const double minDistance = 0.4 * sqrt(hx*hx + hy*hy);
    Point2D border[4]; // массив координат пересечения с границей области
    int edge_table[16][5] = {
        { -1, -1, -1, -1, -1 }, // 0
        {  3,  0, -1, -1, -1 }, // 1
        {  0,  1, -1, -1, -1 }, // 2
        {  3,  1, -1, -1, -1 }, // 3
        {  1,  2, -1, -1, -1 }, // 4
        {  1,  0,  3,  2, -1 }, // 5
        {  0,  2, -1, -1, -1 }, // 6
        {  3,  2, -1, -1, -1 }, // 7
        {  2,  3, -1, -1, -1 }, // 8
        {  2,  0, -1, -1, -1 }, // 9
        {  0,  3,  2,  1, -1 }, // 10
        {  2,  1, -1, -1, -1 }, // 11
        {  1,  3, -1, -1, -1 }, // 12
        {  1,  0, -1, -1, -1 }, // 13
        {  0,  3, -1, -1, -1 }, // 14
        { -1, -1, -1, -1, -1 }  // 15
    };
    int search_table[16] = {
        0, // 0
        1 | 8, // 1
        1 | 2, // 2
        2 | 8, // 3
        2 | 4, // 4
        1 | 2 | 4 | 8, // 5
        1 | 4, // 6
        4 | 8, // 7
        4 | 8, // 8
        1 | 4, // 9
        1 | 2 | 4 | 8, // 10
        2 | 4, // 11
        2 | 8, // 12
        1 | 2, // 13
        1 | 8, // 14
        0 //15
    };
    for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
    {
        pushNode(*cPoint, CHARACTER);
    }
    double x0 = xMin;
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        double y0 = yMin;
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            int index = 0;
            Point2D p0 (x0, y0);
            Point2D p1 (x0 + hx, y0);
            Point2D p2 (x0 + hx, y0 + hy);
            Point2D p3 (x0, y0 + hy);

            if (func(p0.x(), p0.y()) <= 0.0) index |= 1;
            if (func(p1.x(), p1.y()) <= 0.0) index |= 2;
            if (func(p2.x(), p2.y()) <= 0.0) index |= 4;
            if (func(p3.x(), p3.y()) <= 0.0) index |= 8;

            if (search_table[index] & 1) border[0] = binary(p0, p1, func);
            if (search_table[index] & 2) border[1] = binary(p1, p2, func);
            if (search_table[index] & 4) border[2] = binary(p2, p3, func);
            if (search_table[index] & 8) border[3] = binary(p3, p0, func);

            for (int ii = 0; edge_table[index][ii] != -1; ii += 2)
            {
                Point2D prev = border[edge_table[index][ii]];
                Point2D next = border[edge_table[index][ii + 1]];
                UInteger ii0 = addNode(prev, BORDER, minDistance);
                UInteger ii1 = addNode(next, BORDER, minDistance);
                if (ii0 != ii1)
                {
                    addElement(ii0, ii1);
                }
            }

            y0 += hy;
        }
        x0 += hx;
    }
    // учет характерных точек, для которых не нашлось пары
    for (UInteger i = 0; i < charPoint.size(); i++)
    {
        Node2D node = node_[i];
        if (node.adjacent.size() == 0)
        {
            Point2D point = node.point;
            UInteger min_el = 0;
            Point2D p0 = node_[element_[0][0]].point;
            Point2D p1 = node_[element_[0][1]].point;
            double min_d = point.distanceTo(p0) * point.distanceTo(p0) +
                    point.distanceTo(p1) * point.distanceTo(p1);
            for (UInteger j = 1; j < elementsCount(); j++)
            {
                p0 = node_[element_[j][0]].point;
                p1 = node_[element_[j][1]].point;
                double d = point.distanceTo(p0) * point.distanceTo(p0) +
                        point.distanceTo(p1) * point.distanceTo(p1);
                if (d < min_d)
                {
                    min_d = d;
                    min_el = j;
                }
            }
            addElement(i, element_[min_el][1]);
            AdjacentSet a1 = node_[element_[min_el][1]].adjacent;
            a1.erase(min_el);
            element_[min_el][1] = i;
            node_[i].adjacent.insert(min_el);
        }
    }
}

UInteger SegmentMesh2D::elementsCount() const
{
    return element_.size();
}

ElementPointer SegmentMesh2D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

void SegmentMesh2D::addElement(const UInteger &node0, const UInteger &node1)
{
    Segment segment(node0, node1);
    element_.push_back(segment);
    // обновление списка смежных узлов
    node_[node0].adjacent.insert(element_.size() - 1);
    node_[node1].adjacent.insert(element_.size() - 1);
}

void SegmentMesh2D::directionChange()
{
    // to do
}

double SegmentMesh2D::area(const UInteger &number) const
{
    return node_[element_[number].vertexNode(0)].point.distanceTo(node_[element_[number].vertexNode(1)].point);
}

bool SegmentMesh2D::isBorderElement(const UInteger &number) const
{
    return true;
}

Point2D SegmentMesh2D::binary(Point2D p0, Point2D p1, std::function<double (double, double)> func)
{
    double val0 = func(p0.x(), p0.y());
    double val1 = func(p1.x(), p1.y());
    if (val0 * val1 > 0)
    {
        return Point2D(); // значения в узлах отрезка одного знака => нет решения
    }
    if (fabs(val0) < epsilon_) return p0;
    if (fabs(val1) < epsilon_) return p1;
    Point2D center;
    double val;
    do
    {
        center = 0.5 * (p0 + p1);
        val = func(center.x(), center.y());
        if (val0 * val < 0)
        {
            p1 = center;
            val1 = val;
        }
        else
        {
            p0 = center;
            val0 = val;
        }
    } while (!p0.isEqualTo(p1, epsilon_) && fabs(val) > epsilon_);
    return center;
}

}
