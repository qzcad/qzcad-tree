#include "segmentmesh2d.h"
#include "segment.h"

#include <iostream>
#include <math.h>

namespace msh
{
SegmentMesh2D::SegmentMesh2D() : Mesh2D(NULL)
{
}

SegmentMesh2D::SegmentMesh2D(const SegmentMesh2D &mesh) : Mesh2D(&mesh)
{
    node_ = mesh.node_;
    element_ = mesh.element_;
}

SegmentMesh2D::SegmentMesh2D(const SegmentMesh2D *mesh) : Mesh2D(mesh)
{
    element_ = mesh->element_;
    node_ = mesh->node_;
}

void SegmentMesh2D::functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, bool isOptimized, std::function<double(Point2D, Point2D)> distance)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    double minDistance = 0.4 * sqrt(hx*hx + hy*hy);
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
    std::cout << "Character nodes: " << charPoint.size() << std::endl;
    for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
    {
        pushNode(*cPoint, CHARACTER);
    }
    double x0 = xMin_;
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        double y0 = yMin_;
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            int index = 0;
            Point2D p0 (x0, y0);
            Point2D p1 ((i < xCount - 2) ? (x0 + hx) : xMax_, y0);
            Point2D p2 ((i < xCount - 2) ? (x0 + hx) : xMax_, (j < yCount - 2) ? (y0 + hy) : yMax_);
            Point2D p3 (x0, (j < yCount - 2) ? (y0 + hy) : yMax_);

            if (distance != nullptr) minDistance = 0.4 * distance(p0, p3);

            if (func(p0.x(), p0.y()) < epsilon_) index |= 1;
            if (func(p1.x(), p1.y()) < epsilon_) index |= 2;
            if (func(p2.x(), p2.y()) < epsilon_) index |= 4;
            if (func(p3.x(), p3.y()) < epsilon_) index |= 8;

            if (search_table[index] & 1) border[0] = binary(p0, p1, func);
            if (search_table[index] & 2) border[1] = binary(p1, p2, func);
            if (search_table[index] & 4) border[2] = binary(p2, p3, func);
            if (search_table[index] & 8) border[3] = binary(p3, p0, func);

            for (int ii = 0; edge_table[index][ii] != -1; ii += 2)
            {
                Point2D prev = border[edge_table[index][ii]];
                Point2D next = border[edge_table[index][ii + 1]];
                // упорядоченное добавление узлов
                UInteger ii0 = (prev.x() < next.x() || (prev.x() == next.x() && prev.y() < next.y())) ? addNode(prev, BORDER, minDistance, distance) : addNode(next, BORDER, minDistance, distance);
                UInteger ii1 = (prev.x() < next.x() || (prev.x() == next.x() && prev.y() < next.y())) ? addNode(next, BORDER, minDistance, distance) : addNode(prev, BORDER, minDistance, distance);
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
            double dp0 = (distance == nullptr) ? point.distanceTo(p0) : distance(point, p0);
            double dp1 = (distance == nullptr) ? point.distanceTo(p1) : distance(point, p1);
            double min_d = dp0 * dp0 + dp1 * dp1;
            for (UInteger j = 1; j < elementsCount(); j++)
            {
                p0 = node_[element_[j][0]].point;
                p1 = node_[element_[j][1]].point;
                dp0 = (distance == nullptr) ? point.distanceTo(p0) : distance(point, p0);
                dp1 = (distance == nullptr) ? point.distanceTo(p1) : distance(point, p1);
                double d = dp0 * dp0 + dp1 * dp1;
                if (d < min_d)
                {
                    min_d = d;
                    min_el = j;
                }
            }
            addElement(i, element_[min_el][1]);
            node_[element_[min_el][1]].adjacent.erase(min_el);
            element_[min_el][1] = i;
            node_[i].adjacent.insert(min_el);
        }
    }
    // оптимизация по кривизне границы
    for (int it = 0; (it < 4) && isOptimized; ++it)
    {
        std::cout << "Curvature it " << it << " ";
        isOptimized = false;
        UInteger es = element_.size();
        for (UInteger i = 0; i < es; i++)
        {
            Segment s = element_[i];
            Point2D a = node_[s[0]].point;
            Point2D b = node_[s[1]].point;
            Point2D v(a, b);
            Point2D n = v.perpendicular().normalized();
            Point2D c = 0.5 * (a + b);
            double l = v.length();
            Point2D p0 = (-0.25 * l * n) + c;
            Point2D p1 = (0.25 * l * n) + c;
            if (func(p0.x(), p0.y()) * func(p1.x(), p1.y()) < epsilon_)
            {
                Point2D border = binary(p0, p1, func);
                if (!border.isEqualTo(c, 0.08 * l))
                {
                    UInteger j = pushNode(border, BORDER);
                    addElement(j, s[1]);
                    node_[s[1]].adjacent.erase(i);
                    element_[i][1] = j;
                    node_[j].adjacent.insert(i);
                    std::cout << '+';
                    isOptimized = true;
                }
            }
            else
            {
                std::cout << '-';
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Segments mesh: nodes - " << nodesCount() << " elements - " << elementsCount() << std::endl;
}

void SegmentMesh2D::functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func_a, std::function<double (double, double)> func_b, std::list<Point2D> charPoint, double delta)
{
    clear();
    SegmentMesh2D mesh_a, mesh_b;
    std::list<Point2D> points_a, points_b;
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    // формировка массивов особых точек
    for (std::list<Point2D>::iterator p = charPoint.begin(); p != charPoint.end(); ++p)
    {
        Point2D point = *p;
        if (fabs(func_a(point.x(), point.y())) < epsilon_) points_a.push_back(point);
        else if (fabs(func_b(point.x(), point.y())) < epsilon_) points_b.push_back(point); // общие точки добавятся во время обработки узлов зоны контакта
    }
    mesh_a.functionalDomain(xCount, yCount, xMin, yMin, width, height, func_a, points_a, true);
    // обработка узлов зоны контакта
    for (UInteger i = 0 ; i < mesh_a.nodesCount(); i++)
    {
        Point2D point = mesh_a.node_[i].point;
        if (fabs(func_b(point.x(), point.y())) < epsilon_) points_b.push_back(point);
    }
    mesh_b.functionalDomain(xCount, yCount, xMin, yMin, width, height, func_b, points_b, true);
    node_ = mesh_a.node_;
    element_ = mesh_a.element_;
    for (ElementIterator el_b = mesh_b.element_.begin(); el_b != mesh_b.element_.end(); ++el_b)
    {
        Segment seg = *el_b;
        Node2D n0 = mesh_b.node_[seg[0]];
        Node2D n1 = mesh_b.node_[seg[1]];
        seg[0] = addNode(n0.point, n0.type, 100.0 * epsilon_);
        seg[1] = addNode(n1.point, n1.type, 100.0 * epsilon_);
        // Новый элемент добавляем, если такого элемента нет
        bool exist = false;
        for (ElementIterator el_a = element_.begin(); el_a != element_.end(); ++el_a)
            if (seg.isSame(*el_a))
            {
                exist = true;
                break;
            }
        if (!exist)
            addElement(seg);
    }
    // сгущение сетки в окрестности контакта
    UInteger baseElementsCount = elementsCount();
    for (UInteger i = 0; i < baseElementsCount; i++)
    {
        // во все элементы, находящиеся в окрестности контакта добавляем узел в серидину
        Segment seg = element_[i];
        Point2D a = node_[seg[0]].point;
        Point2D b = node_[seg[1]].point;
        double val_a_a = func_a(a.x(), a.y()), val_a_b = func_b(a.x(), a.y());
        double val_b_a = func_a(b.x(), b.y()), val_b_b = func_b(b.x(), b.y());
        if (sqrt(val_a_a*val_a_a + val_a_b*val_a_b) < delta || sqrt(val_b_a*val_b_a + val_b_b*val_b_b) < delta)
        {
            Point2D v(a, b);
            Point2D n = v.perpendicular().normalized();
            Point2D c = 0.5 * (a + b);
            double l = v.length();
            Point2D p0 = (-0.25 * l * n) + c;
            Point2D p1 = (0.25 * l * n) + c;
            if (func_a(p0.x(), p0.y()) * func_a(p1.x(), p1.y()) < epsilon_)
            {
                Point2D border = binary(p0, p1, func_a);
                UInteger ic = pushNode(border, CHARACTER);
                addElement(ic, seg[1]);
                node_[seg[1]].adjacent.erase(i);
                element_[i][1] = ic;
                node_[ic].adjacent.insert(i);
            }
        }
    }

    std::cout << "Contact mesh: nodes - " << nodesCount() << " elements - " << elementsCount() << std::endl;
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

void SegmentMesh2D::addElement(const Segment &segment)
{
    element_.push_back(segment);
    // обновление списка смежных узлов
    node_[segment[0]].adjacent.insert(element_.size() - 1);
    node_[segment[1]].adjacent.insert(element_.size() - 1);
}

void SegmentMesh2D::addElement(const UInteger &node0, const UInteger &node1)
{
    Segment segment(node0, node1);
    addElement(segment);
}

void SegmentMesh2D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0], nodes_ref[1]);
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

Segment SegmentMesh2D::segment(const UInteger &number) const
{
    return element_[number];
}

void SegmentMesh2D::clearElements()
{
    element_.clear();
}

bool SegmentMesh2D::isCrossedElement(const Point2D &p0, const Point2D &p1, UInteger &number)
{
    double x[2][2], y[2][2];
    double A[2], B[2], C[2];
    double det;
    double t[2];
    x[0][0] = p0.x();
    y[0][0] = p0.y();
    x[1][0] = p1.x();
    y[1][0] = p1.y();
    for (ElementIterator el = element_.begin(); el != element_.end(); ++el)
    {
        Segment segment = *el;
        Point2D p0 = node_[segment[0]].point;
        Point2D p1 = node_[segment[1]].point;
        x[0][1] = p0.x();
        y[0][1] = p0.y();
        x[1][1] = p1.x();
        y[1][1] = p1.y();
        A[0] = x[1][0] - x[0][0];
        A[1] = y[1][0] - y[0][0];
        B[0] = x[0][1] - x[1][1];
        B[1] = y[0][1] - y[1][1];
        C[0] = x[0][1] - x[0][0];
        C[1] = y[0][1] - y[0][0];
        det = A[0] * B[1] - B[0] * A[1];
        if (fabs(det) > epsilon_)
        {
            t[0] = (C[0] * B[1] - B[0] * C[1]) / det;
            t[1] = (A[0] * C[1] - C[0] * A[1]) / det;
            if (t[0] > - epsilon_ && t[0] < 1.0 + epsilon_ && t[1] > epsilon_ && t[1] < 1.0 - epsilon_)
            {
                // точка пересечения находится строго можду началом у концом второго отрезка (возможно совпадение с началом или концом первого отрезка)
                number = el - element_.begin();
                return true;
            }
        }
    }
    return false;
}

Point2D SegmentMesh2D::refineMidpoint(const UInteger &number, std::function<double(double, double)> func)
{
    Segment seg = element_[number];
    Point2D a = node_[seg[0]].point;
    Point2D b = node_[seg[1]].point;
    Point2D v(a, b);
    Point2D n = v.perpendicular().normalized();
    Point2D c = 0.5 * (a + b);
    double l = v.length();
    Point2D p0 = ((-0.25 * l) * n) + c;
    Point2D p1 = ((0.25 * l) * n) + c;
    Point2D border = c;
    if (func != nullptr && func(p0.x(), p0.y()) * func(p1.x(), p1.y()) < epsilon_)
    {
        border = binary(p0, p1, func);
    }
    else if (func != nullptr)
    {
        std::cout << "l = " << l << std::endl;
        std::cout << "c: " << c.x() << " " << c.y() << std::endl;
        std::cout << p0.x() << ", " << p0.y() << "; " << p1.x() << ", " << p1.y() << ": "<< func(p0.x(), p0.y()) << " * " << func(p1.x(), p1.y()) << std::endl;
    }
    UInteger ic = pushNode(border, BORDER);
    addElement(ic, seg[1]);
    node_[seg[1]].adjacent.erase(number);
    element_[number][1] = ic;
    node_[ic].adjacent.insert(number);
    return border;
}

bool SegmentMesh2D::isEncroached(const Point2D &point, UInteger &number)
{
    for (ElementIterator el = element_.begin(); el != element_.end(); ++el)
    {
        Segment segment = *el;
        Point2D p0 = node_[segment[0]].point;
        Point2D p1 = node_[segment[1]].point;
        Point2D c = 0.5 * (p0 + p1);
        double r = p0.distanceTo(p1) / 2.0;
        if (c.distanceTo(point) < r - epsilon_)
        {
            number = el - element_.begin();
            return true;
        }
    }
    return false;
}

}
