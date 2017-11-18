#include "segmentmesh2d.h"
#include "segment.h"

#include <iostream>
#include <limits>
#include <math.h>
#include <time.h>

#include "consoleprogress.h"

namespace msh
{
SegmentMesh2D::SegmentMesh2D() : Mesh2D(NULL)
{
}

SegmentMesh2D::SegmentMesh2D(const SegmentMesh2D &mesh) : Mesh2D(&mesh)
{
    element_ = mesh.element_;
}

SegmentMesh2D::SegmentMesh2D(const SegmentMesh2D *mesh) : Mesh2D(mesh)
{
    element_ = mesh->element_;
}

void SegmentMesh2D::MarchingQuads(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, double level, int smooth, int optimize, std::function<double(Point2D, Point2D)> distance)
{
    clear();
    std::cout << "Mesh generation: level = " << level << std::endl;
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;

    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
//    const double tolerance = 0.2;

    ConsoleProgress progress(xCount - 1);

    double x0 = xMin_;
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        double y0 = yMin_;
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            Point2D p0 (x0, y0);
            Point2D p1 ((i < xCount - 2) ? (x0 + hx) : xMax_, y0);
            Point2D p2 ((i < xCount - 2) ? (x0 + hx) : xMax_, (j < yCount - 2) ? (y0 + hy) : yMax_);
            Point2D p3 (x0, (j < yCount - 2) ? (y0 + hy) : yMax_);

            cellContours(p0, p1, p2, p3,
                         func(p0.x(), p0.y()),
                         func(p1.x(), p1.y()),
                         func(p2.x(), p2.y()),
                         func(p3.x(), p3.y()),
                         func, level, distance);

            y0 += hy;
        }
        x0 += hx;
        ++progress;
    }
    // учет характерных точек
    std::cout << "Character nodes: " << charPoint.size() << std::endl;
    for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
    {
        Point2D character = *cPoint;
        double d = (distance == nullptr) ? character.distanceTo(node_[0].point) : distance(character, node_[0].point);
        UInteger num = 0;
        character.print();
        std::cout << "F = " << func(character.x(), character.y()) << std::endl;
        for (UInteger i = 1; i < nodesCount(); i++)
        {
            double c = ((distance == nullptr) ? character.distanceTo(node_[i].point) : distance(character, node_[i].point));
            if (c < d)
            {
                d = c;
                num = i;
            }
        }
        node_[num].point = character;
        node_[num].type = CHARACTER;
    }
    laplacianSmoothing(func, level, smooth);
    curvatureSmoothing(func, level, 0.05, optimize);
    distlenSmoothing(func, level, optimize);
    std::cout << "Segments mesh: nodes - " << nodesCount() << " elements - " << elementsCount() << std::endl;
}

void SegmentMesh2D::MarchingQuads(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func_a, std::function<double (double, double)> func_b, std::list<Point2D> charPoint, double delta)
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
    mesh_a.MarchingQuads(xCount, yCount, xMin, yMin, width, height, func_a, points_a, true);
    // обработка узлов зоны контакта
    for (UInteger i = 0 ; i < mesh_a.nodesCount(); i++)
    {
        Point2D point = mesh_a.node_[i].point;
        if (fabs(func_b(point.x(), point.y())) < epsilon_) points_b.push_back(point);
    }
    mesh_b.MarchingQuads(xCount, yCount, xMin, yMin, width, height, func_b, points_b, true);
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

void SegmentMesh2D::contourGraph(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, int contours, int smooth, int optimize)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    double maxf = 0.0;
    double h = 0.0;

    MarchingQuads(xCount, yCount, xMin, yMin, width, height, func, charPoint, 0.0, smooth, optimize);

    double x0 = xMin_;
    for (UInteger i = 0; i < xCount; i++)
    {
        double y0 = yMin_;
        for (UInteger j = 0; j < yCount; j++)
        {
            double f = func(x0, y0);
            if (f > maxf)
                maxf = f;
            y0 += hy;
        }
        x0 += hx;
    }
    std::cout << "max(f) = " << maxf << std::endl;

    h = maxf / (double)contours;

    for (int i = 1; i < contours; i++)
    {
        SegmentMesh2D S;
        std::list<Point2D> C;
        S.MarchingQuads(xCount, yCount, xMin, yMin, width, height, func, C, 0.0 + (double)i * h, smooth, optimize);
        for (UInteger j = 0; j < S.elementsCount(); j++)
        {
            addElement(pushNode(S.node(S.element_[j][0]), S.nodeType(S.element_[j][0])), pushNode(S.node(S.element_[j][1]), S.nodeType(S.element_[j][1])));
        }
    }

    evalNodalValues(func);

    std::cout << "Segments mesh: nodes - " << nodesCount() << " elements - " << elementsCount() << std::endl;
}

void SegmentMesh2D::frontGraph(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, int contours, int smooth, int optimize)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    double maxf = 0.0;
//    double h = 0.0;

    MarchingQuads(xCount, yCount, xMin, yMin, width, height, func, charPoint, 0.0, smooth, optimize);

    UInteger border_elements_count = elementsCount();

    auto contour_func = [&](double x, double y)
    {
        Point2D point(x, y);
        double min_distance = std::numeric_limits<double>::max();
        double f = func(x, y);
        for (UInteger ie = 0; ie < border_elements_count; ie++)
        {
            Segment s = element_[ie];
            Point2D p1 = node_[s[0]].point;
            Point2D p2 = node_[s[1]].point;
            double a = p2.y() - p1.y();
            double b = p1.x() - p2.x();
            double c = p2.x() * p1.y() - p2.y() * p1.x();
            double xp = (b * (b * x - a * y) - a * c) / (a*a + b*b);
            double yp = (a * (-b * x + a * y) - b * c) / (a*a + b*b);
            double t = (fabs(p2.x() - p1.x()) > epsilon_) ? (xp - p1.x()) / (p2.x() - p1.x()) : (yp - p1.y()) / (p2.y() - p1.y());
            double distance = 0.0;
            if (t < 0.0 || t > 1.0)
                distance = std::min(p1.distanceTo(point), p2.distanceTo(point));
            else
                distance = fabs(a * x + b * y + c) / sqrt(a*a + b*b);
            if (min_distance > distance)
                min_distance = distance;
        }
        if (f < 0.0)
            return -1.0 * min_distance;
        else if (f > 0.0)
            return min_distance;
        return 0.0;
    };

    double x0 = xMin_;
    for (UInteger i = 0; i < xCount; i++)
    {
        double y0 = yMin_;
        for (UInteger j = 0; j < yCount; j++)
        {
            double f = contour_func(x0, y0);
            if (f > maxf)
                maxf = f;
            y0 += hy;
        }
        x0 += hx;
    }
    std::cout << "max(f) = " << maxf << std::endl;

    double h = maxf / (double)contours;

    for (int i = 1; i < contours; i++)
    {
        SegmentMesh2D S;
        std::list<Point2D> C;
        S.MarchingQuads(xCount, yCount, xMin, yMin, width, height, contour_func, C, (double)i * h, smooth, optimize);
        for (UInteger j = 0; j < S.elementsCount(); j++)
        {
            addElement(pushNode(S.node(S.element_[j][0]), S.nodeType(S.element_[j][0])), pushNode(S.node(S.element_[j][1]), S.nodeType(S.element_[j][1])));
        }
    }

    std::vector<double> r(nodesCount());
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point2D p = node_[i].point;
        r[i] = contour_func(p.x(), p.y());
    }
    addDataVector("R", r);

    std::cout << "Segments mesh: nodes - " << nodesCount() << " elements - " << elementsCount() << std::endl;
}

void SegmentMesh2D::parametricDomain(const UInteger &count, const double &tmin, const double &tmax, std::function<Point2D (double)> domainFunction)
{
    double h = (tmax - tmin) / (double)(count - 1);
    double t = tmin;
    std::vector<double> parametric_nodes;
    bool isClosed = (domainFunction(tmin)).isEqualTo(domainFunction(tmax), epsilon_);

    clear();

    for (UInteger i = 0; i < count - 1; i++)
    {
        parametric_nodes.push_back(t);
        t += h;
    }
    if (!isClosed) parametric_nodes.push_back(tmax);
    for (std::vector<double>::iterator it = parametric_nodes.begin(); it != parametric_nodes.end(); )
    {
        double t0 = *it;
        double t1 = ((it + 1) != parametric_nodes.end()) ? *(it + 1) : tmax;
        Point2D p0 = domainFunction(t0);
        Point2D p1 = domainFunction(t1);
        Point2D c = 0.5 * (p0 + p1);
        Point2D s = domainFunction(0.5 * (t0 + t1));
        if (c.distanceTo(s) / p0.distanceTo(p1) > 0.02)
        {
            parametric_nodes.insert(it + 1, 0.5 * (t0 + t1));
            it = parametric_nodes.begin();
            std::cout << c.distanceTo(s) / p0.distanceTo(p1) << std::endl;
        }
        else ++it;
    }
    for (std::vector<double>::iterator it = parametric_nodes.begin(); it != parametric_nodes.end(); ++it)
    {
        pushNode(domainFunction(*it), BORDER);
    }
    for (UInteger i = 0; i < nodesCount() - 1; i++)
    {
        addElement(i, i + 1);
    }
    if (isClosed) addElement(nodesCount() - 1, 0);
    std::cout << "Segments mesh: nodes - " << nodesCount() << " elements - " << elementsCount() << std::endl;
    updateDomain();
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

Point2D SegmentMesh2D::refineMidpoint(const UInteger &number, std::function<double(double, double)> func, double level)
{
    Segment seg = element_[number];
    Point2D a = node_[seg[0]].point;
    Point2D b = node_[seg[1]].point;
    Point2D border = findBorder(a, b, func, 0.5, level);
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

void SegmentMesh2D::laplacianSmoothing(std::function<double(double, double)> func, double level, int iter_num)
{
    if (func == nullptr)
        return;

    std::cout << "Laplacian Smoothing: " << nodesCount() << " nodes, " << elementsCount() << " elements." << std::endl;
    clock_t start = clock();
    for (int iter = 0; iter < iter_num; iter++)
    {
        std::cout << "Laplacian smothing: " << iter << std::endl;
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Node2D node = node_[i];
            Point2D point = node.point;
            AdjacentSet adjacent = node.adjacent;
            if (node.type != CHARACTER && adjacent.size() == 2 && fabs(func(point.x(), point.y()) - level) < epsilon_)
            {
                AdjacentSet::iterator it = adjacent.begin();
                Segment s0 = element_[*(it)];
                Segment s1 = element_[*(++it)];
                Point2D a = ((s0[0] == i) ? node_[s0[1]].point : node_[s0[0]].point);
                Point2D b = ((s1[0] == i) ? node_[s1[1]].point : node_[s1[0]].point);
                Point2D c = 0.5 * (a + b);
                Point2D v(a, b);
                double l = v.length();
                Point2D border = findBorder(c, func, 0.1 * l, level);
                node_[i].point = border;
            }
        }
    }
    std::cout << "Done in " << (double)(clock() - start) / CLOCKS_PER_SEC << "s." << std::endl;
}

void SegmentMesh2D::distlenSmoothing(std::function<double (double, double)> func, double level, int iter_num)
{
    std::cout << "Distance-Length Optimization: " << nodesCount() << " nodes, " << elementsCount() << " elements." << std::endl;

    auto functor = [&](const Point2D &a, const Point2D &o, const Point2D &b)
    {
        const double alpha = 0.001;
        const double beta = 1.0 - alpha;
        double d0 = o.distanceTo(a);
        double d1 = o.distanceTo(b);
        double e0c = distToBorder(o, a, func, 0.5, level);
        double e1c = distToBorder(o, b, func, 0.5, level);
        return beta * (e0c*e0c + e1c*e1c) + alpha * (d0*d0 + d1*d1);
    };
    clock_t start = clock();
    for (int iter = 0; iter < iter_num; ++iter)
    {
        ConsoleProgress progress(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Node2D node = node_[i];
            if (node.adjacent.size() == 2 && node.type != CHARACTER)
            {
                UInteger e0 = *(node.adjacent.begin());
                UInteger e1 = *(++node.adjacent.begin());
                Node2D node0 = (element_[e0][0] == i) ? node_[element_[e0][1]] : node_[element_[e0][0]];
                Node2D node1 = (element_[e1][0] == i) ? node_[element_[e1][1]] : node_[element_[e1][0]];
                double f_current = functor(node0.point, node.point, node1.point);
                bool optimized = false;
                double step = 0.1;
                short iic = 0;
                do {
                    Point2D border = findBorder(node.point, node0.point, func, -step, level);
                    double f_dir = functor(node0.point, border, node1.point);
                    if (f_dir < f_current)
                    {
                        node.point = border;
                        f_current = f_dir;
                        optimized = true;
                    }
                    else
                    {
                        step /= 10.0;
                    }
                    ++iic;
                } while (step >= 0.0001 && iic < 100);
                if (!optimized)
                {
                    step = 0.1;
                    do {
                        Point2D border = findBorder(node.point, node0.point, func, step, level);
                        double f_dir = functor(node0.point, border, node1.point);
                        if (f_dir < f_current)
                        {
                            node.point = border;
                            f_current = f_dir;
                            optimized = true;
                        }
                        else
                        {
                            step /= 10.0;
                        }
                        ++iic;
                    } while (step >= 0.0001 && iic < 100);
                }
                if (!optimized)
                {
                    step = 0.1;
                    do {
                        Point2D border = findBorder(node.point, node1.point, func, -step, level);
                        double f_dir = functor(node0.point, border, node1.point);
                        if (f_dir < f_current)
                        {
                            node.point = border;
                            f_current = f_dir;
                            optimized = true;
                        }
                        else
                        {
                            step /= 10.0;
                        }
                        ++iic;
                    } while (step >= 0.0001 && iic < 100);
                }
                if (!optimized)
                {
                    step = 0.1;
                    do {
                        Point2D border = findBorder(node.point, node1.point, func, step, level);
                        double f_dir = functor(node0.point, border, node1.point);
                        if (f_dir < f_current)
                        {
                            node.point = border;
                            f_current = f_dir;
                            optimized = true;
                        }
                        else
                        {
                            step /= 10.0;
                        }
                        ++iic;
                    } while (step >= 0.0001 && iic < 100);
                }
                node_[i].point = node.point;
//                Point2D border0 = findBorder(node.point, node0.point, func, -0.001, level);
//                bool isOptimized = false;
//                double f0 = functor(node0.point, border0, node1.point);
//                while (f0 < fcurrent)
//                {
//                    fcurrent = f0;
//                    node.point = border0;
//                    border0 = findBorder(node.point, node0.point, func, -0.001, level);
//                    f0 = functor(node0.point, border0, node1.point);
//                    isOptimized = true;
//                }
//                if (!isOptimized)
//                {
//                    border0 = findBorder(node.point, node1.point, func, -0.001, level);
//                    f0 = functor(node0.point, border0, node1.point);
//                    while (f0 < fcurrent)
//                    {
//                        fcurrent = f0;
//                        node.point = border0;
//                        border0 = findBorder(node.point, node1.point, func, -0.001, level);
//                        f0 = functor(node0.point, border0, node1.point);
//                    }
//                }
//                node_[i].point = node.point;
            }
            ++progress;
        }
    }
    std::cout << "Done in " << (double)(clock() - start) / CLOCKS_PER_SEC << "s." << std::endl;
}

void SegmentMesh2D::curvatureSmoothing(std::function<double (double, double)> func, double level, double alpha, int iter_num)
{
    bool optimized = true;
    std::cout << "Curvature Post Smothing: " << nodesCount() << " nodes, " << elementsCount() << " elements." << std::endl;
    // оптимизация по кривизне границы
    clock_t start = clock();
    for (int it = 0; it < iter_num && optimized; ++it)
    {
        std::cout << "Curvature it " << it << " ";
        optimized = false;
        UInteger es = element_.size();
        for (UInteger i = 0; i < es; i++)
        {
            Segment s = element_[i];
            Point2D a = node_[s[0]].point;
            Point2D b = node_[s[1]].point;
            Point2D c = 0.5 * (a + b);
            double l = a.distanceTo(b);
            Point2D border = findBorder(c, func, 0.5 * l, level);
            if (!border.isEqualTo(c, 0.05 * l))
            {
                UInteger j = pushNode(border, BORDER);
                addElement(j, s[1]);
                node_[s[1]].adjacent.erase(i);
                element_[i][1] = j;
                node_[j].adjacent.insert(i);
                std::cout << '+';
                optimized = true;
            }
        }
        std::cout << std::endl;
    }
    std::cout << "Done in " << (double)(clock() - start) / CLOCKS_PER_SEC << "s: " << nodesCount() << " nodes, " << elementsCount() << " elements." << std::endl;
}

double SegmentMesh2D::cfunction(const double &x, const double &y)
{
    Point2D point(x, y);
    double min_distance = std::numeric_limits<double>::max();
    double sign = +1.0;
    for (UInteger ie = 0; ie < elementsCount(); ie++)
    {
        Segment s = element_[ie];
        Point2D p1 = node_[s[0]].point;
        Point2D p2 = node_[s[1]].point;
        double distance = point.distanceTo(p1, p2);
        if (min_distance > distance)
            min_distance = distance;
    }
    // определение знака
    UInteger count = 0;
    Point2D next = point + Point2D(1.0, 0.0);
    for (UInteger je = 0; je < elementsCount(); je++)
    {
        Segment sj = element_[je];
        Point2D q1 = node_[sj[0]].point;
        Point2D q2 = node_[sj[1]].point;
        double tp = -1.0, tq = -1.0;
        if ( isCrossed(point, next, q1, q2, tp, tq) && tp > 0.0 && ((q1.y() < y && y <= q2.y()) || (q2.y() < y && y <= q1.y())) )
            count++;
    }
    sign = (count % 2 == 0 ) ? -1.0 : 1.0;
    return min_distance * sign;
}

void SegmentMesh2D::cellContours(const Point2D &p0, const Point2D &p1, const Point2D &p2, const Point2D &p3, const double &v0, const double &v1, const double &v2, const double &v3, std::function<double (double, double)> func, double level, std::function<double(Point2D, Point2D)> distance)
{
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

    double min_distance = 0.0;
    int index = 0;

    if (distance != nullptr)
        min_distance = 0.2 * distance(p0, p3);
    else
        min_distance = 0.2 * p0.distanceTo(p3);

    if (v0 - level < epsilon_) index |= 1;
    if (v1 - level < epsilon_) index |= 2;
    if (v2 - level < epsilon_) index |= 4;
    if (v3 - level < epsilon_) index |= 8;

    if (index == 5 || index == 10)
    {
        Point2D c = 0.25 * (p0 + p1 + p2 + p3);
        Point2D p01 = 0.5 * (p0 + p1);
        Point2D p12 = 0.5 * (p1 + p2);
        Point2D p23 = 0.5 * (p2 + p3);
        Point2D p30 = 0.5 * (p3 + p0);
        double valc = func(c.x(), c.y());
        double val01 = func(p01.x(), p01.y());
        double val12 = func(p12.x(), p12.y());
        double val23 = func(p23.x(), p23.y());
        double val30 = func(p30.x(), p30.y());
        cellContours(p0, p01, c, p30,
                     v0, val01, valc, val30, func, level, distance);
        cellContours(p1, p12, c, p01,
                     v1, val12, valc, val01, func, level, distance);
        cellContours(p2, p23, c, p12,
                     v2, val23, valc, val12, func, level, distance);
        cellContours(p3, p30, c, p23,
                     v3, val30, valc, val23, func, level, distance);
        return;
    }

    if (search_table[index] & 1) border[0] = binary(p0, p1, func, level);
    if (search_table[index] & 2) border[1] = binary(p1, p2, func, level);
    if (search_table[index] & 4) border[2] = binary(p2, p3, func, level);
    if (search_table[index] & 8) border[3] = binary(p3, p0, func, level);

    for (int ii = 0; edge_table[index][ii] != -1; ii += 2)
    {
        Point2D prev = border[edge_table[index][ii]];
        Point2D next = border[edge_table[index][ii + 1]];
        // упорядоченное добавление узлов
        if ((distance != nullptr && distance(prev, next) >= min_distance) || (distance == nullptr && !prev.isEqualTo(next, min_distance)))
        {
            UInteger ii0 = (prev < next) ? addNode(prev, BORDER, min_distance, distance) : addNode(next, BORDER, min_distance, distance);
            UInteger ii1 = (prev < next) ? addNode(next, BORDER, min_distance, distance) : addNode(prev, BORDER, min_distance, distance);
            if (ii0 != ii1)
            {
                addElement(ii0, ii1);
            }
        }
    }
}

}
