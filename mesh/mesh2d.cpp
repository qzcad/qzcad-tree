#include "mesh2d.h"
#include <ctime>
#include <iostream>
#include <math.h>
namespace msh
{
Mesh2D::Mesh2D(const Mesh2D *mesh) : Mesh(mesh)
{
    if (mesh != NULL)
    {
        xMin_ = mesh->xMin_;
        yMin_ = mesh->yMin_;
        xMax_ = mesh->xMax_;
        yMax_ = mesh->yMax_;
        node_ = mesh->node_;
    }
    else
    {
        xMin_ = yMin_ = -1.0;
        xMax_ = yMax_ = 1.0;
    }
}

UInteger Mesh2D::nodesCount() const
{
    return node_.size();
}

PointPointer Mesh2D::node(const UInteger &number) const
{
    PointPointer pointPtr = &node_[number].point;
    return pointPtr;
}

double Mesh2D::xMin() const
{
    return xMin_;
}

double Mesh2D::xMax() const
{
    return xMax_;
}

double Mesh2D::yMin() const
{
    return yMin_;
}

double Mesh2D::yMax() const
{
    return yMax_;
}

double Mesh2D::zMin() const
{
    return -std::min(std::max(fabs(xMax_), fabs(xMin_)), std::max(fabs(yMax_), fabs(yMin_)));
}

double Mesh2D::zMax() const
{
    return std::min(std::max(fabs(xMax_), fabs(xMin_)), std::max(fabs(yMax_), fabs(yMin_)));
}

int Mesh2D::dimesion() const
{
    return 2;
}

NodeType Mesh2D::nodeType(const UInteger &number) const
{
    return node_[number].type;
}

void Mesh2D::clearNodes()
{
    node_.clear();
}

void Mesh2D::flipVertically()
{
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point2D point = node_[i].point;
        node_[i].point.set(point.x(), yMin_ + yMax_ - point.y());
    }
    directionChange();
}

void Mesh2D::flipHorizontally()
{
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point2D point = node_[i].point;
        node_[i].point.set(xMin_ + xMax_ - point.x(), point.y());
    }
    directionChange();
}

void Mesh2D::mirrorVertically()
{
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point2D point = node_[i].point;
        node_[i].point.set(point.x(), -point.y());
    }
    yMax_ = -yMax_;
    yMin_ = -yMin_;
    std::swap(yMin_, yMax_);
    directionChange();
}

void Mesh2D::mirrorHorizontally()
{
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point2D point = node_[i].point;
        node_[i].point.set(-point.x(), point.y());
    }
    xMin_ = -xMin_;
    xMax_ = -xMax_;
    std::swap(xMin_, xMax_);
    directionChange();
}

double Mesh2D::area() const
{
    double sum = 0.0;
    for (UInteger i = 0; i < elementsCount(); i++)
        sum += area(i);
    return sum;
}

UInteger Mesh2D::pushNode(const Point2D &point, const NodeType &type)
{
    Node2D node;
    node.point = point;
    node.type = type;
    node_.push_back(node);
    return node_.size() - 1;
}

UInteger Mesh2D::pushNode(PointPointer point, const NodeType &type)
{
    return pushNode(Point2D(point->x(), point->y()), type);
}

UInteger Mesh2D::addNode(const Point2D &point, const NodeType &type, double epsilon, std::function<double(Point2D, Point2D)> distance)
{
    UInteger ns = node_.size();
    for (UInteger i = 0; i < ns; i++)
    {
        UInteger ii = ns - i - 1UL;
        if ( (distance == nullptr && point.isEqualTo(node_[ii].point, epsilon)) ||
             (distance != nullptr && distance(point, node_[ii].point) < epsilon))
        {
            if ((node_[ii].type == INNER || node_[ii].type == BORDER) && (type == BORDER || type == CHARACTER)) // обновление типа узла
                node_[ii].type = type;
            return ii;
        }
    }
    return pushNode(point, type);
}

UInteger Mesh2D::addNode(const Node2D &node, double epsilon, std::function<double (Point2D, Point2D)> distance)
{
    return addNode(node.point, node.type, epsilon, distance);
}

Point3D Mesh2D::normal(const UIntegerVector &face) const
{
    return Point3D(0.0, 0.0, -1.0);
}

void Mesh2D::updateDomain()
{
    if (nodesCount() > 0)
    {
        Point2D start = node_[0].point;
        xMin_ = xMax_ = start.x();
        yMin_ = yMax_ = start.y();
        for (UInteger i = 1; i < node_.size(); i++)
        {
            Point2D current = node_[i].point;
            double x = current.x();
            double y = current.y();
            if (xMin_ > x) xMin_ = x;
            if (xMax_ < x) xMax_ = x;
            if (yMin_ > y) yMin_ = y;
            if (yMax_ < y) yMax_ = y;
        }
    }
}

UInteger Mesh2D::adjacentCount(const UInteger &nodeNumber) const
{
    return node_[nodeNumber].adjacent.size();
}

AdjacentSet Mesh2D::adjacent(const UInteger &nodeNumber) const
{
    return node_[nodeNumber].adjacent;
}

Point2D Mesh2D::point2d(const UInteger &number) const
{
    return node_[number].point;
}

void Mesh2D::setPoint(const UInteger &number, const Point2D &p)
{
    node_[number].point = p;
}

void Mesh2D::evalNodalValues(std::function<double (double, double)> func)
{
    if (func != nullptr)
    {
        std::vector<double> r(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Point2D p = node_[i].point;
            r[i] = func(p.x(), p.y());
        }
        addDataVector("F(x,y)", r);
    }
}

Point2D Mesh2D::binary(Point2D p0, Point2D p1, std::function<double (double, double)> func, double level)
{
    double val0 = func(p0.x(), p0.y()) - level;
    double val1 = func(p1.x(), p1.y()) - level;

    if (0.0 <= val0 && val0 < epsilon_) return p0;
    if (0.0 <= val1 && val1 < epsilon_) return p1;

    if (signbit(val0) == signbit(val1))
    {
        return p0; // значения в узлах отрезка одного знака => нет решения
    }

    Point2D center;
    double val;
    do
    {
        center = 0.5 * (p0 + p1);
        val = func(center.x(), center.y()) - level;
        if ( signbit(val0) != signbit(val) )
        {
            p1 = center;
            val1 = val;
        }
        else
        {
            p0 = center;
            val0 = val;
        }
    } while (!(0.0 <= val && val < epsilon_));
    return center;
}

Point2D Mesh2D::findBorder(Point2D point, std::function<double (double, double)> func, const double &h, double level)
{
    if (func != nullptr)
    {
        double f = func(point.x(), point.y()) - level;
        while (fabs(f) >= epsilon_)
        {
            Point2D g = grad(func, point, 0.5*h).normalized();
            Point2D p;
            if (f < 0.0) p = point + h * g;
            else p = point - h * g;
            double fp = func(p.x(), p.y()) - level;
            if (signbit(f) != signbit(fp))
                return binary(point, p, func, level);
            else
            {
                point = p;
                f = fp;
            }
        }
    }
    return point;
}

Point2D Mesh2D::findBorder(const Point2D &a, const Point2D &b, std::function<double(double, double)> func, double alpha, double level)
{
    Point2D v(a, b);
//    double l = v.length();
    Point2D c = a + (alpha * v);
//    if (l < epsilon_)
//        return a;
    if (func == nullptr)
        return c;
    double h = a.distanceTo(c);
    return findBorder(c, func, h, level);
//    Point2D v(a, b);
//    double l = v.length();
//    Point2D n = v.perpendicular().normalized();
//    Point2D c = a + (alpha * v);
//    Point2D p0 = ((-0.25 * l) * n) + c;
//    Point2D p1 = ((0.25 * l) * n) + c;
//    Point2D border = c;
//    if (l < epsilon_)
//        return a;
//    if (func != nullptr)
//    {
//        // R-function is defined
//        if (fabs(func(a.x(), a.y()) - level) < epsilon_ && fabs(func(b.x(), b.y()) - level) < epsilon_)
//        {
//            // zero-level segment
//            if (signbit(func(p0.x(), p0.y()) - level) != signbit(func(p1.x(), p1.y()) - level))
//            {
//                border = binary(p0, p1, func);
//            }
//            else
//            {
//                Point2D h = 0.001 * l * n;
//                p0 = -h + c;
//                p1 = h + c;
////                while (signbit(func(p0.x(), p0.y()) - level) == signbit(func(p1.x(), p1.y()) - level))
////                {
////                    p0 = p0 - h;
////                    p1 = p1 + h;
////                }
//                double val0 = func(p0.x(), p0.y()) - level;
//                double val1 = func(p1.x(), p1.y()) - level;
//                short iic = 0;
//                while (signbit(val0) == signbit(val1) && fabs(val0) >= epsilon_ && fabs(val1) >= epsilon_ && iic < 5000)
//                {
//                    p0 = p0 - h;
//                    p1 = p1 + h;
//                    val0 = func(p0.x(), p0.y()) - level;
//                    val1 = func(p1.x(), p1.y()) - level;
//                    ++iic;
//                }
//                if (fabs(val0) < epsilon_) return p0;
//                if (fabs(val1) < epsilon_) return p1;
//                if (signbit(val0) == signbit(val1))
//                    return 0.5 * (a + b);
//                border = binary(p0, p1, func);
//            } // else
//        } // if zero level
//    } // if func
//    return border;
}

double Mesh2D::distToBorder(const Point2D &a, const Point2D &b, std::function<double (double, double)> func, double alpha, double level)
{
    Point2D border = findBorder(a, b, func, alpha, level);
    return border.distanceTo(a, b);
}

Point2D Mesh2D::grad(std::function<double (double, double)> func, const Point2D &p, const double &h)
{
    double x = p.x();
    double y = p.y();
    return Point2D((func(x + h, y) - func(x - h, y)) / (2.0 * h),
                   (func(x, y + h) - func(x, y - h)) / (2.0 * h));
}

void Mesh2D::laplacianSmoothing(int iter_num)
{
    std::cout << "Laplacian smoothing...";
    std::clock_t start = std::clock();
    // сглаживание Лапласа
    for (int iter = 0; iter < iter_num; iter++)
    {
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            if (node_[i].type == INNER)
            {
                Point2D nn(0.0, 0.0);
                AdjacentSet adjacent = node_[i].adjacent;
                int acount = 0;
                for (UInteger elnum: adjacent)
                {
                    ElementPointer eptr = element(elnum);
                    int index = eptr->index(i);
                    nn = nn + node_[eptr->vertexNode(index - 1)].point;
                    nn = nn + node_[eptr->vertexNode(index + 1)].point;
                    acount += 2;
                }
                node_[i].point = (1.0 /  static_cast<double>(acount)) * nn;
            }
        } // for i
        std::cout << '*';
    }
    double duration = static_cast<double>(std::clock() - start) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
}

double Mesh2D::lengthAspect(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double min = node_[e->vertexNode(0)].point.distanceTo(node_[e->vertexNode(1)].point);
    double max = min;
    for (int i = 1; i < e->verticesCount(); i++)
    {
        double d = node_[e->vertexNode(i)].point.distanceTo(node_[e->vertexNode(i + 1)].point);
        if (d < min) min = d;
        if (d > max) max = d;
    }
    return min / max;
}

double Mesh2D::faceArea(const UIntegerVector &face) const
{
    if (face.size() == 3)
    {
        Point2D A = node_[face[0]].point;
        Point2D B = node_[face[1]].point;
        Point2D C = node_[face[2]].point;
        double a = A.distanceTo(B);
        double b = B.distanceTo(C);
        double c = C.distanceTo(A);
        double p = (a + b + c) / 2.0;
        return sqrt(p * (p - a) * (p - b) * (p - c));
    }
    else if (face.size() == 4)
    {
        Point2D p0 = node_[face[0]].point;
        Point2D p1 = node_[face[1]].point;
        Point2D p2 = node_[face[2]].point;
        Point2D p3 = node_[face[3]].point;
        // стороны
        double a = p0.distanceTo(p1);
        double b = p1.distanceTo(p2);
        double c = p2.distanceTo(p3);
        double d = p3.distanceTo(p0);
        // диагонали
        double d1 = p0.distanceTo(p2);
        double d2 = p1.distanceTo(p3);
        // функция для вычисления квадрата числа (C++0x)
        auto sqr = [](double value) { return value * value; };
        return sqrt(4.0 * sqr(d1) * sqr(d2) - sqr(sqr(b) + sqr(d) - sqr(a) - sqr(c))) / 4.0;
    }
    return  0.0;
}
}
