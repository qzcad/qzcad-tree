#include "mesh3d.h"
#include <math.h>
#include "point3d.h"

namespace msh {
Mesh3D::Mesh3D(const Mesh3D *mesh) : Mesh(mesh)
{
    if (mesh != NULL)
    {
        xMin_ = mesh->xMin_;
        yMin_ = mesh->yMin_;
        zMin_ = mesh->zMin_;
        xMax_ = mesh->xMax_;
        yMax_ = mesh->yMax_;
        zMax_ = mesh->zMax_;
        node_ = mesh->node_;
    }
    else
    {
        xMin_ = yMin_ = zMin_ = -1.0;
        xMax_ = yMax_ = zMax_ = 1.0;
    }
}

UInteger Mesh3D::nodesCount() const
{
    return node_.size();
}

PointPointer Mesh3D::node(const UInteger &number) const
{
    PointPointer pointPtr = &node_[number].point;
    return pointPtr;
}

NodeType Mesh3D::nodeType(const UInteger &number) const
{
    return node_[number].type;
}

double Mesh3D::xMin() const
{
    return xMin_;
}

double Mesh3D::xMax() const
{
    return xMax_;
}

double Mesh3D::yMin() const
{
    return yMin_;
}

double Mesh3D::yMax() const
{
    return yMax_;
}

double Mesh3D::zMin() const
{
    return zMin_;
}

double Mesh3D::zMax() const
{
    return zMax_;
}

int Mesh3D::dimesion() const
{
    return 3;
}

void Mesh3D::clearNodes()
{
    node_.clear();
}

UInteger Mesh3D::pushNode(const Point3D &point, const NodeType &type)
{
    Node3D node;
    node.point = point;
    node.type = type;
    node_.push_back(node);
    return node_.size() - 1;
}

UInteger Mesh3D::pushNode(PointPointer point, const NodeType &type)
{
    return pushNode(Point3D(point->x(), point->y(), point->z()), type);
}

UInteger Mesh3D::addNode(const Point3D &point, const NodeType &type, double epsilon)
{
    UInteger ns = node_.size();
    for (UInteger i = 0; i < ns; i++)
    {
        UInteger ii = ns - i - 1UL;
        if ( point.isEqualTo(node_[ii].point, epsilon) )
        {
            if (node_[ii].type != type) // обновление типа узла
                node_[ii].type = type;
            return ii;
        }
    }
    return pushNode(point, type);
}

UInteger Mesh3D::addNode(const Node3D &node, double epsilon)
{
    return addNode(node.point, node.type, epsilon);
}

Point3D Mesh3D::normal(const UIntegerVector &face) const
{
    int fs = face.size();
    if (fs == 3)
    {
        return normal3(node_[face[0]].point, node_[face[1]].point, node_[face[2]].point);
    }
    // else we use average of normals
    Point3D n(0.0, 0.0, 0.0);
    for (int i = 0; i < fs; i++)
    {
        Point3D prev = (i == 0) ? node_[face[fs - 1]].point : node_[face[i - 1]].point;
        Point3D curr = node_[face[i]].point;
        Point3D next = (i == fs - 1) ? node_[face[0]].point : node_[face[i + 1]].point;
        n = n + normal3(prev, curr, next);
    }
    n.scale(1.0 / (double)fs);
    return n.normalized();
}

void Mesh3D::updateDomain()
{
    if (nodesCount() > 0)
    {
        Point3D start = node_[0].point;
        xMin_ = xMax_ = start.x();
        yMin_ = yMax_ = start.y();
        zMin_ = zMax_ = start.z();
        for (UInteger i = 1; i < node_.size(); i++)
        {
            Point3D current = node_[i].point;
            double x = current.x();
            double y = current.y();
            double z = current.z();
            if (xMin_ > x) xMin_ = x;
            if (xMax_ < x) xMax_ = x;
            if (yMin_ > y) yMin_ = y;
            if (yMax_ < y) yMax_ = y;
            if (zMin_ > z) zMin_ = z;
            if (zMax_ < z) zMax_ = z;
        }
    }
}

AdjacentSet Mesh3D::adjacent(const UInteger &nodeNumber) const
{
    return node_[nodeNumber].adjacent;
}

UInteger Mesh3D::adjacentCount(const UInteger &nodeNumber) const
{
    return node_[nodeNumber].adjacent.size();
}

Point3D Mesh3D::binary(Point3D p0, Point3D p1, std::function<double (double, double, double)> func, double level)
{
    double val0 = func(p0.x(), p0.y(), p0.z()) - level;
    double val1 = func(p1.x(), p1.y(), p1.z()) - level;

    if (fabs(val0) < epsilon_) return p0;
    if (fabs(val1) < epsilon_) return p1;

    if (signbit(val0) == signbit(val1))
    {
        return Point3D(); // значения в узлах отрезка одного знака => нет решения
    }

    Point3D center;
    double val;
    do
    {
        center = 0.5 * (p0 + p1);
        val = func(center.x(), center.y(), center.z()) - level;
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
    } while (!(fabs(val) < epsilon_));
    return center;
}

Point3D Mesh3D::findBorder(Point3D point, std::function<double (double, double, double)> func, double h, double level)
{
    if (func != nullptr)
    {
        double f = func(point.x(), point.y(), point.z()) - level;
        const double dh = 0.5 * h;
        while (fabs(f) >= epsilon_)
        {
            Point3D g = grad(func, point, dh).normalized();
            Point3D p;
            if (f < 0.0) p = point + h * g;
            else p = point - h * g;
            double fp = func(p.x(), p.y(), p.z()) - level;

            if (signbit(f) != signbit(fp))
                return binary(point, p, func, level);

            point = p;
            f = fp;
        }
    }
    return point;
}

Point3D Mesh3D::findBorder(const Point3D &point, const Point3D &normal, std::function<double (double, double, double)> func, const double &h, double level)
{
    if (normal.length() < epsilon_) return point;

    Point3D p0 = point;
    double val0 = func(p0.x(), p0.y(), p0.z()) - level;
    double sign = (val0 > 0.0) ? 1.0 : -1.0;
    Point3D p1 =  point + sign * h * normal;
    double val1 = func(p1.x(), p1.y(), p1.z()) - level;
    short iic = 0;
    double step = h * 0.01;
    while (signbit(val0) == signbit(val1)
           && fabs(val0) >= epsilon_ && fabs(val1) >= epsilon_ && iic < 10000)
    {
        p1 = p0 + sign * step * normal;
        val1 = func(p1.x(), p1.y(), p1.z()) - level;
        step = step + h * 0.01;
        ++iic;
    }
    if (fabs(val0) < epsilon_) return p0;
    if (fabs(val1) < epsilon_) return p1;
    if (signbit(func(p0.x(), p0.y(), p0.z()) - level) == signbit(func(p1.x(), p1.y(), p1.z()) - level))
        return point;
    return binary(p0, p1, func, level);
}

Point3D Mesh3D::findBorder(const Point3D &a, const Point3D &b, const Point3D &c, std::function<double (double, double, double)> func, double lb, double lc, double level)
{
    Point3D point = ((1.0 - lb - lc) * a) + (lb * b) + (lc * c);
    Point3D center = ((1.0 - 0.333 - 0.333) * a) + (0.333 * b) + (0.333 * c);
//    Point3D n = normal3(a, b, c);
//    double h = a.distanceTo(center);
//    Point3D p0 = point - h * n;
//    Point3D p1 = point + h * n;
//    if (fabs(func(point.x(), point.y(), point.z()) - level) < epsilon_)
//        return point;
//    if (signbit(func(p0.x(), p0.y(), p0.z()) - level) == signbit(func(p1.x(), p1.y(), p1.z()) - level))
//        return a;
//    return binary(p0, p1, func, level);
//    if (a.isEqualTo(center, epsilon_) || b.isEqualTo(center, epsilon_) || c.isEqualTo(center, epsilon_))
//        return a;
    return findBorder(point, func, a.distanceTo(center), level);

//    Point3D ab(a, b);
//    Point3D ac(a, c);
//    double lab = ab.length();
//    double lac = ac.length();
//    if (lab < epsilon_ || lac < epsilon_) return a;
//    Point3D normal = ab.product(ac);
//    if (normal.length() < epsilon_) return a;

//    Point3D point = ((1.0 - lb - lc) * a) + (lb * b) + (lc * c);
//    if (fabs(func(point.x(), point.y(), point.z()) - level) < epsilon_) return point;

//    double avr_len = 0.5 * (lab + lac);
//    Point3D h = 0.01 * avr_len * normal.normalized();
//    Point3D p0 = point - h;
//    Point3D p1 = point + h;
//    double val0 = func(p0.x(), p0.y(), p0.z()) - level;
//    double val1 = func(p1.x(), p1.y(), p1.z()) - level;
//    short iic = 0;
//    while (signbit(val0) == signbit(val1)
//           && fabs(val0) >= epsilon_ && fabs(val1) >= epsilon_ && iic < 10000)
//    {
//        p0 = p0 - h;
//        p1 = p1 + h;
//        val0 = func(p0.x(), p0.y(), p0.z()) - level;
//        val1 = func(p1.x(), p1.y(), p1.z()) - level;
//        ++iic;
//    }
//    if (fabs(val0) < epsilon_) return p0;
//    if (fabs(val1) < epsilon_) return p1;
//    if (signbit(func(p0.x(), p0.y(), p0.z()) - level) == signbit(func(p1.x(), p1.y(), p1.z()) - level))
//        return a;
    //    return binary(p0, p1, func, level);
}

Point3D Mesh3D::findBorder(const Point3D &a, const Point3D &b, const Point3D &c, std::function<double (double, double, double)> func, Point3D n, double lb, double lc, double level)
{
    Point3D point = ((1.0 - lb - lc) * a) + (lb * b) + (lc * c);
    Point3D center = ((1.0 - 0.333 - 0.333) * a) + (0.333 * b) + (0.333 * c);
    double h = a.distanceTo(center);
    Point3D p0 = point - h * n;
    Point3D p1 = point + h * n;
    if (fabs(func(point.x(), point.y(), point.z()) - level) < epsilon_)
        return point;
    if (signbit(func(p0.x(), p0.y(), p0.z()) - level) == signbit(func(p1.x(), p1.y(), p1.z()) - level))
        return a;
    return binary(p0, p1, func, level);
}

double Mesh3D::distToBorder(const Point3D &a, const Point3D &b, const Point3D &c, std::function<double (double, double, double)> func, double lb, double lc, double level)
{
    Point3D border = findBorder(a, b, c, func, lb, lc, level);
    return border.distanceTo(((1.0 - lb - lc) * a) + (lb * b) + (lc * c));
}

Point3D Mesh3D::grad(std::function<double (double, double, double)> func, const Point3D &p, const double &h)
{
    double x = p.x();
    double y = p.y();
    double z = p.z();
    return Point3D((func(x + h, y, z) - func(x - h, y, z)) / (2.0 * h),
                   (func(x, y + h, z) - func(x, y - h, z)) / (2.0 * h),
                   (func(x, y, z + h) - func(x, y, z - h)) / (2.0 * h));
//    return Point3D((func(x - h - h, y, z) - 8.0 * func(x - h, y, z) + 8.0 * func(x + h, y, z) - func(x + h + h, y, z)) / (12.0 * h),
//                   (func(x, y - h - h, z) - 8.0 * func(x, y - h, z) + 8.0 * func(x, y + h, z) - func(x, y + h + h, z)) / (12.0 * h),
//                   (func(x, y, z - h - h) - 8.0 * func(x, y, z - h) + 8.0 * func(x, y, z + h) - func(x, y, z + h + h)) / (12.0 * h));
}


void Mesh3D::evalNodalValues(std::function<double (double, double, double)> func)
{
    if (func != nullptr)
    {
        std::vector<double> r(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Point3D p = node_[i].point;
            r[i] = func(p.x(), p.y(), p.z());
        }
        addDataVector("F(x,y,z)", r);
    }
}

Point3D Mesh3D::point3d(const UInteger &nnumber) const
{
    return node_[nnumber].point;
}
}


