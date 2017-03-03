#include "mesh3d.h"
#include <math.h>

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

UInteger Mesh3D::addNode(const Point3D &point, const NodeType &type)
{
    for (UInteger i = 0; i < node_.size(); i++)
    {
        if ( point.isEqualTo(node_[i].point, epsilon_) )
        {
            if (node_[i].type != type) // обновление типа узла
                node_[i].type = type;
            return i;
        }
    }
    return pushNode(point, type);
}

UInteger Mesh3D::addNode(const Node3D &node)
{
    return addNode(node.point, node.type);
}

void Mesh3D::updateDomain()
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

UInteger Mesh3D::adjacentCount(const UInteger &nodeNumber) const
{
    return node_[nodeNumber].adjacent.size();
}

Point3D Mesh3D::binary(Point3D p0, Point3D p1, std::function<double (double, double, double)> func, double level)
{
    double val0 = func(p0.x(), p0.y(), p0.z()) - level;
    double val1 = func(p1.x(), p1.y(), p1.z()) - level;

    if (0.0 <= val0 && val0 < epsilon_) return p0;
    if (0.0 <= val1 && val1 < epsilon_) return p1;

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
    } while (!(0.0 <= val && val < epsilon_) && !p0.isEqualTo(p1, epsilon_ / 10.0));
    return center;
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
}


