#include "mesh3d.h"

namespace msh {
Mesh3D::Mesh3D()
{
    xMin_ = yMin_ = zMin_ = -1.0;
    xMax_ = yMax_ = zMax_ = 1.0;
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

UInteger Mesh3D::pushNode(const Point3D &point, const NodeType &type)
{
    Node3D node;
    node.point = point;
    node.type = type;
    node_.push_back(node);
    return node_.size() - 1;
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
}


