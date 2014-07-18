#include "mesh2d.h"
#include <math.h>
namespace msh
{
Mesh2D::Mesh2D()
{
    xMin_ = yMin_ = -1.0;
    xMax_ = yMax_ = 1.0;
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

Floating Mesh2D::xMin() const
{
    return xMin_;
}

Floating Mesh2D::xMax() const
{
    return xMax_;
}

Floating Mesh2D::yMin() const
{
    return yMin_;
}

Floating Mesh2D::yMax() const
{
    return yMax_;
}

Floating Mesh2D::zMin() const
{
    return -std::min(std::max(fabs(xMax_), fabs(xMin_)), std::max(fabs(yMax_), fabs(yMin_)));
}

Floating Mesh2D::zMax() const
{
    return std::min(std::max(fabs(xMax_), fabs(xMin_)), std::max(fabs(yMax_), fabs(yMin_)));
}

int Mesh2D::dimesion() const
{
    return 2;
}

Floating Mesh2D::elementValue(const UInteger &number) const
{
    return (Floating)number;
}

Floating Mesh2D::nodeValue(const UInteger &number) const
{
    return node_[number].point.length();
}

NodeType Mesh2D::nodeType(const UInteger &number) const
{
    return node_[number].type;
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

Floating Mesh2D::area()
{
    Floating sum = 0.0;
    for (UInteger i = 0; i < elementsCount(); i++)
        sum += area(i);
    return sum;
}

UInteger Mesh2D:: pushNode(const Point2D &point, const NodeType &type)
{
    Node2D node;
    node.point = point;
    node.type = type;
    node_.push_back(node);
    return node_.size() - 1;
}

UInteger Mesh2D::addNode(const Point2D &point, const NodeType &type)
{
    for (UInteger i = 0; i < node_.size(); i++)
    {
        if ( point.isEqualTo(node_[i].point) )
        {
            if (node_[i].type != type) // обновление типа узла
                node_[i].type = type;
            return i;
        }
    }
    return pushNode(point, type);
}
}
