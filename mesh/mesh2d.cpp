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

UInteger Mesh2D::addNode(const Point2D &point, const NodeType &type, double epsilon)
{
    for (UInteger i = 0; i < node_.size(); i++)
    {
        if ( point.isEqualTo(node_[i].point, epsilon) )
        {
            if (node_[i].type != type) // обновление типа узла
                node_[i].type = type;
            return i;
        }
    }
    return pushNode(point, type);
}

void Mesh2D::updateDomain()
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

UInteger Mesh2D::adjacentCount(const UInteger &nodeNumber) const
{
    return node_[nodeNumber].adjacent.size();
}
}
