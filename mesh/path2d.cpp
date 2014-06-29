#include "path2d.h"
#include "segment.h"

namespace msh
{
Path2D::Path2D()
{
    xMax_ = yMax_ = 1.0;
    xMin_ = yMin_ = -1.0;
}

Path2D::Path2D(const Path2D &path2d)
{
    xMin_ = path2d.xMin_;
    xMax_ = path2d.xMax_;
    yMin_ = path2d.yMin_;
    yMax_ = path2d.yMax_;
    node_ = path2d.node_;
}

Floating Path2D::yMax() const
{
    return yMax_;
}

void Path2D::setYMax(const Floating &yMax)
{
    yMax_ = yMax;
}

UInteger Path2D::nodesCount() const
{
    return node_.size();
}

PointPointer Path2D::node(UInteger number) const
{
    PointPointer pointPointer(new Point2D(node_[number]));
    return pointPointer;
}

UInteger Path2D::elementsCount() const
{
    if (node_.size() > 1)
        return node_.size();
    return 0;
}

ElementPointer Path2D::element(UInteger number) const
{
//    Segment *segmentPtr;
//    if (number < node_.size() - 1)
//    {
//        segmentPtr = new Segment(number, number + 1);
//    }
//    else
//    {
//        segmentPtr = new Segment(number, 0);
//    }
//    ElementPointer elementPointer (segmentPtr);
//    return elementPointer;
}

Floating Path2D::zMin() const
{
    return -1.0;
}

Floating Path2D::zMax() const
{
    return 1.0;
}

int Path2D::dimesion() const
{
    return 2;
}

void Path2D::addNode(const Point2D &point)
{
    node_.push_back(point);
    if (xMin_ > point.x()) xMin_ = point.x();
    if (xMax_ < point.x()) xMax_ = point.x();
    if (yMin_ > point.y()) yMin_ = point.y();
    if (yMax_ < point.y()) yMax_ = point.y();
}

void Path2D::removeNode(UInteger number)
{
    node_.erase(node_.begin() + number);
}

void Path2D::insertNode(UInteger number, const Point2D &point)
{
    node_.insert(node_.begin() + number, point);
}

Floating Path2D::yMin() const
{
    return yMin_;
}

void Path2D::setYMin(const Floating &yMin)
{
    yMin_ = yMin;
}

Floating Path2D::xMax() const
{
    return xMax_;
}

void Path2D::setXMax(const Floating &xMax)
{
    xMax_ = xMax;
}

Floating Path2D::xMin() const
{
    return xMin_;
}

void Path2D::setXMin(const Floating &xMin)
{
    xMin_ = xMin;
}
}
