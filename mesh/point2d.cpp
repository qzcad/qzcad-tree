#include <math.h>
#include "point2d.h"

namespace msh
{
Point2D::Point2D() : Point1D()
{
    y_ = 0.0;
}

Point2D::Point2D(const double &x, const double &y) : Point1D(x)
{
    y_ = y;
}

Point2D::Point2D(const Point2D &point) : Point1D(point)
{
    y_ = point.y_;
}

Point2D::Point2D(const Point2D &firstPoint, const Point2D &secondPoint) : Point1D(firstPoint, secondPoint)
{
    y_ = secondPoint.y_ - firstPoint.y_;
}

int Point2D::dimension() const
{
    return 2;
}

double Point2D::y() const
{
    return y_;
}

void Point2D::set(const double &x, const double &y)
{
    Point1D::set(x);
    y_ = y;
}

void Point2D::setY(const double &y)
{
    y_ = y;
}

Point2D &Point2D::operator =(const Point2D &point)
{
    if (this != &point)
    {
        Point1D::set(point.x());
        y_ = point.y_;
    }
    return *this;
}

bool operator ==(const Point2D &leftPoint, const Point2D &rightPoint)
{
    return (leftPoint.x() == rightPoint.x()) && (leftPoint.y() == rightPoint.y());
}

const Point2D operator -(const Point2D &point)
{
    return Point2D (-point.x(), -point.y());
}

const Point2D operator -(const Point2D &leftPoint, const Point2D &rightPoint)
{
    return Point2D (rightPoint, leftPoint);
}

const Point2D operator +(const Point2D &leftPoint, const Point2D &rightPoint)
{
    return Point2D (leftPoint.x() + rightPoint.x(), leftPoint.y() + rightPoint.y());
}

double Point2D::operator *(const Point2D &point) const
{
    return this->x() * point.x() + this->y() * point.y();
}

const Point2D operator *(double dec, const Point2D &point)
{
    return Point2D (dec * point.x(), dec * point.y());
}

const Point2D operator /(const Point2D &point, double dec)
{
    return Point2D (point.x() / dec, point.y() / dec);
}

double Point2D::length() const
{
    return sqrt(x() * x() + y() * y());
}

Point2D Point2D::normalized() const
{
    double l = length();
    // для нулевого вектора направление не определено - возвращаем нулевой вектор
    if (l == 0.0)
        return Point2D(0.0, 0.0);
    return Point2D(x() / l, y() / l);
}

double Point2D::distanceTo(const Point2D &point) const
{
    double dx = point.x() - this->x();
    double dy = point.y() - this->y();
    return sqrt(dx * dx + dy * dy);
}

bool Point2D::isEqualTo(const Point2D &point, double epsilon) const
{
    return distanceTo(point) < epsilon;
}

double Point2D::product(const Point2D &point) const
{
    return x() * point.y() - point.x() * y();
}
}
