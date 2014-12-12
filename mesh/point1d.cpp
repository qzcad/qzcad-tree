#include "point1d.h"

namespace msh
{
Point1D::Point1D()
{
    x_ = 0.0;
}

Point1D::Point1D(const double &x)
{
    x_ = x;
}

Point1D::Point1D(const Point1D &point)
{
    x_ = point.x_;
}

Point1D::Point1D(const Point1D &firstPoint, const Point1D &secondPoint)
{
    x_ = secondPoint.x_ - firstPoint.x_;
}

int Point1D::dimension() const
{
    return 1;
}

double Point1D::x() const
{
    return x_;
}

double Point1D::y() const
{
    return 0.0;
}

double Point1D::z() const
{
    return 0.0;
}

double Point1D::length() const
{
    return x_ >= 0.0 ? x_ : -x_;
}

void Point1D::set(const double &x)
{
    x_ = x;
}

void Point1D::setX(const double &x)
{
    x_ = x;
}

Point1D &Point1D::operator =(const Point1D &point)
{
    if (this != &point)
        x_ = point.x_;
    return *this;
}

bool operator ==(const Point1D &leftPoint, const Point1D &rightPoint)
{
    return (leftPoint.x_ == rightPoint.x_);
}

Point1D operator -(const Point1D &point)
{
    return Point1D(-point.x_);
}

Point1D operator -(const Point1D &leftPoint, const Point1D &rightPoint)
{
    return Point1D(leftPoint.x_ - rightPoint.x_);
}

Point1D operator +(const Point1D &leftPoint, const Point1D &rightPoint)
{
    return Point1D(leftPoint.x_ + rightPoint.x_);
}

double Point1D::operator *(const Point1D &point) const
{
    return x_ * point.x_;
}

Point1D operator *(double d, const Point1D &point)
{
    return Point1D(d * point.x_);
}

Point1D operator /(const Point1D &point, double dec)
{
    return Point1D(point.x_ / dec);
}

double Point1D::distanceTo(const Point1D &point) const
{
    double d = x_ - point.x_;
    return d >= 0.0 ? d : -d;
}

bool Point1D::isEqualTo(const Point1D &point, double epsilon) const
{
    return distanceTo(point) < epsilon;
}


Point1D Point1D::normalized() const
{
    if (x_ > 0.0) return Point1D(1.0);
    if (x_ < 0.0) return Point1D(-1.0);
    // для нулевого вектора направление не определено - возвращаем нулевой вектор
    return Point1D(0.0);
}
}
