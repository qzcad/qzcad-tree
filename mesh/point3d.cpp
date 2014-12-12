#include "point3d.h"
#include <math.h>

namespace msh
{
Point3D::Point3D() : Point2D()
{
    z_ = 0.0;
}

Point3D::Point3D(const double &x, const double &y, const double &z) : Point2D(x, y)
{
    z_ = z;
}

Point3D::Point3D(const Point3D &point) : Point2D(point)
{
    z_ = point.z_;
}

Point3D::Point3D(const Point3D &firstPoint, const Point3D &secondPoint) : Point2D(firstPoint, secondPoint)
{
    z_ = secondPoint.z_ - firstPoint.z_;
}

int Point3D::dimension() const
{
    return 3;
}

double Point3D::z() const
{
    return z_;
}

double Point3D::length() const
{
    return sqrt(x() * x() + y() * y() + z() * z());
}

Point3D Point3D::normalized() const
{
    double l = length();
    // для нулевого вектора направление не определено - возвращаем нулевой вектор
    if (l == 0.0)
        return Point3D(0.0, 0.0, 0.0);
    return Point3D(x() / l, y() / l, z() / l);
}

double Point3D::distanceTo(const Point3D &point) const
{
    double dx = point.x() - x();
    double dy = point.y() - y();
    double dz = point.z() - z();
    return sqrt(dx * dx + dy * dy + dz * dz);
}

bool Point3D::isEqualTo(const Point3D &point, double epsilon) const
{
    return distanceTo(point) < epsilon;
}

Point3D Point3D::product(const Point3D &point) const
{
    return Point3D(y() * point.z() - z() * point.y(),
                   z() * point.x() - x() * point.z(),
                   x() * point.y() - y() * point.x());
}

void Point3D::set(const double &x, const double &y, const double &z)
{
    Point2D::set(x, y);
    z_ = z;
}

void Point3D::setZ(const double &z)
{
    z_ = z;
}

Point2D &Point3D::operator =(const Point3D &point)
{
    if (this != &point)
    {
        set(point.x(), point.y(), point.z());
    }
    return *this;
}

double Point3D::operator *(const Point3D &point) const
{
    return x() * point.x() + y() * point.y() + z() * point.z();
}

const Point3D operator /(const Point3D &point, double dec)
{
    return Point3D(point.x() / dec, point.y() / dec, point.z() / dec);
}

const Point3D operator *(double dec, const Point3D &point)
{
    return Point3D(dec * point.x(), dec * point.y(), dec * point.z());
}

const Point3D operator +(const Point3D &leftPoint, const Point3D &rightPoint)
{
    return Point3D(leftPoint.x() + rightPoint.x(), leftPoint.y() + rightPoint.y(), leftPoint.z() + rightPoint.z());
}

const Point3D operator -(const Point3D &leftPoint, const Point3D &rightPoint)
{
    return Point3D(rightPoint, leftPoint); // left - right
}

const Point3D operator -(const Point3D &point)
{
    return Point3D (-point.x(), -point.y(), -point.z());
}

bool operator ==(const Point3D &leftPoint, const Point3D &rightPoint)
{
    return (leftPoint.x() == rightPoint.x()) && (leftPoint.y() == rightPoint.y()) && (leftPoint.z() == rightPoint.z());
}
}
