#include "point3d.h"
#include <math.h>
#include <iostream>

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

double Point3D::distanceTo(const Point3D &segment0, const Point3D &segment1) const
{
    Point3D v(segment0, segment1);
    Point3D w(segment0, *this);
    double c1 = w * v;
    if (c1 <= 0.0)
        return distanceTo(segment0);
    double c2 = v * v;
    if (c2 <= c1)
        return distanceTo(segment1);
    Point3D b = segment0 + (c1 / c2) * v;
    return distanceTo(b);
}

double Point3D::distanceTo(const Point3D &triangle0, const Point3D &triangle1, const Point3D &triangle2) const
{
    Point3D edge0(triangle0, triangle1); // edge0 = triangle1 - triangle0
    Point3D edge1(triangle1, triangle2); // edge1 = triangle2 - triangle1
    Point3D edge2(triangle2, triangle0); // edge2 = triangle0 - triangle2
    Point3D n = edge0.product(-edge2); // the normal to the triangle
    Point3D q(triangle0, *this); // q = this - triangle0
    double l = (q * n) / n.length(); // the distance to the plane of the tiangle
    Point3D b(l * n.normalized(), *this); // this - l * n / |n|
    Point3D c0(triangle0, b); // q = b - triangle0
    Point3D c1(triangle1, b); // q = b - triangle1
    Point3D c2(triangle2, b); // q = b - triangle2
    double v0 = (n * edge0.product(c0));
    double v1 = (n * edge1.product(c1));
    double v2 = (n * edge2.product(c2));
    if (v0 >= 0.0 && v1 >= 0.0 && v2 >= 0.0)
        return fabs(l); // if b in the triangle
    if (v2 < 0 && v0 < 0)
        return distanceTo(triangle0);
    if (v0 < 0 && v1 < 0)
        return distanceTo(triangle1);
    if (v1 < 0 && v2 < 0)
        return distanceTo(triangle2);
    if (v0 < 0)
        return distanceTo(triangle0, triangle1);
    if (v1 < 0)
        return distanceTo(triangle1, triangle2);
    if (v2 < 0)
        return distanceTo(triangle2, triangle0);
    // otherwise (absolutely theoretical now)
    return std::min(distanceTo(triangle0, triangle1), std::min(distanceTo(triangle1, triangle2), distanceTo(triangle2, triangle0)));
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

void Point3D::print() const
{
    std::cout << '(' << x() << "; " << y() << "; " << z() << ')';
}

void Point3D::println() const
{
    std::cout << '(' << x() << "; " << y() << "; " << z() << ')' << std::endl;
}

Point3D Point3D::inCoordSystem(const Point3D &Vx, const Point3D &Vy, const Point3D &Vz) const
{
    Point3D p(Vx.x() * x() + Vx.y() * y() + Vx.z() * z(),
              Vy.x() * x() + Vy.y() * y() + Vy.z() * z(),
              Vz.x() * x() + Vz.y() * y() + Vz.z() * z());
    return p;
}

void Point3D::scale(const double &d)
{
    Point2D::scale(d);
    z_ *= d;
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

bool operator <(const Point3D &leftPoint, const Point3D &rightPoint)
{
    return  (leftPoint.x() < rightPoint.x()) ||
            (leftPoint.x() == rightPoint.x() && leftPoint.y() < rightPoint.y()) ||
            (leftPoint.x() == rightPoint.x() && leftPoint.y() == rightPoint.y() && leftPoint.z() < rightPoint.z());
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

bool isSkew(const Point3D &P0, const Point3D &P1, const Point3D &Q0, const Point3D &Q1, double &p, double &q)
{
    Point3D AB = P1 - P0;
    Point3D AC = Q1 - P0;
    Point3D N = AB.product(AC);
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vz.product(Vx);

    Point3D p0 = (P0 - P0).inCoordSystem(Vx, Vy, Vz);
    Point3D p1 = (P1 - P0).inCoordSystem(Vx, Vy, Vz);
    Point3D q0 = (Q0 - P0).inCoordSystem(Vx, Vy, Vz);
    Point3D q1 = (Q1 - P0).inCoordSystem(Vx, Vy, Vz);

    return isCrossed(p0, p1, q0, q1, p, q);
}

}
