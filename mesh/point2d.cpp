#undef __STRICT_ANSI__
#include <math.h>
#include <float.h>
#include <iostream>
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

void Point2D::print() const
{
    std::cout << '(' << x() << "; " << y() << ')';
}

void Point2D::println() const
{
    std::cout << '(' << x() << "; " << y() << ')' << std::endl;
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
//    return (leftPoint.x() == rightPoint.x()) && (leftPoint.y() == rightPoint.y());
    return fabs(leftPoint.x() - rightPoint.x()) < DBL_EPSILON && fabs(leftPoint.y() - rightPoint.y()) < DBL_EPSILON;
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

bool operator <(const Point2D &leftPoint, const Point2D &rightPoint)
{
//    return (leftPoint.x() < rightPoint.x()) || (leftPoint.x() == rightPoint.x() && leftPoint.y() < rightPoint.y());
    return (leftPoint.x() < rightPoint.x()) || (fabs(leftPoint.x() - rightPoint.x()) < DBL_EPSILON && leftPoint.y() < rightPoint.y());
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

double Point2D::distanceTo(const Point2D &segment0, const Point2D &segment1) const
{
    Point2D v(segment0, segment1);
    Point2D w(segment0, *this);
    double c1 = w * v;
    if (c1 <= 0.0)
        return distanceTo(segment0);
    double c2 = v * v;
    if (c2 <= c1)
        return distanceTo(segment1);
    Point2D b = segment0 + (c1 / c2) * v;
    return distanceTo(b);
}

bool Point2D::isEqualTo(const Point2D &point, double epsilon) const
{
    return distanceTo(point) < epsilon;
}

double Point2D::product(const Point2D &point) const
{
    return x() * point.y() - point.x() * y();
}

Point2D::PointToSegment Point2D::classify(const Point2D &p0, const Point2D &p1)
{
    Point2D p2 = *this;
    Point2D a = p1 - p0;
    Point2D b = p2 - p0;
    double sa = a.x() * b.y() - b.x() * a.y();
    if (sa > 0.0)
        return LEFT;
    if (sa < 0.0)
        return RIGHT;
    if ((a.x() * b.x() < 0.0) || (a.y() * b.y() < 0.0))
        return BEHIND;
    if (a.length() < b.length())
        return BEYOND;
    if (p0 == p2)
        return ORIGIN;
    if (p1 == p2)
        return DESTINATION;
    return BETWEEN;
}

Point2D Point2D::perpendicular() const
{
    return Point2D(-y(), x());
}

void Point2D::scale(const double &d)
{
    Point1D::scale(d);
    y_ *= d;
}

double Point2D::angle(const Point2D &B, const Point2D &C) const
{
    Point2D AB(*this, B);
    Point2D AC(*this, C);
    double ab = AB.length();
    double ac = AC.length();
    if (fabs(ab) < DBL_EPSILON || fabs(ac) < DBL_EPSILON)
        return 0.0;
    if (AC.product(AB) < 0.0)
        return 2.0 * M_PI - acos((AB * AC) / ab / ac);
    return acos((AB * AC) / ab / ac);
}

bool isCrossed(const Point2D &P0, const Point2D &P1, const Point2D &Q0, const Point2D &Q1, double &p, double &q)
{
    double x[2][2], y[2][2];
    double A[2], B[2], C[2];
    double det;
    x[0][0] = P0.x();
    y[0][0] = P0.y();
    x[1][0] = P1.x();
    y[1][0] = P1.y();
    x[0][1] = Q0.x();
    y[0][1] = Q0.y();
    x[1][1] = Q1.x();
    y[1][1] = Q1.y();
    A[0] = x[1][0] - x[0][0];
    A[1] = y[1][0] - y[0][0];
    B[0] = x[0][1] - x[1][1];
    B[1] = y[0][1] - y[1][1];
    C[0] = x[0][1] - x[0][0];
    C[1] = y[0][1] - y[0][0];
    det = A[0] * B[1] - B[0] * A[1];
    if (fabs(det) < 1.0E-10)
        return false;
    p = (C[0] * B[1] - B[0] * C[1]) / det;
    q = (A[0] * C[1] - C[0] * A[1]) / det;
    return true;
}
}
