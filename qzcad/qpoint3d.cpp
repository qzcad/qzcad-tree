#include "qpoint3d.h"

QPoint3D::QPoint3D(QObject *parent) :
    QObject(parent), Point3D()
{
}

QPoint3D::QPoint3D(const double &x, const double &y, const double &z, QObject *parent) :
    QObject(parent), Point3D(x, y, z)
{
}

QPoint3D::QPoint3D(const QPoint3D &point) :
    QObject(point.parent()), Point3D(point.x(), point.y(), point.z())
{
}

QString QPoint3D::toString() const
{
    return QString("(%1; %2; %3)").arg(x()).arg(y()).arg(z());
}
