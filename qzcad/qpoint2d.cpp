#include "qpoint2d.h"

QPoint2D::QPoint2D(QObject *parent) :
    QObject(parent), Point2D()
{
}

QPoint2D::QPoint2D(double x, double y, QObject *parent) :
    QObject(parent), Point2D(x, y)
{
}

QPoint2D::QPoint2D(const QPoint2D &point) :
    QObject(point.parent()), Point2D(point.x(), point.y())
{
}

QString QPoint2D::toString() const
{
    return QString("(%1; %2)").arg(x()).arg(y());
}
