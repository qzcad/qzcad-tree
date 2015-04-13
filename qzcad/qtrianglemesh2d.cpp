#include "qtrianglemesh2d.h"

QTriangleMesh2D::QTriangleMesh2D(QObject *parent) :
    QObject(parent), TriangleMesh2D()
{
}

QTriangleMesh2D::QTriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, QObject *parent) :
    QObject(parent), TriangleMesh2D(xCount, yCount, xMin, yMin, width, height)
{
}

QTriangleMesh2D::QTriangleMesh2D(const QTriangleMesh2D &qmesh) :
    QObject(qmesh.parent()), TriangleMesh2D(qmesh)
{
}

QString QTriangleMesh2D::toString() const
{
    return tr("Сетка треугольных элементов. Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}

QTriangleMesh2D::QTriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, QObject *parent) :
    QObject(parent), TriangleMesh2D(xCount, yCount, xMin, yMin, width, height, func)
{
}
