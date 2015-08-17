#include "qtrianglemesh3d.h"

QTriangleMesh3D::QTriangleMesh3D(QObject *parent) :
    QObject(parent), TriangleMesh3D()
{
}

QTriangleMesh3D::QTriangleMesh3D(const QTriangleMesh3D &mesh) :
    QObject(mesh.parent()), TriangleMesh3D(&mesh)
{
}

QTriangleMesh3D::QTriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, QObject *parent) :
    QObject(parent), TriangleMesh3D(rCount, lCount, radius, length)
{
}

QTriangleMesh3D::QTriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, QObject *parent):
    QObject(parent), TriangleMesh3D(rCount, lCount, bottom_radius, top_radius, length)
{
}

QString QTriangleMesh3D::toString() const
{
    return tr("Сетка треугольных элементов (оболочка). Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
