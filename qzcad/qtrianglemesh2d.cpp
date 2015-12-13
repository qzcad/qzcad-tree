#include "qtrianglemesh2d.h"

QTriangleMesh2D::QTriangleMesh2D(QObject *parent) :
    QObject(parent), TriangleMesh2D()
{
}

QTriangleMesh2D::QTriangleMesh2D(const QTriangleMesh2D &qmesh) :
    QObject(qmesh.parent()), TriangleMesh2D(qmesh)
{
}

QTriangleMesh2D::QTriangleMesh2D(TriangleMesh2D *mesh, QObject *parent) :
    QObject(parent), TriangleMesh2D(mesh)

{
}

QString QTriangleMesh2D::toString() const
{
    return tr("Сетка треугольных элементов. Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
