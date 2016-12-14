#include "qtrianglemesh3d.h"

QTriangleMesh3D::QTriangleMesh3D(QObject *parent) :
    QObject(parent), TriangleMesh3D()
{
}

QTriangleMesh3D::QTriangleMesh3D(const QTriangleMesh3D &mesh) :
    QObject(mesh.parent()), TriangleMesh3D(&mesh)
{
}

QTriangleMesh3D::QTriangleMesh3D(TriangleMesh3D *mesh, QObject *parent) :
    QObject(parent), TriangleMesh3D(mesh)
{
}

void QTriangleMesh3D::add(const TriangleMesh3D *mesh)
{
    TriangleMesh3D::add(mesh);
}

QString QTriangleMesh3D::toString() const
{
    return tr("Сетка треугольных элементов (оболочка). Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
