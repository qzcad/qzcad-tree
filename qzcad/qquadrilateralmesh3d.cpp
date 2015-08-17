#include "qquadrilateralmesh3d.h"

QQuadrilateralMesh3D::QQuadrilateralMesh3D(QObject *parent) :
    QObject(parent), QuadrilateralMesh3D()
{
}

QQuadrilateralMesh3D::QQuadrilateralMesh3D(const QQuadrilateralMesh3D &mesh) :
    QObject(mesh.parent()), QuadrilateralMesh3D(&mesh)
{
}

QQuadrilateralMesh3D::QQuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, QObject *parent) :
    QObject(parent), QuadrilateralMesh3D(rCount, lCount, radius, length)
{
}

QQuadrilateralMesh3D::QQuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, QObject *parent) :
    QObject(parent), QuadrilateralMesh3D(rCount, lCount, bottom_radius, top_radius, length)
{
}

QString QQuadrilateralMesh3D::toString() const
{
    return tr("Сетка четырехугольных элементов (оболочка). Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
