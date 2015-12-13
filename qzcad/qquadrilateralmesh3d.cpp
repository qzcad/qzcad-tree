#include "qquadrilateralmesh3d.h"

QQuadrilateralMesh3D::QQuadrilateralMesh3D(QObject *parent) :
    QObject(parent), QuadrilateralMesh3D()
{
}

QQuadrilateralMesh3D::QQuadrilateralMesh3D(const QQuadrilateralMesh3D &mesh) :
    QObject(mesh.parent()), QuadrilateralMesh3D(&mesh)
{
}

QQuadrilateralMesh3D::QQuadrilateralMesh3D(QuadrilateralMesh3D *mesh, QObject *parent) :
    QObject(parent), QuadrilateralMesh3D(mesh)
{
}

QString QQuadrilateralMesh3D::toString() const
{
    return tr("Сетка четырехугольных элементов (оболочка). Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}

void QQuadrilateralMesh3D::add(const QQuadrilateralMesh3D *mesh)
{
    QuadrilateralMesh3D::add(mesh);
}

void QQuadrilateralMesh3D::translate(const double &x, const double &y, const double &z)
{
    QuadrilateralMesh3D::translate(x, y, z);
}
