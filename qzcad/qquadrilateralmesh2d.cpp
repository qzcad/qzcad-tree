#include "qquadrilateralmesh2d.h"

QQuadrilateralMesh2D::QQuadrilateralMesh2D(QObject *parent) :
    QObject(parent), QuadrilateralMesh2D()
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const QQuadrilateralMesh2D &mesh) :
    QObject(mesh.parent()), QuadrilateralMesh2D(&mesh)
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(QuadrilateralMesh2D *mesh, QObject *parent):
    QObject(parent), QuadrilateralMesh2D(mesh)
{
}

QString QQuadrilateralMesh2D::toString() const
{
    return tr("Сетка четырехугольных элементов. Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
