#include "qquadrilateralmesh2d.h"

QQuadrilateralMesh2D::QQuadrilateralMesh2D(QObject *parent) :
    QObject(parent), QuadrilateralMesh2D()
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, QObject *parent) :
QObject(parent), QuadrilateralMesh2D(xCount, yCount, xMin, yMin, width, height)
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const Point2D &v0, const Point2D &v1, const Point2D &v2, const Point2D &v3, QObject *parent) :
    QObject(parent), QuadrilateralMesh2D(xCount, yCount, v0, v1, v2, v3)
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2, QObject *parent) :
    QObject(parent), QuadrilateralMesh2D(count, v0, v1, v2)
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const UInteger &count, const Point2D &center, const double &radius, unsigned short part, QObject *parent) :
    QObject(parent), QuadrilateralMesh2D(count, center, radius, part)
{

}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const QQuadrilateralMesh2D &mesh) :
    QObject(mesh.parent()), QuadrilateralMesh2D(&mesh)
{
}

QQuadrilateralMesh2D::QQuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, QObject *parent) :
    QObject(parent), QuadrilateralMesh2D(xCount, yCount, xMin, yMin, width, height, func, charPoint)
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
