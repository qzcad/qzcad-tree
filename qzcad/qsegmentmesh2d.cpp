#include "qsegmentmesh2d.h"

QSegmentMesh2D::QSegmentMesh2D(QObject *parent) : QObject(parent), SegmentMesh2D()
{    
}

QSegmentMesh2D::QSegmentMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, QObject *parent) :
    QObject(parent), SegmentMesh2D(xCount, yCount, xMin, yMin, width, height, func, charPoint)
{
}

QSegmentMesh2D::QSegmentMesh2D(const QSegmentMesh2D &qmesh) : QObject(qmesh.parent()), SegmentMesh2D(&qmesh)
{
}

QSegmentMesh2D::QSegmentMesh2D(SegmentMesh2D *mesh, QObject *parent) :
    QObject(parent), SegmentMesh2D(mesh)
{
}

QString QSegmentMesh2D::toString() const
{
    return tr("Сетка отрезков (двумерных балок-элементов). Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}

