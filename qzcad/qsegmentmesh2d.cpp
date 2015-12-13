#include "qsegmentmesh2d.h"

QSegmentMesh2D::QSegmentMesh2D(QObject *parent) : QObject(parent), SegmentMesh2D()
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

