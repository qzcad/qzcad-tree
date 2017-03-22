#include "qtetrahedralmesh3d.h"

QTetrahedralMesh3D::QTetrahedralMesh3D(QObject *parent) :
    QObject(parent), TetrahedralMesh3D()
{
}

QTetrahedralMesh3D::QTetrahedralMesh3D(const QTetrahedralMesh3D &qmesh) :
    QObject(qmesh.parent()), TetrahedralMesh3D(qmesh)
{
}

QTetrahedralMesh3D::QTetrahedralMesh3D(TetrahedralMesh3D *mesh, QObject *parent) :
    QObject(parent), TetrahedralMesh3D(mesh)
{
}

QString QTetrahedralMesh3D::toString() const
{
    return tr("Сетка тетраэдрических элементов. Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
