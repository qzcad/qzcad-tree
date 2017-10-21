#include "qhexahedralmesh3d.h"

QHexahedralMesh3D::QHexahedralMesh3D(QObject *parent) : QObject(parent), HexahedralMesh3D()
{

}

QHexahedralMesh3D::QHexahedralMesh3D(const QHexahedralMesh3D &qmesh) : QObject(qmesh.parent()), HexahedralMesh3D(qmesh)
{

}

QHexahedralMesh3D::QHexahedralMesh3D(HexahedralMesh3D *mesh, QObject *parent) : QObject(parent), HexahedralMesh3D(mesh)
{

}

QString QHexahedralMesh3D::toString() const
{
    return tr("Сетка шестигранных элементов. Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
