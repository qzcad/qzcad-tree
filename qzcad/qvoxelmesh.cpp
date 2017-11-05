#include "qvoxelmesh.h"

QVoxelMesh::QVoxelMesh(QObject *parent) : QObject(parent), VoxelMesh()
{

}

QVoxelMesh::QVoxelMesh(const QVoxelMesh &qmesh) : QObject(qmesh.parent()), VoxelMesh(qmesh)
{

}

QVoxelMesh::QVoxelMesh(QVoxelMesh *mesh, QObject *parent) : QObject(parent), VoxelMesh(mesh)
{

}

QString QVoxelMesh::toString() const
{
    return tr("Воксельная модель. Узлов: %1; элементов: %2.").arg(nodesCount()).arg(elementsCount());
}
