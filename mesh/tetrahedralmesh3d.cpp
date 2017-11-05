#include "tetrahedralmesh3d.h"

#include <math.h>

namespace msh {

TetrahedralMesh3D::TetrahedralMesh3D() : Mesh3D(NULL)
{

}

TetrahedralMesh3D::TetrahedralMesh3D(const TetrahedralMesh3D &mesh) : Mesh3D(&mesh)
{
    element_ = mesh.element_;
}

TetrahedralMesh3D::TetrahedralMesh3D(const TetrahedralMesh3D *mesh) : Mesh3D(mesh)
{
    element_ = mesh->element_;
}

UInteger TetrahedralMesh3D::elementsCount() const
{
    return element_.size();
}

ElementPointer TetrahedralMesh3D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

double TetrahedralMesh3D::surfaceArea() const
{
    double area = 0.0;
//    for (UInteger i = 0; i < elementsCount(); i++)
//    {
//        Hexahedral hex = element_[i];
//        for (int j = 0; j < hex.facesCount(); j++)
//        {
//            UIntegerVector face = hex.face(j);
//            bool isBorderFace = true;
//            for (int k = 0; k < 4; k++)
//            {
//                if (node_[face[k]].type == INNER)
//                {
//                    isBorderFace = false;
////                    break;
//                }
//            }
//            if (isBorderFace)
//            {
//                area += faceArea(face);
//            }
//        }
//    }
    return area;
}

void TetrahedralMesh3D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3)
{
    Tetrahedron t(node0, node1, node2, node3);
    element_.push_back(t);
    // обновление списка смежных узлов
    node_[node0].adjacent.insert(element_.size() - 1);
    node_[node1].adjacent.insert(element_.size() - 1);
    node_[node2].adjacent.insert(element_.size() - 1);
    node_[node3].adjacent.insert(element_.size() - 1);
}

void TetrahedralMesh3D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0], nodes_ref[1], nodes_ref[2], nodes_ref[3]);
}

void TetrahedralMesh3D::clearElements()
{
    element_.clear();
}

void TetrahedralMesh3D::sweepBaseMesh(TriangleMesh2D *baseMesh, const double &z0, const double &z1, const double &phi0, const double &phi1, const double &k0, const double &k1, const int &zLayersCount)
{
    clear();
    double hz = (z1 - z0) / (double)zLayersCount;
    double hphi = (phi1 - phi0) / (double)zLayersCount;
    double hk = (k1 - k0) / (double)zLayersCount;
    double z = z0;
    double phi = phi0;
    double k = k0;
    for (int i = 0; i <= zLayersCount; i++)
    {
        for (UInteger j = 0; j < baseMesh->nodesCount(); j++)
        {
            Point2D plane_point = baseMesh->point2d(j);
            Point3D space_point (k * plane_point.x() * cos(phi) + k * plane_point.y() * sin(phi),
                                 k * plane_point.y() * cos(phi) - k * plane_point.x() * sin(phi),
                                 z);
            NodeType type = baseMesh->nodeType(j);
            if (i == 0 || i == zLayersCount) type = BORDER;
            pushNode(space_point, type);
        }
        z += hz;
        phi += hphi;
        k += hk;
    }
    // формирование тетраэдров
    for(int i = 0; i < zLayersCount; i++)
    {
        for(UInteger j = 0; j < baseMesh->elementsCount(); j++)
        {
            Triangle triangle = baseMesh->triangle(j);
            UInteger node1 = triangle[0] + (UInteger)i * baseMesh->nodesCount();
            UInteger node2 = triangle[1] + (UInteger)i * baseMesh->nodesCount();
            UInteger node3 = triangle[2] + (UInteger)i * baseMesh->nodesCount();
            UInteger node4 = triangle[0] + (UInteger)(i + 1) * baseMesh->nodesCount();
            UInteger node5 = triangle[1] + (UInteger)(i + 1) * baseMesh->nodesCount();
            UInteger node6 = triangle[2] + (UInteger)(i + 1) * baseMesh->nodesCount();
            Point3D p2 = node_[triangle[1]].point;
            Point3D p3 = node_[triangle[2]].point;
            if (p2.x() < p3.x() || (p2.x() == p3.y() && p2.y() < p3.y()))
            {
                addElement(node1, node2, node3, node6);
                addElement(node1, node2, node6, node5);
                addElement(node1, node5, node6, node4);
            }
            else
            {
                addElement(node1, node2, node3, node5);
                addElement(node1, node5, node3, node6);
                addElement(node1, node5, node6, node4);
            }
        }
    }
    updateDomain();
}

}

