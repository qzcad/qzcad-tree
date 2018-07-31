#include "tetrahedralmesh3d.h"

#include <math.h>
#include <iostream>

#include "integer.h"

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
    double hz = (z1 - z0) / static_cast<double>(zLayersCount);
    double hphi = (phi1 - phi0) / static_cast<double>(zLayersCount);
    double hk = (k1 - k0) / static_cast<double>(zLayersCount);
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
    // Julien Dompierre, Paul Labbé, Marie-Gabrielle, Vallet Ricardo Camarero
    // "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra"
    const int I[6][6] = {{0, 1, 2, 3, 4, 5},
                         {1, 2, 0, 4, 5, 3},
                         {2, 0, 1, 5, 3, 4},
                         {3, 5, 4, 0, 2, 1},
                         {4, 3, 5, 1, 0, 2},
                         {5, 4, 3, 2, 1, 0}};
    for(int i = 0; i < zLayersCount; i++)
    {
        for(UInteger j = 0; j < baseMesh->elementsCount(); j++)
        {
            Triangle triangle = baseMesh->triangle(j);
            UInteger V[6];
            int r = 0;
            V[0] = triangle[0] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            V[1] = triangle[1] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            V[2] = triangle[2] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            V[3] = triangle[0] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            V[4] = triangle[1] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            V[5] = triangle[2] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            if (V[0] < V[1] && V[0] < V[2]) r = 0;
            if (V[1] < V[0] && V[1] < V[2]) r = 1;
            if (V[2] < V[0] && V[2] < V[1]) r = 2;

            if (std::min(V[I[r][1]], V[I[r][5]]) < std::min(V[I[r][2]], V[I[r][4]]))
            {
                addElement(V[I[r][0]], V[I[r][1]], V[I[r][2]], V[I[r][5]]);
                addElement(V[I[r][0]], V[I[r][1]], V[I[r][5]], V[I[r][4]]);
                addElement(V[I[r][0]], V[I[r][4]], V[I[r][5]], V[I[r][3]]);
            }
            else
            {
                addElement(V[I[r][0]], V[I[r][1]], V[I[r][2]], V[I[r][4]]);
                addElement(V[I[r][0]], V[I[r][4]], V[I[r][2]], V[I[r][5]]);
                addElement(V[I[r][0]], V[I[r][4]], V[I[r][5]], V[I[r][3]]);
            }
        }
    }
    printStats();
    updateDomain();
}

void TetrahedralMesh3D::convertHexahedralMesh(const HexahedralMesh3D *mesh)
{
    clear();
    for (UInteger in = 0; in < mesh->nodesCount(); in++)
    {
        pushNode(mesh->node3d(in));
    }
    // Julien Dompierre, Paul Labbé, Marie-Gabrielle, Vallet Ricardo Camarero
    // "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra"
    const int index[8][8] = {{0, 1, 2, 3, 4, 5, 6, 7},
                             {1, 0, 4, 5, 2, 3, 7, 6},
                             {2, 1, 5, 6, 3, 0, 4, 7},
                             {3, 0, 1, 2, 7, 4, 5, 6},
                             {4, 0, 3, 7, 5, 1, 2, 6},
                             {5, 1, 0, 4, 6, 2, 3, 7},
                             {6, 2, 1, 5, 7, 3, 0, 4},
                             {7, 3, 2, 6, 4, 0, 1, 5}};
    for (UInteger ie = 0; ie < mesh->elementsCount(); ie++)
    {
        Hexahedral V = mesh->hexahedron(ie);
        int bits = 0;
        int n = 0;
        int r = 0;
        int a = 0;
        UInteger minV;
        int I[8];
        minV = V[0];
        for (int i = 0; i < 8; i++)
        {
            if (V[i] < minV)
            {
                minV = V[i];
                r = i;
            }
        }
        for (int i = 0; i < 8; i++)
        {
            I[i] = index[r][i];
        }
        if (std::min(V[I[1]], V[I[6]]) < std::min(V[I[2]], V[I[5]]))
        {
            bits |= 4;
            n += 1;
        }
        if (std::min(V[I[3]], V[I[6]]) < std::min(V[I[2]], V[I[7]]))
        {
            bits |= 2;
            n += 1;
        }
        if (std::min(V[I[4]], V[I[6]]) < std::min(V[I[5]], V[I[7]]))
        {
            bits |= 1;
            n += 1;
        }
        if (bits == 0) a = 0;
        if (bits == 1) a = 120;
        if (bits == 2) a = 240;
        if (bits == 3) a = 0;
        if (bits == 4) a = 0;
        if (bits == 5) a = 240;
        if (bits == 6) a = 120;
        if (bits == 7) a = 0;
        if (a == 120)
        {
            int temp;
            temp = I[1];
            I[1] = I[4];
            I[4] = I[3];
            I[3] = temp;

            temp = I[5];
            I[5] = I[7];
            I[7] = I[2];
            I[2] = temp;
        }
        if (a == 240)
        {
            int temp;
            temp = I[1];
            I[1] = I[3];
            I[3] = I[4];
            I[4] = temp;

            temp = I[5];
            I[5] = I[2];
            I[2] = I[7];
            I[7] = temp;
        }
        if (n == 0)
        {
            addElement(V[I[0]], V[I[1]], V[I[2]], V[I[5]]);
            addElement(V[I[0]], V[I[2]], V[I[7]], V[I[5]]);
            addElement(V[I[0]], V[I[2]], V[I[3]], V[I[7]]);
            addElement(V[I[0]], V[I[5]], V[I[7]], V[I[4]]);
            addElement(V[I[2]], V[I[7]], V[I[5]], V[I[6]]);
        }
        if (n == 1)
        {
            addElement(V[I[0]], V[I[5]], V[I[7]], V[I[4]]);
            addElement(V[I[0]], V[I[1]], V[I[7]], V[I[5]]);
            addElement(V[I[1]], V[I[6]], V[I[7]], V[I[5]]);
            addElement(V[I[0]], V[I[7]], V[I[2]], V[I[3]]);
            addElement(V[I[0]], V[I[7]], V[I[1]], V[I[2]]);
            addElement(V[I[1]], V[I[7]], V[I[6]], V[I[2]]);
        }
        if (n == 2)
        {
            addElement(V[I[0]], V[I[4]], V[I[5]], V[I[6]]);
            addElement(V[I[0]], V[I[3]], V[I[7]], V[I[6]]);
            addElement(V[I[0]], V[I[7]], V[I[4]], V[I[6]]);
            addElement(V[I[0]], V[I[1]], V[I[2]], V[I[5]]);
            addElement(V[I[0]], V[I[3]], V[I[6]], V[I[2]]);
            addElement(V[I[0]], V[I[6]], V[I[5]], V[I[2]]);
        }
        if (n == 3)
        {
            addElement(V[I[0]], V[I[2]], V[I[3]], V[I[6]]);
            addElement(V[I[0]], V[I[3]], V[I[7]], V[I[6]]);
            addElement(V[I[0]], V[I[7]], V[I[4]], V[I[6]]);
            addElement(V[I[0]], V[I[5]], V[I[6]], V[I[4]]);
            addElement(V[I[1]], V[I[5]], V[I[6]], V[I[0]]);
            addElement(V[I[1]], V[I[6]], V[I[2]], V[I[0]]);
//            addElement(V[I[0]], V[I[1]], V[I[6]], V[I[5]]);
//            addElement(V[I[0]], V[I[1]], V[I[2]], V[I[6]]);
//            addElement(V[I[0]], V[I[2]], V[I[3]], V[I[6]]);
//            addElement(V[I[0]], V[I[5]], V[I[6]], V[I[4]]);
//            addElement(V[I[0]], V[I[6]], V[I[3]], V[I[7]]);
//            addElement(V[I[0]], V[I[6]], V[I[7]], V[I[4]]);
        }
    }
    printStats();
    updateDomain();
}

void TetrahedralMesh3D::convertHexahedralMeshFC(const HexahedralMesh3D *mesh)
{
    clear();
    for (UInteger in = 0; in < mesh->nodesCount(); in++)
    {
        pushNode(mesh->node3d(in));
    }
    for (UInteger ie = 0; ie < mesh->elementsCount(); ie++)
    {
        Hexahedral hexaheron = mesh->hexahedron(ie);
        Point3D center(0.0, 0.0, 0.0);
        UInteger center_index;
        for (int i = 0; i < 8; i++)
        {
            center = center + mesh->point3d(hexaheron[i]);
        }
        center.scale(1.0 / 8.0);
        center_index = pushNode(center, INNER);
        for (int i = 0; i < 6; i++)
        {
            UIntegerVector face = hexaheron.face(i);
            Point3D p0 = mesh->point3d(face[0]);
            Point3D p1 = mesh->point3d(face[1]);
            Point3D p2 = mesh->point3d(face[2]);
            Point3D p3 = mesh->point3d(face[3]);
            Point3D face_center = 0.25 * (p0 + p1 + p2 + p3);
            UInteger face_center_index = addNode(face_center, isBorderFace(face) == true ? BORDER : INNER);
            addElement(face[0], face[1], center_index, face_center_index);
            addElement(face[1], face[2], center_index, face_center_index);
            addElement(face[2], face[3], center_index, face_center_index);
            addElement(face[3], face[0], center_index, face_center_index);
        }
    }
    printStats();
    updateDomain();
}

void TetrahedralMesh3D::printStats() const
{
    std::cout << "Tetrahedral mesh: " << nodesCount() << " node(s), " << elementsCount() << " element(s)." << std::endl;
}

}

