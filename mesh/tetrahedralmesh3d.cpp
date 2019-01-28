#include "tetrahedralmesh3d.h"

#include <math.h>
#include <iostream>
#include <map>
#include <ctime>

#include "integer.h"
#include "funcopt.h"
#include "consoleprogress.h"

namespace msh {

TetrahedralMesh3D::TetrahedralMesh3D() : Mesh3D(nullptr)
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

void TetrahedralMesh3D::addElement(const Tetrahedron &t)
{
    element_.push_back(t);
    // обновление списка смежных узлов
    node_[t[0]].adjacent.insert(element_.size() - 1);
    node_[t[1]].adjacent.insert(element_.size() - 1);
    node_[t[2]].adjacent.insert(element_.size() - 1);
    node_[t[3]].adjacent.insert(element_.size() - 1);
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

void TetrahedralMesh3D::addElement(const UInteger &i0, const UInteger &i1, const UInteger &i2, const UInteger &i3, const UInteger &i4, const UInteger &i5)
{
    UInteger V[6];
    // Julien Dompierre, Paul Labbé, Marie-Gabrielle, Vallet Ricardo Camarero
    // "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra"
    const int I[6][6] = {{0, 1, 2, 3, 4, 5},
                         {1, 2, 0, 4, 5, 3},
                         {2, 0, 1, 5, 3, 4},
                         {3, 5, 4, 0, 2, 1},
                         {4, 3, 5, 1, 0, 2},
                         {5, 4, 3, 2, 1, 0}};
    int r = 0;
    V[0] = i0;
    V[1] = i1;
    V[2] = i2;
    V[3] = i3;
    V[4] = i4;
    V[5] = i5;
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

void TetrahedralMesh3D::delaunay(const TriangleMesh3D &mesh, std::function<double (double, double, double)> func)
{
    clear();
    TriangleMesh3D tmesh = mesh;
    std::cout << "Delaunay Meshing: " << tmesh.nodesCount() << " nodes, " << tmesh.elementsCount() << " elements in the initial mesh." << std::endl;
    if (func == nullptr) func = std::bind(&TriangleMesh3D::cfunction, &tmesh, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
}

void TetrahedralMesh3D::backgroundGrid(const TetrahedralMesh3D *mesh, std::function<double (double, double, double)> func, double level, int smooth, int optimize)
{
    clear();
    TriangleMesh3D tmesh;
    std::list<ElementPointer> inner = tmesh.backgroundGrid(mesh, func, level, smooth, optimize, false);
    std::map<UInteger, UInteger> iso;
    for (ElementPointer el: inner)
    {
        for (int j = 0; j < el->facesCount(); j++)
        {
            UIntegerVector f = el->face(j);
            Point3D p[3];
            p[0] = mesh->point3d(f[0]);
            p[1] = mesh->point3d(f[1]);
            p[2] = mesh->point3d(f[2]);
            bool isInner = false;
            // check an other element with these vertices
            for (ElementPointer el1: inner)
            {
                if (el1 != el && el1->in(f[0]) && el1->in(f[1]) && el1->in(f[2]))
                {
                    isInner = true;
                    break;
                }
            }
            if (!isInner)
            {
                iso[f[0]] = addNode(p[0], INNER);
                iso[f[1]] = addNode(p[1], INNER);
                iso[f[2]] = addNode(p[2], INNER);
            }
        }
    }
    if (nodesCount() == tmesh.nodesCount())
    {
        std::cout << "Nodes count is equal." << std::endl;
    }
    else
    {
        std::cout << "Nodes count is not equal!!!" << std::endl;
        return;
    }
    for (ElementPointer el: inner)
    {
        std::vector<UInteger> v(4);
        for (int i = 0; i < 4; i++)
        {
            v[i] = addNode(mesh->point3d(el->vertexNode(i)), INNER);
        }
        addElement(v);
    }
    UInteger shift = nodesCount();
    for (UInteger i = 0; i < tmesh.nodesCount(); i++)
    {
        pushNode(tmesh.node3d(i));
    }
    for (ElementPointer el: inner)
    {
        for (int j = 0; j < el->facesCount(); j++)
        {
            UIntegerVector f = el->face(j);
            Point3D p[3];
            p[0] = mesh->point3d(f[0]);
            p[1] = mesh->point3d(f[1]);
            p[2] = mesh->point3d(f[2]);
            bool isInner = false;
            // check an other element with these vertices
            for (ElementPointer el1: inner)
            {
                if (el1 != el && el1->in(f[0]) && el1->in(f[1]) && el1->in(f[2]))
                {
                    isInner = true;
                    break;
                }
            }
            if (!isInner)
            {
                addElement(iso[f[0]], iso[f[1]], iso[f[2]],
                        iso[f[0]] + shift, iso[f[1]] + shift, iso[f[2]] + shift);
            }
        }
    }
    laplacianSmoothing(smooth);
//    localFubctionalOptimization(optimize);
//    std::list<UInteger> ee;
//    ee.push_back(25943); ee.push_back(25942); ee.push_back(25941);
//    subdivide(ee, func);
    printStats();
    updateDomain();
}

void TetrahedralMesh3D::localFubctionalOptimization(int maxiter)
{
    std::cout << "Local optimization smoothing...";
    std::clock_t start = std::clock();
    // проект сглаживания путем локальной минимизации функционала
    for (int iii = 0; iii < maxiter; iii++)
    {
        ConsoleProgress progress_bar(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            if (node_[i].type == INNER)
            {
                Point3D point = node_[i].point;
                std::vector<double> x0(3);
                std::vector<double> x(3);
                AdjacentSet adjacent = node_[i].adjacent;
                double avr_dist = 0.0;

                for (UInteger elnum: adjacent)
                {
                    Tetrahedron el = element_[elnum];
                    int index = el.index(i);
                    Point3D p0 = node_[el[index + 1]].point;
                    Point3D p1 = node_[el[index + 2]].point;
                    Point3D p2 = node_[el[index + 3]].point;
                    avr_dist += point.distanceTo(p0) + point.distanceTo(p1) + point.distanceTo(p2);
                }

                avr_dist /= static_cast<double>(3 * adjacent.size());

                auto functor = [&](const std::vector<double> &vars){
                    double f = 0.0;
                    for (UInteger elnum: adjacent)
                    {
                        Tetrahedron el = element_[elnum];
                        Point3D p[4];
                        for (int j = 0; j < 4; j++)
                        {
                            UInteger nnode = el[j];
                            if(nnode == i)
                            {
                                p[j].set(vars[0], vars[1], vars[2]);
                            }
                            else
                            {
                                p[j] = node_[nnode].point;
                            }
                        }
//                        double alpha = ((p[1] - p[0]).product(p[2] - p[0])) * (p[3] - p[0]);
//                        f += alpha*alpha;
                        double alpha = 0.0;
                        double t  = 10.0;
                        if (el.index(i) == 0)
                        {
                            alpha = ((p[1] - p[0]).product(p[2] - p[0])) * (p[3] - p[0]);
                            f += exp(-t * alpha);
                        }

                        if (el.index(i) == 1)
                        {
                            alpha = ((p[3] - p[1]).product(p[2] - p[1])) * (p[0] - p[1]);
                            f += exp(-t * alpha);
                        }

                        if (el.index(i) == 2)
                        {
                            alpha = ((p[1] - p[2]).product(p[3] - p[2])) * (p[0] - p[2]);
                            f += exp(-t * alpha);
                        }

                        if (el.index(i) == 3)
                        {
                            alpha = ((p[0] - p[3]).product(p[2] - p[3])) * (p[1] - p[3]);
                            f += exp(-t * alpha);
                        }
//                        The length functional
//                        f += (p[1] - p[0]) * (p[1] - p[0]) + (p[3] - p[0]) * (p[3] - p[0]) + (p[4] - p[0]) * (p[4] - p[0]);
//                        f += (p[0] - p[1]) * (p[0] - p[1]) + (p[5] - p[1]) * (p[5] - p[1]) + (p[2] - p[1]) * (p[2] - p[1]);
//                        f += (p[6] - p[2]) * (p[6] - p[2]) + (p[3] - p[2]) * (p[3] - p[2]) + (p[1] - p[2]) * (p[1] - p[2]);
//                        f += (p[2] - p[3]) * (p[2] - p[3]) + (p[7] - p[3]) * (p[7] - p[3]) + (p[0] - p[3]) * (p[0] - p[3]);
//                        f += (p[5] - p[4]) * (p[5] - p[4]) + (p[0] - p[4]) * (p[0] - p[4]) + (p[7] - p[4]) * (p[7] - p[4]);
//                        f += (p[6] - p[5]) * (p[6] - p[5]) + (p[1] - p[5]) * (p[1] - p[5]) + (p[4] - p[5]) * (p[4] - p[5]);
//                        f += (p[7] - p[6]) * (p[7] - p[6]) + (p[2] - p[6]) * (p[2] - p[6]) + (p[5] - p[6]) * (p[5] - p[6]);
//                        f += (p[6] - p[7]) * (p[6] - p[7]) + (p[4] - p[7]) * (p[4] - p[7]) + (p[3] - p[7]) * (p[3] - p[7]);
//                        std::cout<<area;
//                        std::cin>>i;
                    }
                    return f;
                };

                x0[0] = point.x();
                x0[1] = point.y();
                x0[2] = point.z();
                x = conjugateGradient(functor, x0, 0.001 * avr_dist, 0.000001 * avr_dist, 40, false);
//                node_[i].point.print();
                if (!isnan(x[0]) && !isnan(x[1]) && !isnan(x[2]))
                    node_[i].point.set(x[0], x[1], x[2]);
//                node_[i].point.println();
            }
            ++progress_bar;
        } // for i
    }
    double duration = static_cast<double>(std::clock() - start) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
}

double TetrahedralMesh3D::lengthAspect(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double min = node_[e->vertexNode(0)].point.distanceTo(node_[e->vertexNode(1)].point);
    double max = min;
    for (int i = 1; i < e->verticesCount(); i++)
    {
        double d = node_[e->vertexNode(i)].point.distanceTo(node_[e->vertexNode(i + 1)].point);
        if (d < min) min = d;
        if (d > max) max = d;
    }
    double d02 = node_[e->vertexNode(0)].point.distanceTo(node_[e->vertexNode(2)].point);
    double d13 = node_[e->vertexNode(1)].point.distanceTo(node_[e->vertexNode(3)].point);
    if (d02 < min) min = d02;
    if (d02 > max) max = d02;
    if (d13 < min) min = d13;
    if (d13 > max) max = d13;
    return min / max;
}

double TetrahedralMesh3D::minAngle(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double angle = 7.28;
    for (int i = 0; i < e->facesCount(); i++)
    {
        UIntegerVector face = e->face(i);
        Point3D A = node_[face[0]].point;
        Point3D B = node_[face[1]].point;
        Point3D C = node_[face[2]].point;
        double a = A.angle(C, B);
        double b = B.angle(A, C);
        double c = C.angle(B, A);
        if (a < angle) angle = a;
        if (b < angle) angle = b;
        if (c < angle) angle = c;
    }
    return angle;
}

double TetrahedralMesh3D::maxAngle(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double angle = 0.0;
    for (int i = 0; i < e->facesCount(); i++)
    {
        UIntegerVector face = e->face(i);
        Point3D A = node_[face[0]].point;
        Point3D B = node_[face[1]].point;
        Point3D C = node_[face[2]].point;
        double a = A.angle(C, B);
        double b = B.angle(A, C);
        double c = C.angle(B, A);
        if (a > angle) angle = a;
        if (b > angle) angle = b;
        if (c > angle) angle = c;
    }
    return angle;
}

void TetrahedralMesh3D::subdivide(std::list<UInteger> eNumbers, std::function<double (double, double, double)> func)
{
    int table [][33] = {
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 2
        { 0,  4,  2,  3,  4,  1,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 3
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 4
        { 0,  1,  6,  3,  1,  2,  6,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 5
        { 0,  1,  5,  3,  0,  5,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 6
        { 0,  4,  6,  3,  1,  5,  4,  3,  2,  6,  5,  3,  4,  5,  6,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 7
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 8
        { 0,  1,  2,  7,  7,  1,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 9
        { 0,  1,  2,  8,  0,  8,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 10
        { 0,  7,  4,  2,  1,  4,  8,  2,  3,  8,  7,  2,  4,  7,  8,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 11
        { 0,  1,  2,  9,  0,  1,  9,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 12
        { 0,  6,  7,  1,  2,  9,  6,  1,  3,  7,  9,  1,  6,  9,  7,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 13
        { 1,  8,  5,  0,  2,  5,  9,  0,  3,  9,  8,  0,  5,  8,  9,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 14
        { 0,  4,  6,  7,  1,  5,  4,  8,  2,  6,  5,  9,  7,  8,  9,  3,  4,  8,  6,  7,  6,  8,  9,  7,  4,  6,  8,  5,  6,  9,  8,  5, -1}  // 15
    };
    AdjacentSet refined_nodes;
    AdjacentSet refined_elements;
    Node3D local[10];
    for (UInteger elnum: eNumbers)
    {
        Tetrahedron element = element_[elnum];
        for (int i = 0; i < 4; i++)
        {
            AdjacentSet a = node_[element[i]].adjacent;
            refined_elements.insert(a.begin(), a.end());
            refined_nodes.insert(element[i]);
        }
    }
    for (UInteger elnum: refined_elements)
    {
        Tetrahedron el = element_[elnum];
        Node3D n0 = node_[el[0]];
        Node3D n1 = node_[el[1]];
        Node3D n2 = node_[el[2]];
        Node3D n3 = node_[el[3]];
        int code = 0;
        local[0] = n0;
        local[1] = n1;
        local[2] = n2;
        local[3] = n3;
        local[4].point = 0.5 * (n0.point + n1.point); local[4].type = (n0.type == INNER || n1.type == INNER) ? INNER : BORDER;
        local[5].point = 0.5 * (n1.point + n2.point); local[5].type = (n1.type == INNER || n2.type == INNER) ? INNER : BORDER;
        local[6].point = 0.5 * (n2.point + n0.point); local[6].type = (n2.type == INNER || n0.type == INNER) ? INNER : BORDER;
        local[7].point = 0.5 * (n0.point + n3.point); local[7].type = (n0.type == INNER || n3.type == INNER) ? INNER : BORDER;
        local[8].point = 0.5 * (n1.point + n3.point); local[8].type = (n1.type == INNER || n3.type == INNER) ? INNER : BORDER;
        local[9].point = 0.5 * (n2.point + n3.point); local[9].type = (n2.type == INNER || n3.type == INNER) ? INNER : BORDER;
        if (func != nullptr)
        {
            for (int i = 4; i < 10; i++)
                if (local[i].type == BORDER)
                    local[i].point = findBorder(local[i].point, func, 0.5 * (n0.point.distanceTo(n1.point) + n0.point.distanceTo(n2.point)));
        }
        if (refined_nodes.find(el[0]) != refined_nodes.end())
            code |= 1;
        if (refined_nodes.find(el[1]) != refined_nodes.end())
            code |= 2;
        if (refined_nodes.find(el[2]) != refined_nodes.end())
            code |= 4;
        if (refined_nodes.find(el[3]) != refined_nodes.end())
            code |= 8;
        if (code != 0 && code != 1 && code != 2 && code != 4 && code != 8)
        {
            node_[el[0]].adjacent.erase(elnum);
            node_[el[1]].adjacent.erase(elnum);
            node_[el[2]].adjacent.erase(elnum);
            node_[el[3]].adjacent.erase(elnum);
            for(int i = 0; table[code][i] != -1; i += 4)
            {
                el[0] = addNode(local[table[code][i]]);
                el[1] = addNode(local[table[code][i + 1]]);
                el[2] = addNode(local[table[code][i + 2]]);
                el[3] = addNode(local[table[code][i + 3]]);
                if (i == 0)
                {
                    element_[elnum] = el;
                    node_[el[0]].adjacent.insert(elnum);
                    node_[el[1]].adjacent.insert(elnum);
                    node_[el[2]].adjacent.insert(elnum);
                    node_[el[3]].adjacent.insert(elnum);
                }
                else
                {
                    addElement(el);
                }
            }
        }
    }

//    for (UInteger i = 0; i < nodesCount(); i++) node_[i].type = INNER;
//    for (UInteger i = 0; i < elementsCount(); i++)
//    {
//        Tetrahedron el = element_[i];
//        Node3D n0 = node_[el[0]];
//        Node3D n1 = node_[el[1]];
//        Node3D n2 = node_[el[2]];
//        Node3D n3 = node_[el[3]];
//        for (int j = 0; j < el.facesCount(); j++)
//        {
//            UIntegerVector face = el.face(j);
//            int c = 0;
//            for (UInteger k = 0; k < elementsCount(); k++)
//            {
//                Tetrahedron el2 = element_[k];
//                if (k != i && el2.in(face[0]) && el2.in(face[1]) && el2.in(face[2]))
//                    c++;
//            }
//            if (c == 0)
//            {
//                for (int k = 0; k < face.size(); k++)
//                {
//                    node_[face[k]].type = BORDER;
//                    if (func != nullptr)
//                        node_[face[k]].point = findBorder(node_[face[k]].point, func, 0.5 * (n0.point.distanceTo(n1.point) + n0.point.distanceTo(n2.point)));
//                }
//            }
//        }
//    }
}

}

