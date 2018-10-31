#include "hexahedralmesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <iostream>
#include <float.h>

#include <map>
#include <ctime>

//#include <omp.h>

#include "quadrilateralmesh3d.h"
#include "funcopt.h"
#include "consoleprogress.h"

namespace msh {
HexahedralMesh3D::HexahedralMesh3D() : Mesh3D(NULL)
{
}

void HexahedralMesh3D::prismDomain(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth)
{
    clear();
    double hx = width / static_cast<double>(xCount - 1);
    double hy = height / static_cast<double>(yCount - 1);
    double hz = depth / static_cast<double>(zCount - 1);
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        double x = xMin + static_cast<double>(i) * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double y = yMin + static_cast<double>(j) * hy;
            for (UInteger k = 0; k < zCount; k++)
            {
                double z = zMin + static_cast<double>(k) * hz;
                Point3D point(x, y, z);

                //                if ((i == 0 && j == 0) || (i == 0 & j == yCount - 1) || (i == xCount - 1 && j == 0) || (i == xCount - 1 & j == yCount - 1))
                //                    pushNode(point, CHARACTER);
                //                else
                if (i == 0 || j == 0 || i == xCount - 1 || j == yCount - 1 || k == 0 || k == zCount - 1)
                    pushNode(point, BORDER);
                else
                    pushNode(point, INNER);
            }
        }
    }
    // функция для определения индекса в одномерном массиве, который соответствует индексам трехмерной матрицы (C++0x)
    auto toArray = [](UInteger i, UInteger j, UInteger k, UInteger yCount, UInteger zCount) { return i * yCount * zCount + j * zCount + k; };
    // формирование массива элементов
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            for (UInteger k = 0; k < zCount - 1; k++)
            {
                addElement(toArray(i, j, k, yCount, zCount), toArray(i, j, k+1, yCount, zCount), toArray(i+1, j, k+1, yCount, zCount), toArray(i+1, j, k, yCount, zCount),
                           toArray(i, j+1, k, yCount, zCount), toArray(i, j+1, k+1, yCount, zCount), toArray(i+1, j+1, k+1, yCount, zCount), toArray(i+1, j+1, k, yCount, zCount));
            }
        }
    }
//    std::list<UInteger> ee;
//    ee.push_back(0); ee.push_back(1); ee.push_back(20);
//    subdivide(ee, nullptr);
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    zMin_ = zMin;
    zMax_ = xMin + depth;
    printStats();
}

void HexahedralMesh3D::rotateBaseMesh(QuadrilateralMesh2D *baseMesh, const double &xDelta, const double &yDelta, const int &lCount, bool x_axes, bool withLayersInfo)
{
    clear();
    xMin_ = DBL_MAX;
    xMax_ = -DBL_MAX;
    yMin_ = DBL_MAX;
    yMax_ = -DBL_MAX;
    zMin_ = DBL_MAX;
    zMax_ = -DBL_MAX;

    node_.reserve(baseMesh->nodesCount() *(lCount + 1));
    element_.reserve(baseMesh->elementsCount() * lCount);

    UInteger i;
    // текущий угол поворта
    double phi;
    double delta_phi = 2.0 * M_PI / static_cast<double>(lCount);
    // номера вершин шестигранника
    UInteger nodes_pointers[8];
    UInteger baseNodesCount = baseMesh->nodesCount();

    // формирование узлов
    for(int iz = 0; iz < lCount; iz++)
    {
        phi = static_cast<double>(iz) * delta_phi;

        for(i = 0; i < baseNodesCount; i++)
        {
            // текущая вершина на плоскости
            Point2D planeNode(baseMesh->node(i)->x(), baseMesh->node(i)->y());
            // текущая вершина в пространстве
            Point3D spaceNode;
            NodeType type;

            type = baseMesh->nodeType(i);

            if(x_axes)
            {
                spaceNode.set(xDelta + planeNode.x(), (yDelta + planeNode.y()) * cos(phi), (yDelta + planeNode.y()) * sin(phi));
            }
            else
            {
                spaceNode.set((xDelta + planeNode.x()) * cos(phi), yDelta + planeNode.y(), (xDelta + planeNode.x()) * sin(phi));
            }

            if(xMin_ > spaceNode.x())
                xMin_ = spaceNode.x();

            if(xMax_ < spaceNode.x())
                xMax_ = spaceNode.x();

            if(yMin_ > spaceNode.y())
                yMin_ = spaceNode.y();

            if(yMax_ < spaceNode.y())
                yMax_ = spaceNode.y();

            if(zMin_ > spaceNode.z())
                zMin_ = spaceNode.z();

            if(zMax_ < spaceNode.z())
                zMax_ = spaceNode.z();

            pushNode(spaceNode, type);
        }
    }

    // формирование шестигранников
    for (int iz = 0; iz < lCount; iz++)
    {
        for (i = 0; i < baseMesh->elementsCount(); i++)
        {
            Quadrilateral current_quad = baseMesh->quadrilateral(i);
            nodes_pointers[0] = current_quad[0] + iz * baseNodesCount;
            nodes_pointers[1] = current_quad[1] + iz * baseNodesCount;
            nodes_pointers[2] = current_quad[2] + iz * baseNodesCount;
            nodes_pointers[3] = current_quad[3] + iz * baseNodesCount;
            if(iz < lCount - 1)
            {
                nodes_pointers[4] = current_quad[0] + (iz + 1) * baseNodesCount;
                nodes_pointers[5] = current_quad[1] + (iz + 1) * baseNodesCount;
                nodes_pointers[6] = current_quad[2] + (iz + 1) * baseNodesCount;
                nodes_pointers[7] = current_quad[3] + (iz + 1) * baseNodesCount;
            }
            else
            {
                nodes_pointers[4] = current_quad[0];
                nodes_pointers[5] = current_quad[1];
                nodes_pointers[6] = current_quad[2];
                nodes_pointers[7] = current_quad[3];
            }

//            if (x_axes)
                addElement(nodes_pointers[0], nodes_pointers[1], nodes_pointers[2], nodes_pointers[3], nodes_pointers[4], nodes_pointers[5], nodes_pointers[6], nodes_pointers[7]);
//            else
//                addElement(nodes_pointers[1], nodes_pointers[2], nodes_pointers[6], nodes_pointers[5], nodes_pointers[0], nodes_pointers[3], nodes_pointers[7], nodes_pointers[4]);
            if (withLayersInfo && baseMesh->sizeOfLayers() == baseMesh->elementsCount())
                pushLayer(baseMesh->layer(i));
        }
    }
    printStats();
}

void HexahedralMesh3D::rotateBaseMesh(QuadrilateralMesh2D *baseMesh, const double &xDelta, const double &yDelta, const double &angle, const int &lCount, bool x_axes, bool withLayersInfo)
{
    clear();
    xMin_ = DBL_MAX;
    xMax_ = -DBL_MAX;
    yMin_ = DBL_MAX;
    yMax_ = -DBL_MAX;
    zMin_ = DBL_MAX;
    zMax_ = -DBL_MAX;

    node_.reserve(baseMesh->nodesCount() * (lCount + 1L));
    element_.reserve(baseMesh->elementsCount() * lCount);

    UInteger i;
    // текущий угол поворта
    double phi;
    double delta_phi = angle * M_PI / (static_cast<double>(lCount) * 180.0);
    // номера вершин шестигранника
    UInteger nodes_pointers[8];
    UInteger baseNodesCount = baseMesh->nodesCount();

    // формирование узлов
    for(int iz = 0; iz <= lCount; iz++)
    {
        phi = static_cast<double>(iz) * delta_phi;

        for(i = 0; i < baseNodesCount; i++)
        {
            // текущая вершина на плоскости
            Point2D planeNode(baseMesh->node(i)->x(), baseMesh->node(i)->y());
            // текущая вершина в пространстве
            Point3D spaceNode;
            NodeType type;

            type = baseMesh->nodeType(i);

            if(x_axes)
            {
                spaceNode.set(xDelta + planeNode.x(), (yDelta + planeNode.y()) * cos(phi), (yDelta + planeNode.y()) * sin(phi));
            }
            else
            {
                spaceNode.set((xDelta + planeNode.x()) * cos(phi), yDelta + planeNode.y(), (xDelta + planeNode.x()) * sin(phi));
            }

            if(xMin_ > spaceNode.x())
                xMin_ = spaceNode.x();

            if(xMax_ < spaceNode.x())
                xMax_ = spaceNode.x();

            if(yMin_ > spaceNode.y())
                yMin_ = spaceNode.y();

            if(yMax_ < spaceNode.y())
                yMax_ = spaceNode.y();

            if(zMin_ > spaceNode.z())
                zMin_ = spaceNode.z();

            if(zMax_ < spaceNode.z())
                zMax_ = spaceNode.z();

            pushNode(spaceNode, (iz == 0 || iz == lCount) ? BORDER : type);
        }
    }

    // формирование шестигранников
    for(int iz = 0; iz < lCount; iz++)
    {
        for(i = 0; i < baseMesh->elementsCount(); i++)
        {
            Quadrilateral current_quad = baseMesh->quadrilateral(i);
            if((x_axes && angle < 0.0) || (!x_axes && angle > 0.0))
            {
                nodes_pointers[1] = current_quad[0] + iz * baseNodesCount;
                nodes_pointers[2] = current_quad[1] + iz * baseNodesCount;
                nodes_pointers[6] = current_quad[2] + iz * baseNodesCount;
                nodes_pointers[5] = current_quad[3] + iz * baseNodesCount;
                nodes_pointers[0] = current_quad[0] + (iz + 1) * baseNodesCount;
                nodes_pointers[3] = current_quad[1] + (iz + 1) * baseNodesCount;
                nodes_pointers[7] = current_quad[2] + (iz + 1) * baseNodesCount;
                nodes_pointers[4] = current_quad[3] + (iz + 1) * baseNodesCount;
            }
            else
            {
                nodes_pointers[0] = current_quad[0] + iz * baseNodesCount;
                nodes_pointers[1] = current_quad[1] + iz * baseNodesCount;
                nodes_pointers[2] = current_quad[2] + iz * baseNodesCount;
                nodes_pointers[3] = current_quad[3] + iz * baseNodesCount;
                nodes_pointers[4] = current_quad[0] + (iz + 1) * baseNodesCount;
                nodes_pointers[5] = current_quad[1] + (iz + 1) * baseNodesCount;
                nodes_pointers[6] = current_quad[2] + (iz + 1) * baseNodesCount;
                nodes_pointers[7] = current_quad[3] + (iz + 1) * baseNodesCount;
            }

            addElement(nodes_pointers[0], nodes_pointers[1], nodes_pointers[2], nodes_pointers[3], nodes_pointers[4], nodes_pointers[5], nodes_pointers[6], nodes_pointers[7]);
            if (withLayersInfo && baseMesh->sizeOfLayers() == baseMesh->elementsCount())
                pushLayer(baseMesh->layer(i));
        }
    }
    printStats();
}

HexahedralMesh3D::HexahedralMesh3D(const HexahedralMesh3D &mesh) : Mesh3D(&mesh)
{
    element_ = mesh.element_;
}

HexahedralMesh3D::HexahedralMesh3D(const HexahedralMesh3D *mesh): Mesh3D(mesh)
{
    element_ = mesh->element_;
}

UInteger HexahedralMesh3D::elementsCount() const
{
    return element_.size();
}

ElementPointer HexahedralMesh3D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

double HexahedralMesh3D::surfaceArea() const
{
    double area = 0.0;
    for (UInteger i = 0; i < elementsCount(); i++)
    {
        Hexahedral hex = element_[i];
        for (int j = 0; j < hex.facesCount(); j++)
        {
            UIntegerVector face = hex.face(j);
            if (isBorderFace(face)) area += faceArea(face);
        }
    }
    return area;
}

void HexahedralMesh3D::addElement(const Hexahedral &hex)
{
    element_.push_back(hex);
    // обновление списка смежных узлов
    node_[hex[0]].adjacent.insert(element_.size() - 1);
    node_[hex[1]].adjacent.insert(element_.size() - 1);
    node_[hex[2]].adjacent.insert(element_.size() - 1);
    node_[hex[3]].adjacent.insert(element_.size() - 1);
    node_[hex[4]].adjacent.insert(element_.size() - 1);
    node_[hex[5]].adjacent.insert(element_.size() - 1);
    node_[hex[6]].adjacent.insert(element_.size() - 1);
    node_[hex[7]].adjacent.insert(element_.size() - 1);
}

void HexahedralMesh3D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3, const UInteger &node4, const UInteger &node5, const UInteger &node6, const UInteger &node7)
{
    Hexahedral hex(node0, node1, node2, node3, node4, node5, node6, node7);
    element_.push_back(hex);
    // обновление списка смежных узлов
    node_[node0].adjacent.insert(element_.size() - 1);
    node_[node1].adjacent.insert(element_.size() - 1);
    node_[node2].adjacent.insert(element_.size() - 1);
    node_[node3].adjacent.insert(element_.size() - 1);
    node_[node4].adjacent.insert(element_.size() - 1);
    node_[node5].adjacent.insert(element_.size() - 1);
    node_[node6].adjacent.insert(element_.size() - 1);
    node_[node7].adjacent.insert(element_.size() - 1);
}

void HexahedralMesh3D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0],
            nodes_ref[1],
            nodes_ref[2],
            nodes_ref[3],
            nodes_ref[4],
            nodes_ref[5],
            nodes_ref[6],
            nodes_ref[7]);
}

void HexahedralMesh3D::clearElements()
{
    element_.clear();
}

void HexahedralMesh3D::sweepBaseMesh(QuadrilateralMesh2D *baseMesh, const double &z0, const double &z1, const double &phi0, const double &phi1, const double &k0, const double &k1, const int &zLayersCount)
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
    // формирование шестигранников
    for(int i = 0; i < zLayersCount; i++)
    {
        for(UInteger j = 0; j < baseMesh->elementsCount(); j++)
        {
            Quadrilateral quad = baseMesh->quadrilateral(j);
            std::vector<UInteger> nnode(8);
            nnode[0] = quad[0] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            nnode[1] = quad[1] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            nnode[2] = quad[2] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            nnode[3] = quad[3] + static_cast<UInteger>(i) * baseMesh->nodesCount();
            nnode[4] = quad[0] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            nnode[5] = quad[1] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            nnode[6] = quad[2] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            nnode[7] = quad[3] + static_cast<UInteger>(i + 1) * baseMesh->nodesCount();
            addElement(nnode);
        }
    }
    printStats();
    updateDomain();
}

Hexahedral HexahedralMesh3D::hexahedron(const UInteger &ie) const
{
    return element_[ie];
}

void HexahedralMesh3D::printStats() const
{
    std::cout << "Hexahedral mesh: " << nodesCount() << " node(s), " << elementsCount() << " element(s)." << std::endl;
}

void HexahedralMesh3D::delNode(const UInteger &nnumber)
{
    if (nnumber >= nodesCount()) return;
    if (!node_[nnumber].adjacent.empty())
    {
        AdjacentSet A = node_[nnumber].adjacent;
        UInteger shift = 0;
        //        std::cout << A.size();
        for (UInteger elnum: A)
        {
            delElement(elnum - shift);
            shift++;
        }
        //        while (!node_[nnumber].adjacent.empty())
        //        {
        //            UInteger elnum = *(node_[nnumber].adjacent.begin());
        //            std::cout /*<< node_[nnumber].adjacent.size()*/ << '!';
        //            delElement(elnum);
        //        }
    }
    {
        //        for (UInteger i = nnumber + 1; i < nodesCount(); i++)
        //        {
        //            AdjacentSet A = node_[i].adjacent;
        //            for (UInteger elnum: A)
        //            {
        //                Hexahedral h = element_[elnum];
        //                int index = h.index(i);
        //                if (index != -1)
        //                {
        //                    element_[elnum][index] = i - 1;
        //                }
        //            }
        //        }
        for (std::vector<Hexahedral>::iterator el = element_.begin(); el != element_.end(); el++)
        {
            for (int i = 0; i < 8; i++)
            {
                if ((*el)[i] > nnumber)
                    (*el)[i] -= 1;
            }
        }
        //#pragma omp parallel for
        //        for (UInteger j = 0; j < elementsCount(); j++)
        //        {
        //            for (int i = 0; i < 8; i++)
        //            {
        //                if (element_[j][i] > nnumber)
        //                    element_[j][i] -= 1;
        //            }
        //        }
        node_.erase(node_.begin() + nnumber);
    }
}

void HexahedralMesh3D::delElement(const UInteger &elnum)
{
    /*Hexahedral h = element_[elnum];
    node_[h[0]].adjacent.remove(elnum);
    node_[h[1]].adjacent.remove(elnum);
    node_[h[2]].adjacent.remove(elnum);
    node_[h[3]].adjacent.remove(elnum);
    node_[h[4]].adjacent.remove(elnum);
    node_[h[5]].adjacent.remove(elnum);
    node_[h[6]].adjacent.remove(elnum);
    node_[h[7]].adjacent.remove(elnum);
    for (int i = 0; i < 8; i++)
    {
        node_[h[i]].adjacent.remove(elnum);
//        std::cout << h[i] << ' ';
//        if (node_[h[i]].adjacent.empty())
//        {
////            std::cout << ">>";
//            delNode(h[i]);
//        }
    }
    for (std::vector<Node3D>::iterator it = node_.begin(); it != node_.end(); it++)
    {
//        AdjacentSet A;
        for (AdjacentSet::iterator e = (*it).adjacent.begin(); e != (*it).adjacent.end(); e++)
        {
            if ((*e) > elnum)
            {
//                A.insert((*e) - 1);
                (*e) -= 1;
            }
//            else
//            {
//                A.insert((*e));
//            }
        }
//        (*it).adjacent = A;
    }
//#pragma omp parallel for
//    for (UInteger i = 0; i < nodesCount(); i++)
//    {
//        AdjacentSet A;
//        for (AdjacentSet::iterator e = node_[i].adjacent.begin(); e != node_[i].adjacent.end(); e++)
//        {
//            if ((*e) > elnum)
//            {
//                A.insert((*e) - 1);
//            }
//            else
//            {
//                A.insert((*e));
//            }
//        }
//        node_[i].adjacent = A;
//    }
    element_.erase(element_.begin() + elnum);*/
}

void HexahedralMesh3D::clearFuncNodes(std::function<double (double, double, double)> func, UInteger maxCount)
{
    UInteger i = 0;
    UInteger c = nodesCount();
    UInteger nc = c;
    while (i < nc && c - nc < maxCount)
    {
        Point3D p = node_[i].point;
        if (func(p.x(), p.y(), p.z()) < 0.0)
        {
            delNode(i);
            //            p.print();
            //            std::cout << ' ' << func(p.x(), p.y(), p.z()) << std::endl;
        }
        else
        {
            i += 1;
        }
        nc = nodesCount();
    }
}

void HexahedralMesh3D::backgroundGrid(const HexahedralMesh3D *mesh, std::function<double (double, double, double)> func, double level, int smooth, int optimize)
{
    clear();
    QuadrilateralMesh3D qmesh;
    std::list<ElementPointer> inner = qmesh.backgroundGrid(mesh, func, level, smooth, optimize, false);
    std::map<UInteger, UInteger> iso;
    for (ElementPointer el: inner)
    {
        for (int j = 0; j < el->facesCount(); j++)
        {
            UIntegerVector f = el->face(j);
            Point3D p[4];
            p[0] = mesh->point3d(f[0]);
            p[1] = mesh->point3d(f[1]);
            p[2] = mesh->point3d(f[2]);
            p[3] = mesh->point3d(f[3]);
            bool isInner = false;
            // check an other element with these vertices
            for (ElementPointer el1: inner)
            {
                if (el1 != el && el1->in(f[0]) && el1->in(f[1]) && el1->in(f[2]) && el1->in(f[3]))
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
                iso[f[3]] = addNode(p[3], INNER);
            }
        }
    }
    if (nodesCount() == qmesh.nodesCount())
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
        std::vector<UInteger> v(8);
        for (int i = 0; i < 8; i++)
        {
            v[i] = addNode(mesh->point3d(el->vertexNode(i)), INNER);
        }
        addElement(v);
    }
    UInteger shift = nodesCount();
    for (UInteger i = 0; i < qmesh.nodesCount(); i++)
    {
        pushNode(qmesh.node3d(i));
    }
    for (ElementPointer el: inner)
    {
        for (int j = 0; j < el->facesCount(); j++)
        {
            UIntegerVector f = el->face(j);
            Point3D p[4];
            p[0] = mesh->point3d(f[0]);
            p[1] = mesh->point3d(f[1]);
            p[2] = mesh->point3d(f[2]);
            p[3] = mesh->point3d(f[3]);
            bool isInner = false;
            // check an other element with these vertices
            for (ElementPointer el1: inner)
            {
                if (el1 != el && el1->in(f[0]) && el1->in(f[1]) && el1->in(f[2]) && el1->in(f[3]))
                {
                    isInner = true;
                    break;
                }
            }
            if (!isInner)
            {
                addElement(iso[f[0]], iso[f[1]], iso[f[2]], iso[f[3]],
                        iso[f[0]] + shift, iso[f[1]] + shift, iso[f[2]] + shift, iso[f[3]] + shift);
            }
        }
    }
    laplacianSmoothing(smooth);
    localFubctionalOptimization(optimize);
//    std::list<UInteger> ee;
//    ee.push_back(4415); ee.push_back(4414); ee.push_back(4413);
//    subdivide(ee, func);
    printStats();
    updateDomain();
}

void HexahedralMesh3D::localFubctionalOptimization(int maxiter, double t)
{
    std::cout << "Local optimization smoothing...";
    std::clock_t start = std::clock();
    // проект сглаживания путем локальной минимизации функционала
    for (int iii = 0; iii < maxiter; iii++)
    {
        //        double sum = 0.0;
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
                    Hexahedral el = element_[elnum];
                    int index = el.index(i);
                    Point3D p0 = node_[el[index + 3]].point;
                    Point3D p1 = node_[el[index + 1]].point;
                    Point3D p2 = node_[el[index + 4]].point;
                    avr_dist += point.distanceTo(p0) + point.distanceTo(p1) + point.distanceTo(p2);
                }

                avr_dist /= static_cast<double>(3 * adjacent.size());

                auto functor = [&](const std::vector<double> &vars){
                    double f = 0.0;
                    for (UInteger elnum: adjacent)
                    {
                        Hexahedral el = element_[elnum];
                        Point3D p[8];
                        for (int j = 0; j < 8; j++)
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
                        // Рассмотрим 4х-угольник как 4 треугольника, определенных на его углах
                        double area = ((p[1] - p[0]).product(p[3] - p[0])) * (p[4] - p[0]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 0)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[5] - p[1]).product(p[2] - p[1])) * (p[0] - p[1]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 1)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[6] - p[2]).product(p[3] - p[2])) * (p[1] - p[2]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 2)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[2] - p[3]).product(p[7] - p[3])) * (p[0] - p[3]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 3)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[5] - p[4]).product(p[0] - p[4])) * (p[7] - p[4]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 4)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[6] - p[5]).product(p[1] - p[5])) * (p[4] - p[5]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 5)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[5] - p[6]).product(p[7] - p[6])) * (p[2] - p[6]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 6)
                            f += exp(-t * area);
                        //                        f += area*area;
                        area = ((p[6] - p[7]).product(p[4] - p[7])) * (p[3] - p[7]);
                        //                        std::cout << area << std::endl;
                        if (el.index(i) == 7)
                            f += exp(-t * area);
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
                //                sum += functor(x0);
                //                x = descentGradient(functor, x0, 0.001 * avr_dist, 0.00001 * avr_dist, 50, false);
                x = conjugateGradient(functor, x0, 0.001 * avr_dist, 0.00001 * avr_dist, 40, false);
                //                node_[i].point.print();
                if (!isnan(x[0]) && !isnan(x[1]) && !isnan(x[2]))
                    node_[i].point.set(x[0], x[1], x[2]);
                //                node_[i].point.println();
            }
            ++progress_bar;
        } // for i
        //        std::cout << "F = " << sum;
    }
    double duration = static_cast<double>(std::clock() - start) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
}

double HexahedralMesh3D::lengthAspect(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double min = node_[e->vertexNode(0)].point.distanceTo(node_[e->vertexNode(1)].point);
    double max = min;
    for (int i = 0; i < 3; i++)
    {
        double d0 = node_[e->vertexNode(i)].point.distanceTo(node_[e->vertexNode(i + 1)].point);
        double d1 = node_[e->vertexNode(i + 4)].point.distanceTo(node_[e->vertexNode(i + 4 + 1)].point);
        if (d0 < min) min = d0;
        if (d0 > max) max = d0;
        if (d1 < min) min = d1;
        if (d1 > max) max = d1;
    }
    double d30 = node_[e->vertexNode(3)].point.distanceTo(node_[e->vertexNode(0)].point);
    double d74 = node_[e->vertexNode(7)].point.distanceTo(node_[e->vertexNode(4)].point);
    if (d30 < min) min = d30;
    if (d30 > max) max = d30;
    if (d74 < min) min = d74;
    if (d74 > max) max = d74;
    for (int i = 0; i < 4; i++)
    {
        double d = node_[e->vertexNode(i)].point.distanceTo(node_[e->vertexNode(i + 4)].point);
        if (d < min) min = d;
        if (d > max) max = d;
    }
    return min / max;
}

double HexahedralMesh3D::minAngle(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double angle = 7.28;
    for (int i = 0; i < e->facesCount(); i++)
    {
        UIntegerVector face = e->face(i);
        Point3D A = node_[face[0]].point;
        Point3D B = node_[face[1]].point;
        Point3D C = node_[face[2]].point;
        Point3D D = node_[face[3]].point;
        double a = A.angle(D, B);
        double b = B.angle(A, C);
        double c = C.angle(B, D);
        double d = D.angle(C, A);
        if (a < angle) angle = a;
        if (b < angle) angle = b;
        if (c < angle) angle = c;
        if (d < angle) angle = d;
    }
    return angle;
}

double HexahedralMesh3D::maxAngle(const UInteger &elnum) const
{
    ElementPointer e = element(elnum);
    double angle = 0.0;
    for (int i = 0; i < e->facesCount(); i++)
    {
        UIntegerVector face = e->face(i);
        Point3D A = node_[face[0]].point;
        Point3D B = node_[face[1]].point;
        Point3D C = node_[face[2]].point;
        Point3D D = node_[face[3]].point;
        double a = A.angle(D, B);
        double b = B.angle(A, C);
        double c = C.angle(B, D);
        double d = D.angle(C, A);
        if (a > angle) angle = a;
        if (b > angle) angle = b;
        if (c > angle) angle = c;
        if (d > angle) angle = d;
    }
    return angle;
}

void HexahedralMesh3D::subdivide(std::list<UInteger> eNumbers, std::function<double (double, double, double)> func)
{
    const int HexTable [][3][105] = {
        /*00000000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00000001*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00000010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00000011*/{{0,0,2,3, 0,0,2,3, 0,0,2,2, 0,0,2,2, 0,0,3,2, 0,0,3,2, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,3,3, -1},
                     {0,0,0,0, 3,2,2,3, 0,0,0,0, 2,2,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, -1},
                     {0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 1,2,2,1, 0,3,3,0, -1}},
        /*00000100*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00000101*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00000110*/{{0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,2,2,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,3,1, 0,3,3,1, 1,3,3,1, 1,3,3,1, 1,3,3,0, 1,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 0,3,3,0, -1}},
        /*00000111*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00001000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00001001*/{{0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,0,3,2, 1,0,3,2, 1,1,2,2, 0,0,3,3, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,2,2,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 2,3,3,2, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,2,0, 0,3,2,0, 0,2,2,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,3,3,0, -1}},
        /*00001010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00001011*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00001100*/{{0,0,1,1, 0,0,1,1, 0,1,3,3, 0,1,3,3, 1,1,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 3,2,2,3, 0,0,0,0, 2,2,2,2, 0,0,0,0, 2,3,3,2, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, -1}},
        /*00001101*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00001110*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00001111*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 3, 3, 1, 1, 2, 2, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 2, 1, 2, 0, 0, 0, 0, 2, 2, 1, 1, 0, 0, 0, 0, 2, 3, 2, 1, 0, 0, 0, 0, 2, 1, 1, 2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 2, 2, 1, 0, 0, 0, 0, 2, 1, 2, 3, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 1, 2, 3, 2, 2, 2, 1, 1, 3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, -1},
                     {0, 1, 1, 0, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 0, 1, 1, 0, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 0, 1, 1, 0, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 1, 2, 2, 1, 0, 3, 3, 0, 1, 2, 2, 1, 0, 3, 3, 0, 1, 2, 2, 1, 0, 3, 3, 0, 0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00010000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00010001*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 2,3,3,2, 2,3,3,2, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,2,2,0, 0,2,2,0, -1}},
        /*00010010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00010011*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00010100*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00010101*/{{-1},
                     {-1},
                     {-1}},
        /*00010110*/{{-1},
                     {-1},
                     {-1}},
        /*00010111*/{{-1},
                     {-1},
                     {-1}},
        /*00011000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00011001*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*00011010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00011011*/{{-1},
                     {-1},
                     {-1}},
        /*00011100*/{{-1},
                     {-1},
                     {-1}},
        /*00011101*/{{-1},
                     {-1},
                     {-1}},
        /*00011110*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00011111*/{{-1},
                     {-1},
                     {-1}},
        /*00100000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00100001*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00100010*/{{0,0,3,3, 0,0,2,2, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, 0,0,2,2, 0,0,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,3,0, 1,3,3,1, 0,1,1,0, 0,1,1,0, 1,3,3,0, 1,3,3,0, 1,3,3,1, 1,3,3,1, 1,3,3,1, 0,3,3,0, -1}},
        /*00100011*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00100100*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00100101*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00100110*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*00100111*/{{-1},
                     {-1},
                     {-1}},
        /*00101000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*00101001*/{{-1},
                     {-1},
                     {-1}},
        /*00101010*/{{-1},
                     {-1},
                     {-1}},
        /*00101011*/{{-1},
                     {-1},
                     {-1}},
        /*00101100*/{{-1},
                     {-1},
                     {-1}},
        /*00101101*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00101110*/{{-1},
                     {-1},
                     {-1}},
        /*00101111*/{{-1},
                     {-1},
                     {-1}},
        /*00110000*/{{0,0,3,3, 0,0,2,2, 0,0,2,3, 0,0,2,3, 0,0,2,2, 0,0,2,2, 0,0,3,2, 0,0,3,2, 2,2,3,3, 2,2,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,3,3,3, 1,1,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, -1},
                     {0,3,3,0, 1,2,2,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, -1}},
        /*00110001*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00110010*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00110011*/{{0,0,2,3, 0,0,1,2, 0,0,2,2, 0,0,1,1, 0,0,3,2, 0,0,2,1, 2,2,3,3, 1,1,2,2, 0,0,1,2, 0,0,1,2, 0,0,1,1, 0,0,1,1, 0,0,2,1, 0,0,2,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 0,0,1,2, 0,0,2,3, 0,0,1,1, 0,0,2,2, 0,0,2,1, 0,0,3,2, 1,1,2,2, 2,2,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, -1},
                     {0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,3,3,0, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, -1}},
        /*00110100*/{{-1},
                     {-1},
                     {-1}},
        /*00110101*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00110110*/{{-1},
                     {-1},
                     {-1}},
        /*00110111*/{{-1},
                     {-1},
                     {-1}},
        /*00111000*/{{-1},
                     {-1},
                     {-1}},
        /*00111001*/{{-1},
                     {-1},
                     {-1}},
        /*00111010*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*00111011*/{{-1},
                     {-1},
                     {-1}},
        /*00111100*/{{-1},
                     {-1},
                     {-1}},
        /*00111101*/{{-1},
                     {-1},
                     {-1}},
        /*00111110*/{{-1},
                     {-1},
                     {-1}},
        /*00111111*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, -1}},
        /*01000000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01000001*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01000010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01000011*/{{-1},
                     {-1},
                     {-1}},
        /*01000100*/{{0,0,1,1, 0,0,1,1, 0,0,3,3, 1,1,3,3, 0,1,3,3, 0,1,3,3, 1,1,3,3, 1,1,3,3, 1,1,3,3, 0,0,3,3, -1},
                     {0,0,1,1, 3,3,2,2, 0,0,0,0, 1,1,1,1, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,3,1, 0,3,3,1, 0,3,3,0, 1,3,3,1, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, 1,3,3,1, 0,3,3,0, -1}},
        /*01000101*/{{-1},
                     {-1},
                     {-1}},
        /*01000110*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*01000111*/{{-1},
                     {-1},
                     {-1}},
        /*01001000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01001001*/{{-1},
                     {-1},
                     {-1}},
        /*01001010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01001011*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*01001100*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*01001101*/{{-1},
                     {-1},
                     {-1}},
        /*01001110*/{{-1},
                     {-1},
                     {-1}},
        /*01001111*/{{-1},
                     {-1},
                     {-1}},
        /*01010000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01010001*/{{-1},
                     {-1},
                     {-1}},
        /*01010010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01010011*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*01010100*/{{-1},
                     {-1},
                     {-1}},
        /*01010101*/{{-1},
                     {-1},
                     {-1}},
        /*01010110*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*01010111*/{{-1},
                     {-1},
                     {-1}},
        /*01011000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01011001*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*01011010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*01011011*/{{-1},
                     {-1},
                     {-1}},
        /*01011100*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*01011101*/{{-1},
                     {-1},
                     {-1}},
        /*01011110*/{{-1},
                     {-1},
                     {-1}},
        /*01011111*/{{-1},
                     {-1},
                     {-1}},
        /*01100000*/{{0,0,3,3, 1,1,2,2, 0,1,2,3, 0,1,2,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,3,3,3, 0,0,1,1, 3,3,3,3, 1,1,1,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, -1},
                     {0,3,3,0, 1,3,3,1, 0,1,1,0, 0,1,1,0, 0,3,3,1, 0,3,3,1, 1,3,3,1, 1,3,3,1, 1,3,3,0, 1,3,3,0, -1}},
        /*01100001*/{{-1},
                     {-1},
                     {-1}},
        /*01100010*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*01100011*/{{-1},
                     {-1},
                     {-1}},
        /*01100100*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*01100101*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*01100110*/{{0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,1,1,0, 3,2,2,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,3,1, 1,3,3,2, 1,3,3,1, 2,3,3,2, 1,3,3,0, 2,3,3,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, 0,1,2,1, 0,1,2,1, 1,2,1,0, 1,2,1,0, 1,3,3,2, 1,3,3,2, 2,3,3,2, 2,3,3,2, 2,3,3,1, 2,3,3,1, 1,3,3,2, 0,3,3,1, 2,3,3,2, 1,3,3,1, 2,3,3,1, 1,3,3,0, -1}},
        /*01100111*/{{-1},
                     {-1},
                     {-1}},
        /*01101000*/{{-1},
                     {-1},
                     {-1}},
        /*01101001*/{{-1},
                     {-1},
                     {-1}},
        /*01101010*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*01101011*/{{-1},
                     {-1},
                     {-1}},
        /*01101100*/{{-1},
                     {-1},
                     {-1}},
        /*01101101*/{{-1},
                     {-1},
                     {-1}},
        /*01101110*/{{-1},
                     {-1},
                     {-1}},
        /*01101111*/{{0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, -1}},
        /*01110000*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*01110001*/{{-1},
                     {-1},
                     {-1}},
        /*01110010*/{{-1},
                     {-1},
                     {-1}},
        /*01110011*/{{-1},
                     {-1},
                     {-1}},
        /*01110100*/{{-1},
                     {-1},
                     {-1}},
        /*01110101*/{{-1},
                     {-1},
                     {-1}},
        /*01110110*/{{-1},
                     {-1},
                     {-1}},
        /*01110111*/{{-1},
                     {-1},
                     {-1}},
        /*01111000*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*01111001*/{{-1},
                     {-1},
                     {-1}},
        /*01111010*/{{-1},
                     {-1},
                     {-1}},
        /*01111011*/{{-1},
                     {-1},
                     {-1}},
        /*01111100*/{{-1},
                     {-1},
                     {-1}},
        /*01111101*/{{-1},
                     {-1},
                     {-1}},
        /*01111110*/{{-1},
                     {-1},
                     {-1}},
        /*01111111*/{{-1},
                     {-1},
                     {-1}},
        /*10000000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10000001*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10000010*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10000011*/{{-1},
                     {-1},
                     {-1}},
        /*10000100*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10000101*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10000110*/{{-1},
                     {-1},
                     {-1}},
        /*10000111*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 0,0,3,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,2,2, 0,0,0,0, 2,3,3,2, 0,0,0,0, 2,2,3,3, 0,0,0,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 0,0,0,0, 2,2,2,2, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*10001000*/{{0,0,1,1, 0,0,1,1, 1,0,3,3, 1,0,3,3, 0,0,3,3, 1,1,3,3, 1,1,3,3, 1,1,3,3, 1,1,3,3, 0,0,3,3, -1},
                     {0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 0,0,0,0, 1,1,1,1, 1,1,1,1, 2,2,2,2, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,3,3,0, 0,2,2,0, 0,2,2,0, 0,2,2,0, 0,2,2,0, 0,3,3,0,-1}},
        /*10001001*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*10001010*/{{-1},
                     {-1},
                     {-1}},
        /*10001011*/{{-1},
                     {-1},
                     {-1}},
        /*10001100*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*10001101*/{{-1},
                     {-1},
                     {-1}},
        /*10001110*/{{-1},
                     {-1},
                     {-1}},
        /*10001111*/{{-1},
                     {-1},
                     {-1}},
        /*10010000*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,0,3,2, 1,0,3,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,3,3, 1,1,1,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 1,0,0,1, 3,3,3,3, -1},
                     {0,3,3,0, 0,2,2,0, 0,3,2,0, 0,3,2,0, 0,2,2,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 2,3,3,2, 2,3,3,2, -1}},
        /*10010001*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*10010010*/{{-1},
                     {-1},
                     {-1}},
        /*10010011*/{{-1},
                     {-1},
                     {-1}},
        /*10010100*/{{-1},
                     {-1},
                     {-1}},
        /*10010101*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*10010110*/{{-1},
                     {-1},
                     {-1}},
        /*10010111*/{{-1},
                     {-1},
                     {-1}},
        /*10011000*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*10011001*/{{0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0, 3,3,3,3, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 1,1,1,1, 2,2,2,2, 1,0,0,1, 2,3,3,2, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,2,0, 0,2,1,0, 0,2,2,0, 0,1,1,0, 0,2,3,0, 0,1,2,0, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, 0,2,1,0, 0,2,1,0, 0,1,1,0, 0,1,1,0, 0,1,2,0, 0,1,2,0, 2,3,2,1, 2,3,2,1, 1,2,3,2, 1,2,3,2, 0,2,1,0, 0,3,2,0, 0,1,1,0, 0,2,2,0, 0,1,2,0, 0,2,3,0, -1}},
        /*10011010*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*10011011*/{{-1},
                     {-1},
                     {-1}},
        /*10011100*/{{-1},
                     {-1},
                     {-1}},
        /*10011101*/{{-1},
                     {-1},
                     {-1}},
        /*10011110*/{{-1},
                     {-1},
                     {-1}},
        /*10011111*/{{0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, -1}},
        /*10100000*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10100001*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10100010*/{{-1},
                     {-1},
                     {-1}},
        /*10100011*/{{0,0,3,3, 0,0,2,2, 0,0,3,2, 0,0,3,2, 0,0,2,2, 0,0,3,3, 0,0,2,3, 0,0,2,3, 2,2,3,3, 2,2,3,3, 0,0,2,2, 0,0,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,0,0, 2,2,3,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,2,3,0, 1,2,3,0, 1,2,2,1, 1,2,2,1, -1}},
        /*10100100*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10100101*/{{0, 0, 3, 3, 0, 0, 3, 3, -1},
                     {0, 0, 0, 0, 3, 3, 3, 3, -1},
                     {0, 3, 3, 0, 0, 3, 3, 0, -1}},
        /*10100110*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,2,2, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,3,3,1, 0,3,3,1, 0,3,3,1, 1,3,3,0, 1,3,3,0, 1,3,3,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 1,3,3,1, 1,3,3,1, -1}},
        /*10100111*/{{-1},
                     {-1},
                     {-1}},
        /*10101000*/{{-1},
                     {-1},
                     {-1}},
        /*10101001*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,2,2,0, 0,2,3,0, 0,2,3,0, 0,2,2,0, 0,3,3,0, 0,3,2,0, 0,3,2,0, 2,3,3,2, 2,3,3,2, 0,2,2,0, 0,2,2,0, -1}},
        /*10101010*/{{-1},
                     {-1},
                     {-1}},
        /*10101011*/{{-1},
                     {-1},
                     {-1}},
        /*10101100*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*10101101*/{{-1},
                     {-1},
                     {-1}},
        /*10101110*/{{-1},
                     {-1},
                     {-1}},
        /*10101111*/{{-1},
                     {-1},
                     {-1}},
        /*10110000*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*10110001*/{{-1},
                     {-1},
                     {-1}},
        /*10110010*/{{-1},
                     {-1},
                     {-1}},
        /*10110011*/{{-1},
                     {-1},
                     {-1}},
        /*10110100*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*10110101*/{{-1},
                     {-1},
                     {-1}},
        /*10110110*/{{-1},
                     {-1},
                     {-1}},
        /*10110111*/{{-1},
                     {-1},
                     {-1}},
        /*10111000*/{{-1},
                     {-1},
                     {-1}},
        /*10111001*/{{-1},
                     {-1},
                     {-1}},
        /*10111010*/{{-1},
                     {-1},
                     {-1}},
        /*10111011*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, -1}},
        /*10111100*/{{-1},
                     {-1},
                     {-1}},
        /*10111101*/{{-1},
                     {-1},
                     {-1}},
        /*10111110*/{{-1},
                     {-1},
                     {-1}},
        /*10111111*/{{-1},
                     {-1},
                     {-1}},
        /*11000000*/{{0,0,3,3, 1,1,3,3, 0,0,1,1, 0,0,1,1, 1,0,3,3, 1,0,3,3, 1,1,3,3, 1,1,3,3, 0,1,3,3, 0,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,1,1, 3,3,3,3, 0,1,1,0, 3,3,3,3, -1},
                     {0,3,3,0, 1,2,2,1, 0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, 0,1,1,0, 0,1,1,0, -1}},
        /*11000001*/{{-1},
                     {-1},
                     {-1}},
        /*11000010*/{{-1},
                     {-1},
                     {-1}},
        /*11000011*/{{-1},
                     {-1},
                     {-1}},
        /*11000100*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11000101*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11000110*/{{-1},
                     {-1},
                     {-1}},
        /*11000111*/{{-1},
                     {-1},
                     {-1}},
        /*11001000*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11001001*/{{-1},
                     {-1},
                     {-1}},
        /*11001010*/{{0,0,3,3, 1,1,3,3, 1,0,3,3, 1,0,3,3, 1,1,3,3, 0,0,3,3, 0,1,3,3, 0,1,3,3, 0,0,1,1, 0,0,1,1, 1,1,3,3, 1,1,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 1,0,0,1, 2,3,3,2, 2,2,2,2, 3,3,3,3, 0,1,1,0, 3,2,2,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 2,3,3,2, 2,3,3,2, 1,2,2,1, 0,3,3,0, 0,1,1,0, 0,1,1,0, 0,3,2,1, 0,3,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11001011*/{{-1},
                     {-1},
                     {-1}},
        /*11001100*/{{0, 0, 1, 1, 1, 1, 2, 2, 1, 0, 3, 3, 2, 1, 3, 3, 1, 1, 3, 3, 2, 2, 3, 3, 0, 1, 3, 3, 1, 2, 3, 3, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 1, 3, 3, 2, 1, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 1, 2, 3, 3, 1, 2, 3, 3, 2, 1, 3, 3, 1, 0, 3, 3, 2, 2, 3, 3, 1, 1, 3, 3, 1, 2, 3, 3, 0, 1, 3, 3, 1, 1, 2, 2, 0, 0, 1, 1, -1},
                     {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, -1},
                     {0, 3, 2, 1, 0, 3, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 1, 2, 2, 1, 1, 2, 2, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 3, 3, 0, 0, 3, 3, 0, 0, 3, 2, 1, 0, 3, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 1, 2, 2, 1, 1, 2, 2, 1, 0, 1, 1, 0, 0, 1, 1, 0, 2, 3, 3, 2, 2, 3, 3, 2, 1, 2, 2, 1, 1, 2, 2, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, 0, 3, 2, 1, -1}},
        /*11001101*/{{-1},
                     {-1},
                     {-1}},
        /*11001110*/{{-1},
                     {-1},
                     {-1}},
        /*11001111*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, -1}},
        /*11010000*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11010001*/{{-1},
                     {-1},
                     {-1}},
        /*11010010*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11010011*/{{-1},
                     {-1},
                     {-1}},
        /*11010100*/{{-1},
                     {-1},
                     {-1}},
        /*11010101*/{{-1},
                     {-1},
                     {-1}},
        /*11010110*/{{-1},
                     {-1},
                     {-1}},
        /*11010111*/{{-1},
                     {-1},
                     {-1}},
        /*11011000*/{{-1},
                     {-1},
                     {-1}},
        /*11011001*/{{-1},
                     {-1},
                     {-1}},
        /*11011010*/{{-1},
                     {-1},
                     {-1}},
        /*11011011*/{{-1},
                     {-1},
                     {-1}},
        /*11011100*/{{-1},
                     {-1},
                     {-1}},
        /*11011101*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, -1}},
        /*11011110*/{{-1},
                     {-1},
                     {-1}},
        /*11011111*/{{-1},
                     {-1},
                     {-1}},
        /*11100000*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11100001*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 0,0,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,1,1, 3,3,3,3, 1,0,0,1, 3,3,3,3, 1,1,0,0, 3,3,3,3, 0,1,1,0, 3,3,3,3, 0,0,0,0, 1,1,1,1, 1,1,1,1, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, -1}},
        /*11100010*/{{-1},
                     {-1},
                     {-1}},
        /*11100011*/{{-1},
                     {-1},
                     {-1}},
        /*11100100*/{{-1},
                     {-1},
                     {-1}},
        /*11100101*/{{-1},
                     {-1},
                     {-1}},
        /*11100110*/{{-1},
                     {-1},
                     {-1}},
        /*11100111*/{{-1},
                     {-1},
                     {-1}},
        /*11101000*/{{-1},
                     {-1},
                     {-1}},
        /*11101001*/{{-1},
                     {-1},
                     {-1}},
        /*11101010*/{{-1},
                     {-1},
                     {-1}},
        /*11101011*/{{-1},
                     {-1},
                     {-1}},
        /*11101100*/{{-1},
                     {-1},
                     {-1}},
        /*11101101*/{{-1},
                     {-1},
                     {-1}},
        /*11101110*/{{0,0,1,1, 0,0,1,1, 1,0,3,2, 1,0,3,2, 2,2,3,3, 2,2,3,3, 0,1,2,3, 0,1,2,3, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, 0,0,0,0, 3,3,3,3, -1},
                     {0,3,2,1, 0,3,2,1, 2,3,3,2, 2,3,3,2, 1,2,3,0, 1,2,3,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, -1}},
        /*11101111*/{{-1},
                     {-1},
                     {-1}},
        /*11110000*/{{0,0,3,3, 1,1,2,2, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 0,0,1,1, 0,0,1,1, 0,0,1,1, 0,0,1,1, 0,0,1,1, 0,0,1,1, 2,2,3,3, 2,2,3,3, 2,2,3,3, 2,2,3,3, 1,1,2,2, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 1,1,2,2, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 0,0,1,1, 1,1,2,2, 1,1,1,1, 2,2,2,2, 1,1,0,0, 2,2,1,1, 1,1,2,2, 3,3,3,3, 1,0,1,2, 3,3,3,3, 0,1,2,1, 3,3,3,3, 2,1,0,1, 3,3,3,3, 2,2,1,1, 3,3,3,3, 2,1,1,2, 3,3,3,3, 1,2,1,0, 3,3,3,3, 1,2,2,1, 3,3,3,3, 2,2,2,2, 3,3,3,3, -1},
                     {0,3,3,0, 0,3,3,0, 0,3,3,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 0,3,3,0, 1,2,2,1, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 0,1,1,0, 0,1,1,0, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, 2,3,3,2, 2,3,3,2, 0,1,1,0, 0,1,1,0, 0,1,1,0, 0,1,1,0, 1,2,2,1, 1,2,2,1, -1}},
        /*11110001*/{{-1},
                     {-1},
                     {-1}},
        /*11110010*/{{-1},
                     {-1},
                     {-1}},
        /*11110011*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, -1}},
        /*11110100*/{{-1},
                     {-1},
                     {-1}},
        /*11110101*/{{-1},
                     {-1},
                     {-1}},
        /*11110110*/{{0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, -1}},
        /*11110111*/{{-1},
                     {-1},
                     {-1}},
        /*11111000*/{{-1},
                     {-1},
                     {-1}},
        /*11111001*/{{0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, 0,0,3,3, -1},
                     {0,0,0,0, 1,1,1,1, 0,1,1,0, 3,2,2,3, 2,2,2,2, 3,3,3,3, 1,0,0,1, 2,3,3,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 1,2,2,1, 0,1,1,0, 0,1,1,0, 1,2,2,1, 0,3,3,0, 2,3,3,2, 2,3,3,2, 1,2,2,1, 1,2,2,1, -1}},
        /*11111010*/{{-1},
                     {-1},
                     {-1}},
        /*11111011*/{{-1},
                     {-1},
                     {-1}},
        /*11111100*/{{0,0,3,3, 1,1,2,2, 2,2,3,3, 2,2,3,3, 1,1,2,2, 0,0,3,3, 0,0,1,1, 0,0,1,1, 1,1,2,2, 1,1,2,2, -1},
                     {0,0,0,0, 1,1,1,1, 1,1,0,0, 2,2,3,3, 2,2,2,2, 3,3,3,3, 0,0,1,1, 3,3,2,2, 1,1,1,1, 2,2,2,2, -1},
                     {0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, 0,3,3,0, -1}},
        /*11111101*/{{-1},
                     {-1},
                     {-1}},
        /*11111110*/{{-1},
                     {-1},
                     {-1}},
        /*11111111*/{{-1},
                     {-1},
                     {-1}}
    };

    const int subCodesTable [][105] = {
        /*00000000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00000111*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001011*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001101*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001110*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00001111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010011*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00010111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011001*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011110*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00011111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100011*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100110*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00100111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101101*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00101111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110001*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110010*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110101*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00110111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111010*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*00111111*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000110*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01000111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001011*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001100*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01001111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010011*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010110*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01010111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011001*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011100*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01011111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100010*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100100*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100101*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01100111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101010*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01101111*/{1,2,4,8, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110000*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01110111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//{1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111000*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*01111111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10000111*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001001*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001100*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10001111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010001*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010101*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10010111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011000*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011010*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10011111*/{1,2,4,8, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100011*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100110*/{ 1,2,4,8, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10100111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101001*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101100*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10101111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110000*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110100*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10110111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111011*/{1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*10111111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000100*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000101*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11000111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001000*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001010*/{ 1,2,4,8, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11001111*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010000*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010010*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11010111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011101*/{1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11011111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100000*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100001*/{ 1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11100111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101110*/{1,2,16,32, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11101111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110001*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110011*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110100*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110110*/{1,2,4,8, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11110111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111000*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111001*/{1,2,4,8, -1, -1, -1, -1, 1,8,16,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 2,4,32,64, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111010*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111011*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111100*/{ 1,2,4,8, -1, -1, -1, -1, 4,8,64,128, -1, -1, -1, -1, 16,32,64,128, -1, -1, -1, -1, 1,2,16,32, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111101*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111110*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        /*11111111*/{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
    };
//    for (int i = 0; i < 256; i++)
//    {
//        std::cout << "V" << i << " <--" << "{";
//        for (int j = 0; HexTable[i][0][j] != -1; j += 8)
//        {
//            std::cout << "(v" << HexTable[i][0][j] << HexTable[i][1][j] << HexTable[i][2][j] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 1] << HexTable[i][1][j + 1] << HexTable[i][2][j + 1] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 2] << HexTable[i][1][j + 2] << HexTable[i][2][j + 2] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 3] << HexTable[i][1][j + 3] << HexTable[i][2][j + 3] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 4] << HexTable[i][1][j + 4] << HexTable[i][2][j + 4] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 5] << HexTable[i][1][j + 5] << HexTable[i][2][j + 5] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 6] << HexTable[i][1][j + 6] << HexTable[i][2][j + 6] << ", ";
//            std::cout << "v" << HexTable[i][0][j + 7] << HexTable[i][1][j + 7] << HexTable[i][2][j + 7] << "), ";
//        }
//        if (HexTable[i][0][0] == -1)
//        {
//            for (int ii = 0; ii < 3; ii++)
//                for (int j = 0; j < 3; j++)
//                    for (int k = 0; k < 3; k++)
//                    {
//                        std::cout << "(v" << ii << j << k << ", ";
//                        std::cout << "v" << ii << j << k + 1 << ", ";
//                        std::cout << "v" << ii + 1 << j << k + 1 << ", ";
//                        std::cout << "v" << ii + 1 << j << k << ", ";
//                        std::cout << "v" << ii << j + 1 << k << ", ";
//                        std::cout << "v" << ii << j + 1 << k + 1 << ", ";
//                        std::cout << "v" << ii + 1 << j + 1 << k + 1 << ", ";
//                        std::cout << "v" << ii + 1 << j + 1 << k << "), ";
//                    }
//        }
//        std::cout << std::endl;
//    }

    AdjacentSet refined_nodes;
    AdjacentSet refined_elements;
    AdjacentSet additional_elemetns;;
    // build a set of nodes that connect refined elements
    for (UInteger elnum: eNumbers)
    {
        Hexahedral element = element_[elnum];
        for (int i = 0; i < 8; i++)
        {
            AdjacentSet a = node_[element[i]].adjacent;
            refined_elements.insert(a.begin(), a.end());
            refined_nodes.insert(element[i]);
        }
    }
    bool is_exists = true;
    while (is_exists)
    {
        is_exists = false;
        for (UInteger elnum: refined_elements)
        {
            Hexahedral element = element_[elnum];
            int code = 0;
            if (refined_nodes.find(element[0]) != refined_nodes.end())
                code |= 1;
            if (refined_nodes.find(element[1]) != refined_nodes.end())
                code |= 2;
            if (refined_nodes.find(element[2]) != refined_nodes.end())
                code |= 4;
            if (refined_nodes.find(element[3]) != refined_nodes.end())
                code |= 8;
            if (refined_nodes.find(element[4]) != refined_nodes.end())
                code |= 16;
            if (refined_nodes.find(element[5]) != refined_nodes.end())
                code |= 32;
            if (refined_nodes.find(element[6]) != refined_nodes.end())
                code |= 64;
            if (refined_nodes.find(element[7]) != refined_nodes.end())
                code |= 128;
            if (code != 0 && code != 255 && (HexTable[code][0][0] == -1 /*|| subCodesTable[code][0] != -1*/))
            {
                for (int i = 0; i < 8; i++)
                {
                    AdjacentSet a = node_[element[i]].adjacent;
                    refined_elements.insert(a.begin(), a.end());
                    refined_nodes.insert(element[i]);
                    is_exists = true;
                }
            }
        }
    }
    do
    {
        additional_elemetns.clear();
        ConsoleProgress progress(refined_elements.size());
        for (UInteger elnum: refined_elements)
        {
            ++progress;
            Hexahedral element = element_[elnum];
            Node3D n[8];
            Node3D local_nodes[4][4][4];
            for (int i = 0; i < 8; i++)
            {
                n[i] = node_[element[i]];
            }
            int code = 0;
            if (refined_nodes.find(element[0]) != refined_nodes.end())
                code |= 1;
            if (refined_nodes.find(element[1]) != refined_nodes.end())
                code |= 2;
            if (refined_nodes.find(element[2]) != refined_nodes.end())
                code |= 4;
            if (refined_nodes.find(element[3]) != refined_nodes.end())
                code |= 8;
            if (refined_nodes.find(element[4]) != refined_nodes.end())
                code |= 16;
            if (refined_nodes.find(element[5]) != refined_nodes.end())
                code |= 32;
            if (refined_nodes.find(element[6]) != refined_nodes.end())
                code |= 64;
            if (refined_nodes.find(element[7]) != refined_nodes.end())
                code |= 128;

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                    {
                        local_nodes[i][j][k].type = INNER;
                    }

            local_nodes[0][0][0] = n[0];
            local_nodes[0][0][3] = n[1];
            local_nodes[3][0][3] = n[2];
            local_nodes[3][0][0] = n[3];
            local_nodes[0][3][0] = n[4];
            local_nodes[0][3][3] = n[5];
            local_nodes[3][3][3] = n[6];
            local_nodes[3][3][0] = n[7];

            const double step = 1.0 / 3.0;
            for (int i = 1; i < 3; i++)
            {
                // ребра: "нижняя" грань
                local_nodes[0][0][i].point = local_nodes[0][0][0].point + static_cast<double>(i) * step * (local_nodes[0][0][3].point - local_nodes[0][0][0].point);
                local_nodes[i][0][3].point = local_nodes[0][0][3].point + static_cast<double>(i) * step * (local_nodes[3][0][3].point - local_nodes[0][0][3].point);
                local_nodes[3][0][i].point = local_nodes[3][0][0].point + static_cast<double>(i) * step * (local_nodes[3][0][3].point - local_nodes[3][0][0].point);
                local_nodes[i][0][0].point = local_nodes[0][0][0].point + static_cast<double>(i) * step * (local_nodes[3][0][0].point - local_nodes[0][0][0].point);
                // ребра: "верхняя" грань
                local_nodes[0][3][i].point = local_nodes[0][3][0].point + static_cast<double>(i) * step * (local_nodes[0][3][3].point - local_nodes[0][3][0].point);
                local_nodes[i][3][3].point = local_nodes[0][3][3].point + static_cast<double>(i) * step * (local_nodes[3][3][3].point - local_nodes[0][3][3].point);
                local_nodes[3][3][i].point = local_nodes[3][3][0].point + static_cast<double>(i) * step * (local_nodes[3][3][3].point - local_nodes[3][3][0].point);
                local_nodes[i][3][0].point = local_nodes[0][3][0].point + static_cast<double>(i) * step * (local_nodes[3][3][0].point - local_nodes[0][3][0].point);
                // "боковые" ребра
                local_nodes[0][i][0].point = local_nodes[0][0][0].point + static_cast<double>(i) * step * (local_nodes[0][3][0].point - local_nodes[0][0][0].point);
                local_nodes[0][i][3].point = local_nodes[0][0][3].point + static_cast<double>(i) * step * (local_nodes[0][3][3].point - local_nodes[0][0][3].point);
                local_nodes[3][i][3].point = local_nodes[3][0][3].point + static_cast<double>(i) * step * (local_nodes[3][3][3].point - local_nodes[3][0][3].point);
                local_nodes[3][i][0].point = local_nodes[3][0][0].point + static_cast<double>(i) * step * (local_nodes[3][3][0].point - local_nodes[3][0][0].point);
            }
            // узлы на гранях
            for (int i = 1; i < 3; i++)
            {
                for (int j = 1; j < 3; j++)
                {
                    // низ
                    local_nodes[j][0][i].point = local_nodes[0][0][i].point + static_cast<double>(j) * step * (local_nodes[3][0][i].point - local_nodes[0][0][i].point);
                    // верх
                    local_nodes[j][3][i].point = local_nodes[0][3][i].point + static_cast<double>(j) * step * (local_nodes[3][3][i].point - local_nodes[0][3][i].point);
                    // лево
                    local_nodes[0][j][i].point = local_nodes[0][0][i].point + static_cast<double>(j) * step * (local_nodes[0][3][i].point - local_nodes[0][0][i].point);
                    // право
                    local_nodes[3][j][i].point = local_nodes[3][0][i].point + static_cast<double>(j) * step * (local_nodes[3][3][i].point - local_nodes[3][0][i].point);
                    // зад
                    local_nodes[j][i][0].point = local_nodes[0][i][0].point + static_cast<double>(j) * step * (local_nodes[3][i][0].point - local_nodes[0][i][0].point);
                    // перед
                    local_nodes[j][i][3].point = local_nodes[0][i][3].point + static_cast<double>(j) * step * (local_nodes[3][i][3].point - local_nodes[0][i][3].point);
                }
            }

            for (int i = 1; i < 3; i++)
            {
                for (int j = 1; j < 3; j++)
                {
                    for (int k = 1; k < 3; k++)
                    {
                        local_nodes[i][k][j].point = local_nodes[i][0][j].point + static_cast<double>(k) * step * (local_nodes[i][3][j].point - local_nodes[i][0][j].point);
                    }
                }
            }

            // низ.
            if(n[0].type == BORDER && n[1].type == BORDER && n[2].type == BORDER && n[3].type == BORDER)
            {
                for (int i = 0; i <= 3; i++)
                    for (int j = 0; j <= 3; j++)
                    {
                        local_nodes[i][0][j].type = BORDER;
//                        if (func != nullptr)
//                        {
//                            local_nodes[i][0][j].point = findBorder(local_nodes[i][0][j].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
//                        }
                    }
            }
            // верх
            if(n[4].type == BORDER && n[5].type == BORDER && n[6].type == BORDER && n[7].type == BORDER)
            {
                for (int i = 0; i <= 3; i++)
                    for (int j = 0; j <= 3; j++)
                    {
                        local_nodes[i][3][j].type = BORDER;
//                        if (func != nullptr)
//                        {
//                            local_nodes[i][3][j].point = findBorder(local_nodes[i][3][j].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
//                        }
                    }
            }
            // лево
            if(n[0].type == BORDER && n[1].type == BORDER && n[4].type == BORDER && n[5].type == BORDER)
            {
                for (int i = 0; i <= 3; i++)
                    for (int j = 0; j <= 3; j++)
                    {
                        local_nodes[0][i][j].type = BORDER;
//                        if (func != nullptr)
//                        {
//                            local_nodes[0][i][j].point = findBorder(local_nodes[0][i][j].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
//                        }
                    }
            }
            // право
            if(n[2].type == BORDER && n[3].type == BORDER && n[6].type == BORDER && n[7].type == BORDER)
            {
                for (int i = 0; i <= 3; i++)
                    for (int j = 0; j <= 3; j++)
                    {
                        local_nodes[3][i][j].type = BORDER;
//                        if (func != nullptr)
//                        {
//                            local_nodes[3][i][j].point = findBorder(local_nodes[3][i][j].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
//                        }
                    }
            }
            // зад
            if(n[0].type == BORDER && n[3].type == BORDER && n[4].type == BORDER && n[7].type == BORDER)
            {
                for (int i = 0; i <= 3; i++)
                    for (int j = 0; j <= 3; j++)
                    {
                        local_nodes[i][j][0].type = BORDER;
//                        if (func != nullptr)
//                        {
//                            local_nodes[i][j][0].point = findBorder(local_nodes[i][j][0].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
//                        }
                    }
            }
            // перед
            if(n[1].type == BORDER && n[2].type == BORDER && n[5].type == BORDER && n[6].type == BORDER)
            {
                for (int i = 0; i <= 3; i++)
                    for (int j = 0; j <= 3; j++)
                    {
                        local_nodes[i][j][3].type = BORDER;
//                        if (func != nullptr)
//                        {
//                            local_nodes[i][j][3].point = findBorder(local_nodes[i][j][3].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
//                        }
                    }
            }
            // узлы на ребрах
            if (local_nodes[0][0][0].type == BORDER && local_nodes[0][0][3].type == BORDER)
            {
                local_nodes[0][0][1].type = BORDER;
                local_nodes[0][0][2].type = BORDER;
            }
            if (local_nodes[0][0][3].type == BORDER && local_nodes[3][0][3].type == BORDER)
            {
                local_nodes[1][0][3].type = BORDER;
                local_nodes[2][0][3].type = BORDER;
            }
            if (local_nodes[3][0][3].type == BORDER && local_nodes[3][0][0].type == BORDER)
            {
                local_nodes[3][0][1].type = BORDER;
                local_nodes[3][0][2].type = BORDER;
            }
            if (local_nodes[0][0][0].type == BORDER && local_nodes[3][0][0].type == BORDER)
            {
                local_nodes[1][0][0].type = BORDER;
                local_nodes[2][0][0].type = BORDER;
            }
            if (local_nodes[0][3][0].type == BORDER && local_nodes[0][3][3].type == BORDER)
            {
                local_nodes[0][3][1].type = BORDER;
                local_nodes[0][3][2].type = BORDER;
            }
            if (local_nodes[0][3][3].type == BORDER && local_nodes[3][3][3].type == BORDER)
            {
                local_nodes[1][3][3].type = BORDER;
                local_nodes[2][3][3].type = BORDER;
            }
            if (local_nodes[3][3][3].type == BORDER && local_nodes[3][3][0].type == BORDER)
            {
                local_nodes[3][3][1].type = BORDER;
                local_nodes[3][3][2].type = BORDER;
            }
            if (local_nodes[0][3][0].type == BORDER && local_nodes[3][3][0].type == BORDER)
            {
                local_nodes[1][3][0].type = BORDER;
                local_nodes[2][3][0].type = BORDER;
            }
            if (local_nodes[0][0][0].type == BORDER && local_nodes[0][3][0].type == BORDER)
            {
                local_nodes[0][1][0].type = BORDER;
                local_nodes[0][2][0].type = BORDER;
            }
            if (local_nodes[0][0][3].type == BORDER && local_nodes[0][3][3].type == BORDER)
            {
                local_nodes[0][1][3].type = BORDER;
                local_nodes[0][2][3].type = BORDER;
            }
            if (local_nodes[3][0][3].type == BORDER && local_nodes[3][3][3].type == BORDER)
            {
                local_nodes[3][1][3].type = BORDER;
                local_nodes[3][2][3].type = BORDER;
            }
            if (local_nodes[3][0][0].type == BORDER && local_nodes[3][3][0].type == BORDER)
            {
                local_nodes[3][1][0].type = BORDER;
                local_nodes[3][2][0].type = BORDER;
            }

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                    {
                        if (local_nodes[i][j][k].type == BORDER)
                            local_nodes[i][j][k].point = findBorder(local_nodes[i][j][k].point, func, 0.5 * (n[0].point.distanceTo(n[6].point) + n[1].point.distanceTo(n[7].point)));
                    }

            if (code != 0)
            {
                node_[element[0]].adjacent.erase(elnum);
                node_[element[1]].adjacent.erase(elnum);
                node_[element[2]].adjacent.erase(elnum);
                node_[element[3]].adjacent.erase(elnum);
                node_[element[4]].adjacent.erase(elnum);
                node_[element[5]].adjacent.erase(elnum);
                node_[element[6]].adjacent.erase(elnum);
                node_[element[7]].adjacent.erase(elnum);
                if (code == 255)
                {
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            for (int k = 0; k < 3; k++)
                            {
                                element[0] = addNode(local_nodes[i][j][k]);
                                element[1] = addNode(local_nodes[i][j][k + 1]);
                                element[2] = addNode(local_nodes[i + 1][j][k + 1]);
                                element[3] = addNode(local_nodes[i + 1][j][k]);
                                element[4] = addNode(local_nodes[i][j + 1][k]);
                                element[5] = addNode(local_nodes[i][j + 1][k + 1]);
                                element[6] = addNode(local_nodes[i + 1][j + 1][k + 1]);
                                element[7] = addNode(local_nodes[i + 1][j + 1][k]);
                                if (i == 0 && j == 0 && k == 0)
                                {
                                    element_[elnum] = element;
                                    node_[element[0]].adjacent.insert(elnum);
                                    node_[element[1]].adjacent.insert(elnum);
                                    node_[element[2]].adjacent.insert(elnum);
                                    node_[element[3]].adjacent.insert(elnum);
                                    node_[element[4]].adjacent.insert(elnum);
                                    node_[element[5]].adjacent.insert(elnum);
                                    node_[element[6]].adjacent.insert(elnum);
                                    node_[element[7]].adjacent.insert(elnum);
                                }
                                else
                                {
                                    addElement(element);
                                }
                            }
                }
                else
                {
                    for (int i = 0; HexTable[code][0][i] != -1; i += 8)
                    {
                        for (int j = 0; j < 8; j++)
                        {
                            element[j] = addNode(local_nodes[HexTable[code][0][i + j]][HexTable[code][1][i + j]][HexTable[code][2][i + j]]);
                        }
                        if (i == 0)
                        {
                            element_[elnum] = element;
                            node_[element[0]].adjacent.insert(elnum);
                            node_[element[1]].adjacent.insert(elnum);
                            node_[element[2]].adjacent.insert(elnum);
                            node_[element[3]].adjacent.insert(elnum);
                            node_[element[4]].adjacent.insert(elnum);
                            node_[element[5]].adjacent.insert(elnum);
                            node_[element[6]].adjacent.insert(elnum);
                            node_[element[7]].adjacent.insert(elnum);
                            if (subCodesTable[code][i] != -1)
                            {
                                additional_elemetns.insert(elnum);
                            }
                        }
                        else
                        {
                            addElement(element);
                            if (subCodesTable[code][i] != -1)
                            {
                                additional_elemetns.insert(elementsCount() - 1);
                            }
                        }
                    }
                }
            }
        }
        refined_elements = additional_elemetns;
    }
    while (!additional_elemetns.empty());
}

}


