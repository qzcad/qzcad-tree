#include "hexahedralmesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <iostream>
#include <float.h>

namespace msh {
HexahedralMesh3D::HexahedralMesh3D() : Mesh3D(NULL)
{
}

void HexahedralMesh3D::prismDomain(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth)
{
    clear();
    double hx = width / (double)(xCount - 1);
    double hy = height / (double)(yCount - 1);
    double hz = depth / (double)(zCount - 1);
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        double x = xMin + (double) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double y = yMin + (double) j * hy;
            for (UInteger k = 0; k < zCount; k++)
            {
                double z = zMin + (double)k * hz;
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
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    zMin_ = zMin;
    zMax_ = xMin + depth;
    std::cout << "Создана равномерная сетка шестигранных элементов: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
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
    double delta_phi = 2.0 * M_PI / (double) lCount;
    // номера вершин шестигранника
    UInteger nodes_pointers[8];
    UInteger baseNodesCount = baseMesh->nodesCount();

    // формирование узлов
    for(int iz = 0; iz < lCount; iz++)
    {
        phi = (double)iz * delta_phi;

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

            if (x_axes)
                addElement(nodes_pointers[0], nodes_pointers[1], nodes_pointers[2], nodes_pointers[3], nodes_pointers[4], nodes_pointers[5], nodes_pointers[6], nodes_pointers[7]);
            else
                addElement(nodes_pointers[1], nodes_pointers[2], nodes_pointers[6], nodes_pointers[5], nodes_pointers[0], nodes_pointers[3], nodes_pointers[7], nodes_pointers[4]);
            if (withLayersInfo && baseMesh->sizeOfLayers() == baseMesh->elementsCount())
                pushLayer(baseMesh->layer(i));
        }
    }
    std::cout << "Создана сетка шестигранных элементов вращением плоского профиля: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
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
    double delta_phi = angle * M_PI / ((double) lCount * 180.0);
    // номера вершин шестигранника
    UInteger nodes_pointers[8];
    UInteger baseNodesCount = baseMesh->nodesCount();

    // формирование узлов
    for(int iz = 0; iz <= lCount; iz++)
    {
        phi = (double)iz * delta_phi;

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
    std::cout << "Создана сетка шестигранных элементов вращением плоского профиля: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

HexahedralMesh3D::HexahedralMesh3D(const HexahedralMesh3D &mesh) : Mesh3D(&mesh)
{
    element_ = mesh.element_;
    node_ = mesh.node_;
}

HexahedralMesh3D::HexahedralMesh3D(const HexahedralMesh3D *mesh): Mesh3D(mesh)
{
    element_ = mesh->element_;
    node_ = mesh->node_;
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

double HexahedralMesh3D::faceArea(const UIntegerVector &face) const
{
    Point3D p0 = node_[face[0]].point;
    Point3D p1 = node_[face[1]].point;
    Point3D p2 = node_[face[2]].point;
    Point3D p3 = node_[face[3]].point;
    // стороны
    double a = p0.distanceTo(p1);
    double b = p1.distanceTo(p2);
    double c = p2.distanceTo(p3);
    double d = p3.distanceTo(p0);
    // диагонали
    double d1 = p0.distanceTo(p2);
    double d2 = p1.distanceTo(p3);
    // функция для вычисления квадрата числа (C++0x)
    auto sqr = [](double value) { return value * value; };
    return sqrt(4.0 * sqr(d1) * sqr(d2) - sqr(sqr(b) + sqr(d) - sqr(a) - sqr(c))) / 4.0;
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
    // формирование шестигранников
    for(int i = 0; i < zLayersCount; i++)
    {
        for(UInteger j = 0; j < baseMesh->elementsCount(); j++)
        {
            Quadrilateral quad = baseMesh->quadrilateral(j);
            std::vector<UInteger> nnode(8);
            nnode[0] = quad[0] + (UInteger)i * baseMesh->nodesCount();
            nnode[1] = quad[1] + (UInteger)i * baseMesh->nodesCount();
            nnode[2] = quad[2] + (UInteger)i * baseMesh->nodesCount();
            nnode[3] = quad[3] + (UInteger)i * baseMesh->nodesCount();
            nnode[4] = quad[0] + (UInteger)(i + 1) * baseMesh->nodesCount();
            nnode[5] = quad[1] + (UInteger)(i + 1) * baseMesh->nodesCount();
            nnode[6] = quad[2] + (UInteger)(i + 1) * baseMesh->nodesCount();
            nnode[7] = quad[3] + (UInteger)(i + 1) * baseMesh->nodesCount();
            addElement(nnode);
        }
    }
    updateDomain();
}

}


