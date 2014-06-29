#include "hexahedralmesh3d.h"
#include <math.h>

namespace msh {
HexahedralMesh3D::HexahedralMesh3D()
{
}

HexahedralMesh3D::HexahedralMesh3D(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const Floating &xMin, const Floating &yMin, const Floating &zMin, const Floating &width, const Floating &height, const Floating &depth)
{
    Floating hx = width / (Floating)(xCount - 1);
    Floating hy = height / (Floating)(yCount - 1);
    Floating hz = depth / (Floating)(zCount - 1);
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        Floating x = xMin + (Floating) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            Floating y = yMin + (Floating) j * hy;
            for (UInteger k = 0; k < zCount; k++)
            {
                Floating z = zMin + (Floating)k * hz;
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
    //    minimizeFunctional();
}

HexahedralMesh3D::HexahedralMesh3D(QuadrilateralMesh2D *baseMesh, const Floating &xDelta, const Floating &yDelta, const int &lCount, bool x_axes)
{
    xMin_ = baseMesh->xMin();
    xMax_ = baseMesh->xMax();
    yMin_ = baseMesh->yMin();
    yMax_ = baseMesh->yMax();
    zMin_ = 0.0;
    zMax_ = 0.0;

    node_.reserve(baseMesh->nodesCount() *(lCount + 1));
    element_.reserve(baseMesh->elementsCount() * lCount);

    UInteger i;
    // текущий угол поворта
    Floating phi;
    Floating delta_phi = 2.0 * M_PI / (double) lCount;
    // номера вершин шестигранника
    UInteger nodes_pointers[8];
    UInteger baseNodesCount = baseMesh->nodesCount();

    // формирование узлов
    for(int iz = 0; iz < lCount; iz++)
    {
        phi = (Floating)iz * delta_phi;

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
    for(int iz = 0; iz < lCount; iz++)
    {
        for(i = 0; i < baseMesh->elementsCount(); i++)
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

            addElement(nodes_pointers[0], nodes_pointers[1], nodes_pointers[2], nodes_pointers[3], nodes_pointers[4], nodes_pointers[5], nodes_pointers[6], nodes_pointers[7]);
        }
    }
}

HexahedralMesh3D::HexahedralMesh3D(QuadrilateralMesh2D *baseMesh, const Floating &xDelta, const Floating &yDelta, const Floating &angle, const int &lCount, bool x_axes)
{
    xMin_ = baseMesh->xMin();
    xMax_ = baseMesh->xMax();
    yMin_ = baseMesh->yMin();
    yMax_ = baseMesh->yMax();
    zMin_ = 0.0;
    zMax_ = 0.0;

    node_.reserve(baseMesh->nodesCount() * lCount);
    element_.reserve(baseMesh->elementsCount() * lCount);

    UInteger i;
    // текущий угол поворта
    Floating phi;
    Floating delta_phi = angle * M_PI / ((Floating) lCount * 180.0);
    // номера вершин шестигранника
    UInteger nodes_pointers[8];
    UInteger baseNodesCount = baseMesh->nodesCount();

    // формирование узлов
    for(int iz = 0; iz < lCount; iz++)
    {
        phi = (Floating)iz * delta_phi;

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
    for(int iz = 0; iz < lCount - 1; iz++)
    {
        for(i = 0; i < baseMesh->elementsCount(); i++)
        {
            Quadrilateral current_quad = baseMesh->quadrilateral(i);
//            nodes_pointers[0] = current_quad[0] + iz * baseNodesCount;
//            nodes_pointers[1] = current_quad[1] + iz * baseNodesCount;
//            nodes_pointers[2] = current_quad[2] + iz * baseNodesCount;
//            nodes_pointers[3] = current_quad[3] + iz * baseNodesCount;
//            nodes_pointers[4] = current_quad[0] + (iz + 1) * baseNodesCount;
//            nodes_pointers[5] = current_quad[1] + (iz + 1) * baseNodesCount;
//            nodes_pointers[6] = current_quad[2] + (iz + 1) * baseNodesCount;
//            nodes_pointers[7] = current_quad[3] + (iz + 1) * baseNodesCount;

            nodes_pointers[1] = current_quad[0] + iz * baseNodesCount;
            nodes_pointers[2] = current_quad[1] + iz * baseNodesCount;
            nodes_pointers[6] = current_quad[2] + iz * baseNodesCount;
            nodes_pointers[5] = current_quad[3] + iz * baseNodesCount;
            nodes_pointers[0] = current_quad[0] + (iz + 1) * baseNodesCount;
            nodes_pointers[3] = current_quad[1] + (iz + 1) * baseNodesCount;
            nodes_pointers[7] = current_quad[2] + (iz + 1) * baseNodesCount;
            nodes_pointers[4] = current_quad[3] + (iz + 1) * baseNodesCount;

            addElement(nodes_pointers[0], nodes_pointers[1], nodes_pointers[2], nodes_pointers[3], nodes_pointers[4], nodes_pointers[5], nodes_pointers[6], nodes_pointers[7]);
        }
    }
}

UInteger HexahedralMesh3D::elementsCount() const
{
    return element_.size();
}

ElementPointer HexahedralMesh3D::element(const UInteger &number) const
{
    Hexahedral *hptr = new Hexahedral(element_[number]);
    ElementPointer elementPtr(hptr);
    return elementPtr;
}

bool HexahedralMesh3D::isBorderElement(const UInteger &number) const
{
    for (int i = 0; i < 8; i++)
    {
        if (node_[element_[number].vertexNode(i)].type == BORDER || node_[element_[number].vertexNode(i)].type == CHARACTER)
            return true;
    }
    return false;
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
}


