#include "trianglemesh3d.h"

#include <math.h>

namespace msh {

TriangleMesh3D::TriangleMesh3D()
{
    xMin_ = -1.0;
    xMax_ = 1.0;
    yMin_ = -1.0;
    yMax_ = 1.0;
    zMin_ = -1.0;
    zMax_ = 1.0;
}

TriangleMesh3D::TriangleMesh3D(const TriangleMesh3D &mesh)
{
    node_ = mesh.node_;
    element_ = mesh.element_;
    xMin_ = mesh.xMin_;
    xMax_ = mesh.xMax_;
    yMin_ = mesh.yMin_;
    yMax_ = mesh.yMax_;
    zMin_ = mesh.zMin_;
    zMax_ = mesh.zMax_;
}

TriangleMesh3D::TriangleMesh3D(const TriangleMesh3D *mesh)
{
    node_ = mesh->node_;
    element_ = mesh->element_;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
    zMin_ = mesh->zMin_;
    zMax_ = mesh->zMax_;
}

TriangleMesh3D::TriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length)
{
    double hphi = 2.0 * M_PI / (double)rCount;
    double hl = length / (double)lCount;
    double phi = 0.0;
    // формирование массива узлов
    for (UInteger i = 0; i < rCount; i++)
    {
        double l = 0.0;
        for (UInteger j = 0; j <= lCount; j++)
        {
            Point3D point(radius * cos(phi), l, radius * sin(phi));
            if (j == 0 || j == lCount)
                pushNode(point, CHARACTER);
            else
                pushNode(point, BORDER);
            l += hl;
        }
        phi += hphi;
    }
    // формирование массива элементов
    for (UInteger i = 0; i < rCount; i++)
    {
        for (UInteger j = 0; j < lCount; j++)
        {
            if (i < rCount - 1)
            {
                addElement(i * (lCount + 1) + j, (i + 1) * (lCount + 1) + j, (i + 1) * (lCount + 1) + j + 1);
                addElement(i * (lCount + 1) + j, (i + 1) * (lCount + 1) + j + 1, i * (lCount + 1) + j + 1);
            }
            else
            {
                addElement(i * (lCount + 1) + j, j, j + 1);
                addElement(i * (lCount + 1) + j, j + 1, i * (lCount + 1) + j + 1);
            }
        }
    }
    // размеры области
    xMin_ = zMin_ = -radius;
    xMax_ = zMax_ = radius;
    yMin_ = 0.0;
    yMax_ = length;
}

UInteger TriangleMesh3D::elementsCount() const
{
    return element_.size();
}

ElementPointer TriangleMesh3D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

bool TriangleMesh3D::isBorderElement(const UInteger &number) const
{
    for (int i = 0; i < 3; i++)
    {
        if (node_[element_[number].vertexNode(i)].type == BORDER || node_[element_[number].vertexNode(i)].type == CHARACTER)
            return true;
    }
    return false;
}

double TriangleMesh3D::surfaceArea() const
{
    double s = 0.0;
    for (UInteger i = 0; i < elementsCount(); i++)
        s += area(i);
    return s;
}

void TriangleMesh3D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2)
{
    Triangle triangle(node0, node1, node2);
    addElement(triangle);
}

void TriangleMesh3D::addElement(const Triangle &triangle)
{
    element_.push_back(triangle);
    // обновление списка смежных узлов
    node_[triangle[0]].adjacent.insert(element_.size() - 1);
    node_[triangle[1]].adjacent.insert(element_.size() - 1);
    node_[triangle[2]].adjacent.insert(element_.size() - 1);
}

double TriangleMesh3D::area(const UInteger &number) const
{
    Triangle triangle = element_[number];
    Point3D p0 = node_[triangle[0]].point;
    Point3D p1 = node_[triangle[1]].point;
    Point3D p2 = node_[triangle[2]].point;
    double a = p0.distanceTo(p1);
    double b = p1.distanceTo(p2);
    double c = p2.distanceTo(p0);
    double p = (a + b + c) / 2.0;
    return sqrt(p * (p - a) * (p - b) * (p - c)); // формула Герона
}

}