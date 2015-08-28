#include "quadrilateralmesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>

namespace msh {

QuadrilateralMesh3D::QuadrilateralMesh3D()
{
    xMin_ = -1.0;
    xMax_ = 1.0;
    yMin_ = -1.0;
    yMax_ = 1.0;
    zMin_ = -1.0;
    zMax_ = 1.0;
}

QuadrilateralMesh3D::QuadrilateralMesh3D(const QuadrilateralMesh3D &mesh)
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

QuadrilateralMesh3D::QuadrilateralMesh3D(const QuadrilateralMesh3D *mesh)
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

QuadrilateralMesh3D::QuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length)
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
                addElement(i * (lCount + 1) + j, (i + 1) * (lCount + 1) + j, (i + 1) * (lCount + 1) + j + 1, i * (lCount + 1) + j + 1);
            else
                addElement(i * (lCount + 1) + j, j, j + 1, i * (lCount + 1) + j + 1);
        }
    }
    // размеры области
    xMin_ = zMin_ = -radius;
    xMax_ = zMax_ = radius;
    yMin_ = 0.0;
    yMax_ = length;
}

QuadrilateralMesh3D::QuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length)
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
            double radius = bottom_radius + l / length * (top_radius - bottom_radius);
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
                addElement(i * (lCount + 1) + j, (i + 1) * (lCount + 1) + j, (i + 1) * (lCount + 1) + j + 1, i * (lCount + 1) + j + 1);
            else
                addElement(i * (lCount + 1) + j, j, j + 1, i * (lCount + 1) + j + 1);
        }
    }
    // размеры области
    xMin_ = zMin_ = -(bottom_radius > top_radius ? bottom_radius : top_radius);
    xMax_ = zMax_ = (bottom_radius > top_radius ? bottom_radius : top_radius);
    yMin_ = 0.0;
    yMax_ = length;
}

UInteger QuadrilateralMesh3D::elementsCount() const
{
    return element_.size();
}

ElementPointer QuadrilateralMesh3D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

bool QuadrilateralMesh3D::isBorderElement(const UInteger &number) const
{
    for (int i = 0; i < 4; i++)
    {
        if (node_[element_[number].vertexNode(i)].type == BORDER || node_[element_[number].vertexNode(i)].type == CHARACTER)
            return true;
    }
    return false;
}

double QuadrilateralMesh3D::surfaceArea() const
{
    double s = 0.0;
    for (UInteger i = 0; i < elementsCount(); i++)
        s += area(i);
    return s;
}

void QuadrilateralMesh3D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3)
{
    Quadrilateral quad(node0, node1, node2, node3);
    addElement(quad);
}

void QuadrilateralMesh3D::addElement(const Quadrilateral &quad)
{
    element_.push_back(quad);
    // обновление списка смежных узлов
    node_[quad[0]].adjacent.insert(element_.size() - 1);
    node_[quad[1]].adjacent.insert(element_.size() - 1);
    node_[quad[2]].adjacent.insert(element_.size() - 1);
    node_[quad[3]].adjacent.insert(element_.size() - 1);
}

double QuadrilateralMesh3D::area(const UInteger &number) const
{
    Quadrilateral quad = element_[number];
    Point3D p0 = node_[quad[0]].point;
    Point3D p1 = node_[quad[1]].point;
    Point3D p2 = node_[quad[2]].point;
    Point3D p3 = node_[quad[3]].point;
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

void QuadrilateralMesh3D::add(const QuadrilateralMesh3D *mesh)
{
    std::vector<UInteger> nodesPointers(mesh->nodesCount());
    for (UInteger i = 0; i < mesh->nodesCount(); i++)
    {
        nodesPointers[i] = addNode(mesh->node_[i]);
    }
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        Quadrilateral q = mesh->element_[i];
        addElement(nodesPointers[q[0]], nodesPointers[q[1]], nodesPointers[q[2]], nodesPointers[q[3]]);
    }
    updateDomain();
}

void QuadrilateralMesh3D::translate(const double &x, const double &y, const double &z)
{
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point3D p = node_[i].point;
        node_[i].point.set(p.x() + x, p.y() + y, p.z() + z);
    }
}
}
