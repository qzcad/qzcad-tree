#include "trianglemesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#include <iostream>

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

TriangleMesh3D::TriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, std::function<double (double, double, double)> func, std::list<Point3D> charPoints)
{
    double hphi = 2.0 * M_PI / (double)rCount;
    double hl = length / (double)lCount;
    Point3D diag0 (radius * cos(0), 0, radius * sin(0));
    Point3D diag1 (radius * cos(hphi), hl, radius * sin(hphi));
    const double minDistance = 0.4 * diag0.distanceTo(diag1);
    double phi = 0.0;
    std::map<UInteger, UInteger> nodesMap;
    std::vector<Point2D> localCoordinates;
    // формирование массива узлов
    for (UInteger i = 0; i < rCount; i++)
    {
        double l = 0.0;
        for (UInteger j = 0; j <= lCount; j++)
        {
            Point3D point(radius * cos(phi), l, radius * sin(phi));
            if (func(point.x(), point.y(), point.z()) >= 0.0)
            {
                nodesMap[i * (lCount + 1) + j] = pushNode(point, BORDER);
                localCoordinates.push_back(Point2D(l, phi));
            }
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
                std::map<UInteger, UInteger>::iterator iter0 = nodesMap.find(i * (lCount + 1) + j);
                std::map<UInteger, UInteger>::iterator iter1 = nodesMap.find((i + 1) * (lCount + 1) + j);
                std::map<UInteger, UInteger>::iterator iter2 = nodesMap.find((i + 1) * (lCount + 1) + j + 1);
                std::map<UInteger, UInteger>::iterator iter3 = nodesMap.find(i * (lCount + 1) + j + 1);
                if (iter0 != nodesMap.end() && iter1 != nodesMap.end() && iter2 != nodesMap.end() && iter3 != nodesMap.end())
                {
                    addElement(iter0->second, iter1->second, iter2->second);
                    addElement(iter0->second, iter2->second, iter3->second);
                }
            }
            else
            {
                std::map<UInteger, UInteger>::iterator iter0 = nodesMap.find(i * (lCount + 1) + j);
                std::map<UInteger, UInteger>::iterator iter1 = nodesMap.find(j);
                std::map<UInteger, UInteger>::iterator iter2 = nodesMap.find(j + 1);
                std::map<UInteger, UInteger>::iterator iter3 = nodesMap.find(i * (lCount + 1) + j + 1);
                if (iter0 != nodesMap.end() && iter1 != nodesMap.end() && iter2 != nodesMap.end() && iter3 != nodesMap.end())
                {
                    addElement(iter0->second, iter1->second, iter2->second);
                    addElement(iter0->second, iter2->second, iter3->second);
                }
            }
        }
    }
    std::vector<Point2D> normal(nodesCount()); // нормали к узлам
    std::vector<UInteger> iso(nodesCount()); // изо-точки (номера)
    UInteger baseElementCount = elementsCount();
    // построение нормалей для всех узлов начальной сетки
    for (UInteger i = 0; i < normal.size(); i++)
    {
        Point2D currentPoint = localCoordinates[i];
        AdjacentSet triangles = node_[i].adjacent;
        // обработка смежных элементов
        for (AdjacentSet::iterator it = triangles.begin(); it != triangles.end(); ++it)
        {
            Triangle tri = element_[*it];
            Point2D prevPoint;
            Point2D nextPoint;
            Point2D prevNormal;
            Point2D nextNormal;
            if (tri[0] == i)
            {
                prevPoint = localCoordinates[tri[2]];
                nextPoint = localCoordinates[tri[1]];
            }
            else if (tri[1] == i)
            {
                prevPoint = localCoordinates[tri[0]];
                nextPoint = localCoordinates[tri[2]];
            }
            else // if (tri[2] == i)
            {
                prevPoint = localCoordinates[tri[1]];
                nextPoint = localCoordinates[tri[0]];
            }
            prevNormal.set((currentPoint.y() - prevPoint.y()), -(currentPoint.x() - prevPoint.x()));
            nextNormal.set((nextPoint.y() - currentPoint.y()), -(nextPoint.x() - currentPoint.x()));
            normal[i] = normal[i] + prevNormal.normalized() + nextNormal.normalized();
        }
    }
    // поиск изо-точек на границе
    for (UInteger i = 0; i < normal.size(); i++)
    {
        iso[i] = ULONG_MAX;
        if (normal[i].length() > epsilon_)
        {
            Point2D n = normal[i].normalized();
//            std::cout << i << ": " << n.x() << " " << n.y() << std::endl;
            Point2D current = localCoordinates[i];
            Point2D local_outter = current - sqrt(hl*hl + hphi*hphi) * n;
            // двоичный поиск граничной точки
            Point3D inner (radius * cos(current.y()), current.x(), radius * sin(current.y()));
            Point3D outer (radius * cos(local_outter.y()), local_outter.x(), radius * sin(local_outter.y()));
            Point3D mid;
            double val = 0.0;
            if (func(outer.x(), outer.y(), outer.z()) > 0.0)
            {
                continue;
            }
            do
            {
                mid = 0.5 * (inner + outer);
                val = func(mid.x(), mid.y(), mid.z());
                if (val <= 0.0)
                    outer = mid;
                else //if (val > 0.0)
                    inner = mid;
            } while(fabs(val) > epsilon_);
            // поиск соответствующей характерной точки
            for (std::list<Point3D>::iterator cPoint = charPoints.begin(); cPoint != charPoints.end(); ++cPoint)
            {
                if (mid.distanceTo(*cPoint) < minDistance)
                {
                    mid = *cPoint;
                    break;
                }
            }
            // сравнение с существующими изо-точками перед всатвкой
            bool isExist = false;
            for (UInteger j = 0; j < i; j++)
            {
                if (iso[j] < ULONG_MAX)
                {
                    Point3D border = node_[iso[j]].point;
                    if (border.distanceTo(mid) < minDistance)
                    {
                        iso[i] = iso[j];
                        isExist = true;
                        break;
                    }
                }
            }
            if (!isExist)
            {
                if (current.distanceTo(mid) < minDistance)
                {
                    node_[i].point = mid;
                    node_[i].type = BORDER;
                    iso[i] = i;
                }
                else
                {
                    iso[i] = pushNode(mid, BORDER);
                }
            }
        }
    }
    // для характерных точек, которым не нашлась пара на этапе формирования нормалей, используем метод близжайшего узла
    for (std::list<Point3D>::iterator cPoint = charPoints.begin(); cPoint != charPoints.end(); ++cPoint)
    {
        std::vector<Node3D>::iterator min_n;
        double min_d = 10.0 * minDistance;
        for (std::vector<Node3D>::iterator n = node_.begin(); n != node_.end(); ++n)
        {
            if (n->type == BORDER)
            {
                double d = (n->point).distanceTo(*cPoint);
                if (d < min_d)
                {
                    min_d = d;
                    min_n = n;
                }
            }
        }
        if (min_n != node_.end())
        {
            min_n->point = *cPoint;
            min_n->type = BORDER;
        }
    }
    // Форимрование приграничного слоя элементов
    for (UInteger i = 0; i < baseElementCount; i++)
    {
        Triangle triangle = element_[i];
        for (int j = 0; j < 3; j++)
        {
            if (iso[triangle[j]] < ULONG_MAX && iso[triangle[j + 1]] < ULONG_MAX)
            {
                AdjacentSet adjacent = node_[triangle[j]].adjacent;
                int adjacnetCount = 0;
                // обработка смежных элементов
                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                {
                    Triangle atri = element_[*it];
                    if (atri.in(triangle[j]) && atri.in(triangle[j + 1]))
                        ++adjacnetCount;
                }
                if (adjacnetCount == 1)
                {
                    // из двух возможных способов гегенрации треугольников выбирается тот, для которого минимальные углы будут максимальными
//                    Point2D p0 = node_[triangle[j + 1]].point;
//                    Point2D p1 = node_[triangle[j]].point;
//                    Point2D p2 = node_[iso[triangle[j]]].point;
//                    Point2D p3 = node_[iso[triangle[j+ 1]]].point;
//                    double ma1 = minAngle(p0, p1, p3);
//                    double mb1 = minAngle(p1, p2, p3);
//                    double ma2 = minAngle(p0, p1, p2);
//                    double mb2 = minAngle(p0, p2, p3);
//                    if (std::min(ma1, mb1) > std::min(ma2, mb2))
//                    {
                        if (triangle[j + 1] != iso[triangle[j + 1]] && triangle[j] != iso[triangle[j + 1]])
                            addElement(triangle[j + 1], triangle[j], iso[triangle[j + 1]]);
                        if (triangle[j] != iso[triangle[j]] && iso[triangle[j]] != iso[triangle[j + 1]])
                            addElement(triangle[j], iso[triangle[j]], iso[triangle[j + 1]]);
//                    } // if
//                    else
//                    {
//                        if (triangle[j] != iso[triangle[j]] && triangle[j + 1] != iso[triangle[j]])
//                            addElement(triangle[j + 1], triangle[j], iso[triangle[j]]);
//                        if (triangle[j + 1] != iso[triangle[j + 1]] && iso[triangle[j]] != iso[triangle[j + 1]])
//                            addElement(triangle[j + 1], iso[triangle[j]], iso[triangle[j + 1]]);
//                    } // else
                }
            }
        }
    }
    // размеры области
    xMin_ = zMin_ = -radius;
    xMax_ = zMax_ = radius;
    yMin_ = 0.0;
    yMax_ = length;
}

TriangleMesh3D::TriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length)
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
    xMin_ = zMin_ = -(bottom_radius > top_radius ? bottom_radius : top_radius);
    xMax_ = zMax_ = (bottom_radius > top_radius ? bottom_radius : top_radius);
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
