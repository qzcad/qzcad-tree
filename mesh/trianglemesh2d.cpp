#include "trianglemesh2d.h"
#include <iostream>
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#include <float.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace msh
{

TriangleMesh2D::TriangleMesh2D()
{
    xMin_ = -1.0;
    xMax_ = 1.0;
    yMin_ = -1.0;
    yMax_ = 1.0;
}

TriangleMesh2D::TriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height)
{
    double hx = width / (double)(xCount - 1);
    double hy = height / (double)(yCount - 1);
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double xCenter = (xMax_ + xMin_) / 2.0;
    const double yCenter = (yMax_ + yMin_) / 2.0;
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        double x = xMin + (double) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double y = yMin + (double) j * hy;
            Point2D point(x, y);

            if ((i == 0 && j == 0) || (i == 0 && j == yCount - 1) || (i == xCount - 1 && j == 0) || (i == xCount - 1 && j == yCount - 1))
                pushNode(point, CHARACTER);
            else if (i == 0 || j == 0 || i == xCount - 1 || j == yCount - 1)
                pushNode(point, BORDER);
            else
                pushNode(point, INNER);
        }
    }
    // формирование массива элементов
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            // Для симметричности сетки необходимо смена направления диагоналей в зависимости от четверти области
            Point2D p0 = node_[i * yCount + j].point;
            Point2D p1 = node_[(i + 1) * yCount + j].point;
            Point2D p2 = node_[(i + 1) * yCount + j + 1].point;
            Point2D p3 = node_[i * yCount + j + 1].point;
            if ((xCenter <= p0.x() && yCenter <= p0.y() && xCenter <= p1.x() && yCenter <= p1.y() && xCenter <= p2.x() && yCenter <= p2.y() && xCenter <= p3.x() && yCenter <= p3.y())
                    ||
                    (xCenter >= p0.x() && yCenter >= p0.y() && xCenter >= p1.x() && yCenter >= p1.y() && xCenter >= p2.x() && yCenter >= p2.y() && xCenter >= p3.x() && yCenter >= p3.y()))
            {
                addElement(i * yCount + j, (i + 1) * yCount + j, i * yCount + j + 1);
                addElement((i + 1) * yCount + j, (i + 1) * yCount + j + 1, i * yCount + j + 1);
            } // if
            else
            {
                addElement(i * yCount + j, (i + 1) * yCount + j, (i + 1) * yCount + j + 1);
                addElement(i * yCount + j, (i + 1) * yCount + j + 1, i * yCount + j + 1);
            } // else
        }
    }

    std::cout << "Создана равномерная сетка треугольных элементов: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

TriangleMesh2D::TriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint)
{
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    const double minDistance = 0.4 * sqrt(hx*hx + hy*hy);
    std::map<UInteger, UInteger> nodesMap;
    // формирование массива узлов
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for (UInteger i = 0; i < xCount; i++)
    {
        double x = xMin + (double) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double y = yMin + (double) j * hy;
            Point2D point(x, y);

            if (func(x, y) >= 0.0)
            {
#ifdef WITH_OPENMP
#pragma omp critical
#endif
                nodesMap[i * yCount + j] = pushNode(point, INNER);
            }
        }
    }

    // формирование начальной сетки
    const double xCenter = (xMax_ + xMin_) / 2.0;
    const double yCenter = (yMax_ + yMin_) / 2.0;
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            std::map<UInteger, UInteger>::iterator iter0 = nodesMap.find(i * yCount + j);
            std::map<UInteger, UInteger>::iterator iter1 = nodesMap.find((i + 1) * yCount + j);
            std::map<UInteger, UInteger>::iterator iter2 = nodesMap.find((i + 1) * yCount + j + 1);
            std::map<UInteger, UInteger>::iterator iter3 = nodesMap.find(i * yCount + j + 1);
            if (iter0 != nodesMap.end() && iter1 != nodesMap.end() && iter2 != nodesMap.end() && iter3 != nodesMap.end())
            {
                // Для симметричности сетки необходимо смена направления диагоналей в зависимости от четверти области
                Point2D p0 = node_[iter0->second].point;
                Point2D p1 = node_[iter1->second].point;
                Point2D p2 = node_[iter2->second].point;
                Point2D p3 = node_[iter3->second].point;
                if ((xCenter <= p0.x() && yCenter <= p0.y() && xCenter <= p1.x() && yCenter <= p1.y() && xCenter <= p2.x() && yCenter <= p2.y() && xCenter <= p3.x() && yCenter <= p3.y())
                        ||
                        (xCenter >= p0.x() && yCenter >= p0.y() && xCenter >= p1.x() && yCenter >= p1.y() && xCenter >= p2.x() && yCenter >= p2.y() && xCenter >= p3.x() && yCenter >= p3.y()))
                {
#ifdef WITH_OPENMP
#pragma omp critical
                    {
#endif
                        addElement(iter0->second, iter1->second, iter3->second);
                        addElement(iter1->second, iter2->second, iter3->second);
#ifdef WITH_OPENMP
                    } // omp critical
#endif
                } // if
                else
                {
#ifdef WITH_OPENMP
#pragma omp critical
                    {
#endif
                    addElement(iter0->second, iter1->second, iter2->second);
                    addElement(iter0->second, iter2->second, iter3->second);
#ifdef WITH_OPENMP
                    } // omp critical
#endif
                } // else
            }
            else if (iter1 != nodesMap.end() && iter2 != nodesMap.end() && iter3 != nodesMap.end())
            {
                addElement(iter1->second, iter2->second, iter3->second);
            }
            else if (iter0 != nodesMap.end() && iter2 != nodesMap.end() && iter3 != nodesMap.end())
            {
                addElement(iter0->second, iter2->second, iter3->second);
            }
            else if (iter0 != nodesMap.end() && iter1 != nodesMap.end() && iter3 != nodesMap.end())
            {
                addElement(iter0->second, iter1->second, iter3->second);
            }
            else if (iter0 != nodesMap.end() && iter1 != nodesMap.end() && iter2 != nodesMap.end())
            {
                addElement(iter0->second, iter1->second, iter2->second);
            }
        }
    }

    std::vector<Point2D> normal(nodesCount()); // нормали к узлам
    std::vector<UInteger> iso(nodesCount()); // изо-точки (номера)
    UInteger baseElementCount = elementsCount();
    // построение нормалей для всех узлов начальной сетки
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for (UInteger i = 0; i < normal.size(); i++)
    {
        Point2D currentPoint = node_[i].point;
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
                prevPoint = node_[tri[2]].point;
                nextPoint = node_[tri[1]].point;
            }
            else if (tri[1] == i)
            {
                prevPoint = node_[tri[0]].point;
                nextPoint = node_[tri[2]].point;
            }
            else // if (tri[2] == i)
            {
                prevPoint = node_[tri[1]].point;
                nextPoint = node_[tri[0]].point;
            }
            prevNormal.set((currentPoint.y() - prevPoint.y()), -(currentPoint.x() - prevPoint.x()));
            nextNormal.set((nextPoint.y() - currentPoint.y()), -(nextPoint.x() - currentPoint.x()));
            normal[i] = normal[i] + prevNormal.normalized() + nextNormal.normalized();
        }
    }
    // поиск изо-точек на границе
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
    for (UInteger i = 0; i < normal.size(); i++)
    {
        iso[i] = ULONG_MAX;
        if (normal[i].length() > epsilon_)
        {
            Point2D n = normal[i].normalized();
            Point2D current = node_[i].point;
            // двоичный поиск граничной точки
            Point2D inner = current;
            Point2D outer = current + sqrt(hx*hx + hy*hy) * n;
            Point2D mid;
            double val = 0.0;
            if (func(outer.x(), outer.y()) > 0.0)
            {
                // внешняя точка перескачила через границу и попала внутрь
                // ищем внешнюю точку пошагово сканированием
                outer = current;
                while (func(outer.x(), outer.y()) >= 0.0)
                    outer = outer + (0.0001 * sqrt(hx*hx + hy*hy)) * n;
            }
            do
            {
                mid = 0.5 * (inner + outer);
                val = func(mid.x(), mid.y());
                if (val <= 0.0)
                    outer = mid;
                else //if (val > 0.0)
                    inner = mid;
            } while(fabs(val) > epsilon_);

            // поиск соответствующей характерной точки
            for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
            {
                if (mid.distanceTo(*cPoint) < minDistance)
                {
                    mid = *cPoint;
                    //charPoint.erase(cPoint);
                    break;
                }
            }

            // сравнение с существующими изо-точками перед всатвкой
            bool isExist = false;
            for (UInteger j = 0; j < i; j++)
            {
                if (iso[j] < ULONG_MAX)
                {
                    Point2D border = node_[iso[j]].point;
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
#ifdef WITH_OPENMP
#pragma omp critical
#endif
                    iso[i] = pushNode(mid, BORDER);
                }
            }
        }
    }
    // для характерных точек, которым не нашлась пара на этапе формирования нормалей, используем метод близжайшего узла
    for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
    {
        std::vector<Node2D>::iterator min_n;
        double min_d = 10.0 * minDistance;
        for (std::vector<Node2D>::iterator n = node_.begin(); n != node_.end(); ++n)
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
                    Point2D p0 = node_[triangle[j + 1]].point;
                    Point2D p1 = node_[triangle[j]].point;
                    Point2D p2 = node_[iso[triangle[j]]].point;
                    Point2D p3 = node_[iso[triangle[j+ 1]]].point;
                    double ma1 = minAngle(p0, p1, p3);
                    double mb1 = minAngle(p1, p2, p3);
                    double ma2 = minAngle(p0, p1, p2);
                    double mb2 = minAngle(p0, p2, p3);
                    if (std::min(ma1, mb1) > std::min(ma2, mb2))
                    {
                        if (triangle[j + 1] != iso[triangle[j + 1]] && triangle[j] != iso[triangle[j + 1]])
                            addElement(triangle[j + 1], triangle[j], iso[triangle[j + 1]]);
                        if (triangle[j] != iso[triangle[j]] && iso[triangle[j]] != iso[triangle[j + 1]])
                            addElement(triangle[j], iso[triangle[j]], iso[triangle[j + 1]]);
                    } // if
                    else
                    {
                        if (triangle[j] != iso[triangle[j]] && triangle[j + 1] != iso[triangle[j]])
                            addElement(triangle[j + 1], triangle[j], iso[triangle[j]]);
                        if (triangle[j + 1] != iso[triangle[j + 1]] && iso[triangle[j]] != iso[triangle[j + 1]])
                            addElement(triangle[j + 1], iso[triangle[j]], iso[triangle[j + 1]]);
                    } // else
                }
            }
        }
    }

    // the end.
    std::cout << "Создана сетка треугольных элементов для функционального объекта: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

TriangleMesh2D::TriangleMesh2D(const TriangleMesh2D &mesh)
{
    element_ = mesh.element_;
    node_ = mesh.node_;
    xMin_ = mesh.xMin_;
    xMax_ = mesh.xMax_;
    yMin_ = mesh.yMin_;
    yMax_ = mesh.yMax_;
}

TriangleMesh2D::TriangleMesh2D(const TriangleMesh2D *mesh)
{
    element_ = mesh->element_;
    node_ = mesh->node_;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
}

TriangleMesh2D::TriangleMesh2D(const SegmentMesh2D *mesh)
{
    Triangulation triangulation = SuperDelaunay(mesh);
    for (UInteger i = 0; i < triangulation.nodes.size(); i++) pushNode(triangulation.nodes[i], BORDER);
    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
         triangle != triangulation.triangles.end(); ++triangle)
    {
        addElement(triangle->vertexNode(0), triangle->vertexNode(1), triangle->vertexNode(2));
        if (jacobian(elementsCount() - 1) < 0)
        {
            UInteger t = element_[elementsCount() - 1][0];
            element_[elementsCount() - 1][0] = element_[elementsCount() - 1][1];
            element_[elementsCount() - 1][1] = t;
        }
    }
    xMin_ = mesh->xMin(); xMax_ = mesh->xMax();
    yMin_ = mesh->yMin(); yMax_ = mesh->yMax();
}

UInteger TriangleMesh2D::elementsCount() const
{
    return element_.size();
}

ElementPointer TriangleMesh2D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

bool TriangleMesh2D::isBorderElement(const UInteger &number) const
{
    for (int i = 0; i < 3; i++)
    {
        if (node_[element_[number].vertexNode(i)].type == BORDER || node_[element_[number].vertexNode(i)].type == CHARACTER)
            return true;
    }
    return false;
}

void TriangleMesh2D::directionChange()
{
    for (UInteger i = 0; i < elementsCount(); i++)
    {
        std::swap(element_[i][1], element_[i][2]); // изменение порядка обхода: 0-1-2 <= 0-2-1
    }
}

double TriangleMesh2D::area(const UInteger &number) const
{
    Triangle triangle = element_[number];
    Point2D a = node_[triangle[0]].point;
    Point2D b = node_[triangle[1]].point;
    Point2D c = node_[triangle[2]].point;
    //           | xa ya 1 |
    // S = 0.5 * | xb yb 1 |
    //           | xc yc 1 |
    return 0.5 * fabs( (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()) );
}

void TriangleMesh2D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2)
{
    Triangle triangle(node0, node1, node2);
    element_.push_back(triangle);
    // обновление списка смежных узлов
    node_[node0].adjacent.insert(element_.size() - 1);
    node_[node1].adjacent.insert(element_.size() - 1);
    node_[node2].adjacent.insert(element_.size() - 1);
}

double TriangleMesh2D::jacobian(const UInteger &elementNum)
{
    const Triangle tri = element_[elementNum];
    const Point2D p0 = node_[tri[0]].point;
    const Point2D p1 = node_[tri[1]].point;
    const Point2D p2 = node_[tri[2]].point;
    double dNdXi[3] = {
        -1.0,
        1.0,
        0.0
    };
    double dNdEta[3] = {
        -1.0,
        0.0,
        1.0
    };
    double j[2][2];
    j[0][0] = p0.x() * dNdXi[0] + p1.x() * dNdXi[1] + p2.x() * dNdXi[2];    j[0][1] = p0.y() * dNdXi[0] + p1.y() * dNdXi[1] + p2.y() * dNdXi[2];
    j[1][0] = p0.x() * dNdEta[0] + p1.x() * dNdEta[1] + p2.x() * dNdEta[2]; j[1][1] = p0.y() * dNdEta[0] + p1.y() * dNdEta[1] + p2.y() * dNdEta[2];
    return j[0][0] * j[1][1] - j[0][1] * j[1][0];
}

double TriangleMesh2D::lengthAspect(const UInteger &elNum)
{
    const Triangle tri = element_[elNum];
    const Point2D p0 = node_[tri[0]].point;
    const Point2D p1 = node_[tri[1]].point;
    const Point2D p2 = node_[tri[2]].point;
    const double d01 = p0.distanceTo(p1);
    const double d12 = p1.distanceTo(p2);
    const double d20 = p2.distanceTo(p0);
    double min = std::min(d01, std::min(d12, d20));
    double max = std::max(d01, std::max(d12, d20));
    return min / max;
}

double TriangleMesh2D::minAngle(const UInteger &elNum)
{
    const Triangle tri = element_[elNum];
    const Point2D p0 = node_[tri[0]].point;
    const Point2D p1 = node_[tri[1]].point;
    const Point2D p2 = node_[tri[2]].point;
    return minAngle(p0, p1, p2);
}

double TriangleMesh2D::angleAspect(const UInteger &elNum)
{
    const Triangle tri = element_[elNum];
    const Point2D p0 = node_[tri[0]].point;
    const Point2D p1 = node_[tri[1]].point;
    const Point2D p2 = node_[tri[2]].point;
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    double min, max;
    if (angles(p0, p1, p2, alpha, beta, gamma))
    {
        // Треугольник невырожденный
        min = std::min(alpha, std::min(beta, gamma));
        max = std::max(alpha, std::max(beta, gamma));
        return min / max;
    }
    return 0.0;
}

void TriangleMesh2D::delaunay(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint)
{
    SegmentMesh2D mesh(xCount, yCount, xMin, yMin, width, height, func, charPoint);
    Triangulation triangulation = SuperDelaunay(&mesh);
    node_.clear();
    element_.clear();
    for (UInteger i = 0; i < triangulation.nodes.size(); i++) pushNode(triangulation.nodes[i], BORDER);
    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
         triangle != triangulation.triangles.end(); ++triangle)
    {
        Point2D p0 = triangulation.nodes[triangle->vertexNode(0)];
        Point2D p1 = triangulation.nodes[triangle->vertexNode(1)];
        Point2D p2 = triangulation.nodes[triangle->vertexNode(2)];
        Point2D c = (1.0 / 3.0) * (p0 + p1 + p2);
        if (func(c.x(), c.y()) > 0.0)
        {
            addElement(triangle->vertexNode(0), triangle->vertexNode(1), triangle->vertexNode(2));
            if (jacobian(elementsCount() - 1) < 0)
            {
                UInteger t = element_[elementsCount() - 1][0];
                element_[elementsCount() - 1][0] = element_[elementsCount() - 1][1];
                element_[elementsCount() - 1][1] = t;
            }
        }
    }
    xMin_ = mesh.xMin(); xMax_ = mesh.xMax();
    yMin_ = mesh.yMin(); yMax_ = mesh.yMax();
}

void TriangleMesh2D::ruppert(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint)
{
    SegmentMesh2D mesh(xCount, yCount, xMin, yMin, width, height, func, charPoint);
    Triangulation triangulation = SuperDelaunay(&mesh);
    node_.clear();
    element_.clear();
    // удаление треугольников, не попавших в область
    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
         triangle != triangulation.triangles.end(); )
    {
        Point2D p0 = triangulation.nodes[triangle->vertexNode(0)];
        Point2D p1 = triangulation.nodes[triangle->vertexNode(1)];
        Point2D p2 = triangulation.nodes[triangle->vertexNode(2)];
        Point2D c = (1.0 / 3.0) * (p0 + p1 + p2);
        if (func(c.x(), c.y()) < 0.0)
            triangle = triangulation.triangles.erase(triangle);
        else
            ++triangle;
    }
    // обработка ребер
//    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
//         triangle != triangulation.triangles.end(); ++triangle) std::cout << (*triangle)[0] << " " << (*triangle)[1] << " " << (*triangle)[2] << "\n";
    int iii = 0;
//    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
//    while (triangle != triangulation.triangles.end() && iii < 1000)
//    {
//        Segment e0(triangle->vertexNode(0), triangle->vertexNode(1));
//        Segment e1(triangle->vertexNode(1), triangle->vertexNode(2));
//        Segment e2(triangle->vertexNode(2), triangle->vertexNode(0));
//        Segment segments[3] = {e0, e1, e2};
//        bool isRefined = false;
//        std::cout << iii << ":    ";
//        for (int i = 0; i < 3; i++)
//        {
//            bool isEdge = false;
//            for (UInteger j = 0; j < mesh.elementsCount(); j++)
//            {
//                Segment edge = mesh.segment(j);
//                if (edge.isSame(segments[i]))
//                {
//                    isEdge = true;
//                    break;
//                }
//            }
//            if (!isEdge)
//            {
//                Point2D A = triangulation.nodes[segments[i][0]];
//                Point2D B = triangulation.nodes[segments[i][1]];
//                double r = 0.5 * A.distanceTo(B);
//                Point2D center = 0.5 * (A + B);
//                bool shellRefine = false;
//                for (UInteger j = 0; j < triangulation.nodes.size(); j++)
//                {
//                    if (j != segments[i][0] && j != segments[i][1])
//                    {
//                        Point2D p = triangulation.nodes[j];
//                        if ((r*r - (p.x() - center.x())*(p.x() - center.x()) - (p.y() - center.y())*(p.y() - center.y())) > 0.0)
//                        {
//                            shellRefine = true;
//                            std::cout << A.x() << " " << A.y() << " -- " << B.x() << " " << B.y() << " -- " << center.x() << " " << center.y() << "\n";
//                            break;
//                        }
//                    }
//                }
//                if (shellRefine)
//                {
//                    insertDelaunayNode(center, triangulation.nodes, triangulation.triangles);
//                    isRefined = true;
//                    triangle = triangulation.triangles.begin();
////                    continue;
//                }
//            }
//        }

//        ++iii;
//        if (!isRefined) ++triangle;
//    }

    iii = 0;
    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
    while (triangle != triangulation.triangles.end() && iii < 1000)
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        if (minAngle(A, B, C) > 0.1 * M_PI)
        {
            Point2D center = circumcenter(A, B, C);
            std::cout << center.x() << " " << center.y() << "\n";
//            if (func(center.x(), center.y()) > 0)
//            {
                if (insertDelaunayNode(center, triangulation.nodes, triangulation.triangles))
                    triangle = triangulation.triangles.begin();
//            }
            else
                ++triangle;
        }
        else
            ++triangle;
        iii++;
    }

    for (UInteger i = 0; i < triangulation.nodes.size(); i++) pushNode(triangulation.nodes[i], BORDER);
    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
         triangle != triangulation.triangles.end(); ++triangle)
    {
       addElement(triangle->vertexNode(0), triangle->vertexNode(1), triangle->vertexNode(2));
       if (jacobian(elementsCount() - 1) < 0)
       {
           UInteger t = element_[elementsCount() - 1][0];
           element_[elementsCount() - 1][0] = element_[elementsCount() - 1][1];
           element_[elementsCount() - 1][1] = t;
       }
    }
    xMin_ = mesh.xMin(); xMax_ = mesh.xMax();
    yMin_ = mesh.yMin(); yMax_ = mesh.yMax();
}


double TriangleMesh2D::minAngle(const Point2D &A, const Point2D &B, const Point2D &C)
{
    double alpha, beta, gamma;
    angles(A, B, C, alpha, beta, gamma);
    return std::min(alpha, std::min(beta, gamma));
}

bool TriangleMesh2D::angles(const Point2D &A, const Point2D &B, const Point2D &C, double &alpha, double &beta, double &gamma)
{
    const double a = B.distanceTo(C); // сторона, противолежащяя вершине A (BC)
    const double b = A.distanceTo(C); // сторона, противолежащяя вершине B (AC)
    const double c = A.distanceTo(B); // сторона, противолежащяя вершине C (AB)
    if (a < epsilon_ || b < epsilon_ || c < epsilon_)
    {
        alpha = beta = gamma = 0.0;
        return false;
    }
    // Теорема косинусов
    alpha = acos((b*b + c*c - a*a) / (2.0 * b * c)); // Угол в вершине A
    // Теорема синусов
    // обеспечение вычислительной устойчисвости: значение может на бесконечно малую отклоняться от единицы для угла 90 градусов
    double sinBeta = sin(alpha) * b / a;
    beta = (sinBeta > 1.0) ? M_PI_2 : asin(sinBeta); // Угол в вершине B
    // Теорема о сумме углов треугольника
    gamma = M_PI - (alpha + beta); // Угол в вершине C
    return true;
}

bool TriangleMesh2D::inCircumcircle(const Point2D &A, const Point2D &B, const Point2D &C, const Point2D &p)
{
    auto det3 = [](double m[3][3])
    {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    };
    double ma[3][3] = {
        { A.x(), A.y(), 1.0 },
        { B.x(), B.y(), 1.0 },
        { C.x(), C.y(), 1.0 }
    };
    double a = det3(ma);
    double mb[3][3] = {
        { A.x()*A.x() + A.y()*A.y(), A.y(), 1.0 },
        { B.x()*B.x() + B.y()*B.y(), B.y(), 1.0 },
        { C.x()*C.x() + C.y()*C.y(), C.y(), 1.0 }
    };
    double b = det3(mb);
    double mc[3][3] = {
        { A.x()*A.x() + A.y()*A.y(), A.x(), 1.0 },
        { B.x()*B.x() + B.y()*B.y(), B.x(), 1.0 },
        { C.x()*C.x() + C.y()*C.y(), C.x(), 1.0 }
    };
    double c = det3(mc);
    double md[3][3] = {
        { A.x()*A.x() + A.y()*A.y(), A.x(), A.y() },
        { B.x()*B.x() + B.y()*B.y(), B.x(), B.y() },
        { C.x()*C.x() + C.y()*C.y(), C.x(), C.y() }
    };
    double d = det3(md);
    double v = a * (p.x()*p.x() + p.y()*p.y()) - b * p.x() + c * p.y() - d;

    return (v * (a > 0.0 ? 1.0 : -1.0) <= 0.0);
}

Point2D TriangleMesh2D::circumcenter(const Point2D &A, const Point2D &B, const Point2D &C)
{
    auto det3 = [](double m[3][3])
    {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    };
    double ma[3][3] = {
        { A.x(), A.y(), 1.0 },
        { B.x(), B.y(), 1.0 },
        { C.x(), C.y(), 1.0 }
    };
    double a = det3(ma);
    double mb[3][3] = {
        { A.x()*A.x() + A.y()*A.y(), A.y(), 1.0 },
        { B.x()*B.x() + B.y()*B.y(), B.y(), 1.0 },
        { C.x()*C.x() + C.y()*C.y(), C.y(), 1.0 }
    };
    double b = det3(mb);
    double mc[3][3] = {
        { A.x()*A.x() + A.y()*A.y(), A.x(), 1.0 },
        { B.x()*B.x() + B.y()*B.y(), B.x(), 1.0 },
        { C.x()*C.x() + C.y()*C.y(), C.x(), 1.0 }
    };
    double c = det3(mc);

    return Point2D(0.5 * b / a, -0.5 * c / a);
}

TriangleMesh2D::Triangulation TriangleMesh2D::SuperDelaunay(const SegmentMesh2D *mesh)
{
    double super_width = mesh->xMax() - mesh->xMin();
    double super_height = mesh->yMax() - mesh->yMin();
    Point2D super0(mesh->xMin() - 0.5 * super_width, mesh->yMin() - 0.5 * super_height);
    Point2D super1(mesh->xMax() + 0.5 * super_width, mesh->yMin() - 0.5 * super_height);
    Point2D super2(mesh->xMax() + 0.5 * super_width, mesh->yMax() + 0.5 * super_height);
    Point2D super3(mesh->xMin() - 0.5 * super_width, mesh->yMax() + 0.5 * super_height);
    Triangulation triangulation;
    std::vector<Point2D> nodes;
    std::list<Triangle> triangles;
    nodes.push_back(super0);
    nodes.push_back(super1);
    nodes.push_back(super2);
    nodes.push_back(super3);
    nodes.push_back(mesh->point2d(0));
    triangles.push_back(Triangle(0, 1, 4));
    triangles.push_back(Triangle(1, 2, 4));
    triangles.push_back(Triangle(2, 3, 4));
    triangles.push_back(Triangle(3, 0, 4));
    for (UInteger i = 1; i < mesh->nodesCount(); i++)
    {
        Point2D point = mesh->point2d(i);
        insertDelaunayNode(point, nodes, triangles);
    }
//    triangulation.nodes = nodes; triangulation.triangles = triangles;
    for (UInteger i = 4; i < nodes.size(); i++) triangulation.nodes.push_back(nodes[i]);
    for (std::list<Triangle>::iterator triangle = triangles.begin(); triangle != triangles.end(); ++triangle)
    {
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
            triangulation.triangles.push_back(Triangle(triangle->vertexNode(0) - 4UL,
                                                       triangle->vertexNode(1) - 4UL,
                                                       triangle->vertexNode(2) - 4UL));
        }
    }
    return triangulation;
}

bool TriangleMesh2D::insertDelaunayNode(const Point2D &point, std::vector<Point2D> &nodes, std::list<Triangle> &triangles)
{
    UInteger number = nodes.size();
    std::vector<int> power;
    std::vector<Segment> edges;
    for (std::vector<Point2D>::iterator p = nodes.begin(); p != nodes.end(); ++p)
        if (p->isEqualTo(point, epsilon_)) return false;
    nodes.push_back(point);
    for (std::list<Triangle>::iterator triangle = triangles.begin(); triangle != triangles.end(); )
    {
        Point2D A = nodes[triangle->vertexNode(0)];
        Point2D B = nodes[triangle->vertexNode(1)];
        Point2D C = nodes[triangle->vertexNode(2)];
        Segment e0(triangle->vertexNode(0), triangle->vertexNode(1));
        Segment e1(triangle->vertexNode(1), triangle->vertexNode(2));
        Segment e2(triangle->vertexNode(2), triangle->vertexNode(0));
        bool flags[] = {false, false, false};
        if (inCircumcircle(A, B, C, point))
        {
            for (UInteger j = 0; j < edges.size(); j++)
            {
                if (e0.isSame(edges[j]))
                {
                    power[j]++;
                    flags[0] = true;
                }
                else if (e1.isSame(edges[j]))
                {
                    power[j]++;
                    flags[1] = true;
                }
                else if (e2.isSame(edges[j]))
                {
                    power[j]++;
                    flags[2] = true;
                }
            }

            if (!flags[0])
            {
                edges.push_back(e0);
                power.push_back(1);
            }
            if (!flags[1])
            {
                edges.push_back(e1);
                power.push_back(1);
            }
            if (!flags[2])
            {
                edges.push_back(e2);
                power.push_back(1);
            }
            triangle = triangles.erase(triangle);
        }
        else
        {
            ++triangle;
        }
    }

    for (UInteger j = 0; j < edges.size(); j++)
    {
        if (power[j] == 1)
        {
            triangles.push_back(Triangle(edges[j].vertexNode(0), edges[j].vertexNode(1), number));
        }
    }
    return true;
}

}
