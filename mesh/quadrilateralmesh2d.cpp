#include "quadrilateralmesh2d.h"
#include <iostream>
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "funcopt.h"
#include "rfunctions.h"

namespace msh
{
QuadrilateralMesh2D::QuadrilateralMesh2D() : Mesh2D(NULL)
{
}

void QuadrilateralMesh2D::rectangleDomain(const UInteger &xCount,
                                          const UInteger &yCount,
                                          const double &xMin,
                                          const double &yMin,
                                          const double &width,
                                          const double &height)


{
    clear();
    double hx = width / (double)(xCount - 1);
    double hy = height / (double)(yCount - 1);
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
            addElement(i * yCount + j, (i + 1) * yCount + j, (i + 1) * yCount + j + 1, i * yCount + j + 1);
        }
    }
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
//    minimizeFunctional();
    std::cout << "Создана равномерная сетка четырехугольных элементов: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

void QuadrilateralMesh2D::quadDomain(const UInteger &xCount, const UInteger &yCount,
                                         const Point2D &v0, const Point2D &v1,
                                         const Point2D &v2, const Point2D &v3)
{
    clear();
    double hx = 2.0 / (double)(xCount - 1);
    double hy = 2.0 / (double)(yCount - 1);
    for (UInteger i = 0; i < xCount; i++)
    {
        double xi = -1.0 + (double) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double eta = -1.0 + (double) j * hy;
            Point2D point = isoFunc(0, xi, eta) * v0  + isoFunc(1, xi, eta) * v1 + isoFunc(2, xi, eta) * v2 + isoFunc(3, xi, eta) * v3;

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
            addElement(i * yCount + j, (i + 1) * yCount + j, (i + 1) * yCount + j + 1, i * yCount + j + 1);
        }
    }
    xMin_ = std::min(std::min(v0.x(), v1.x()), std::min(v2.x(), v3.x()));
    xMax_ = std::max(std::max(v0.x(), v1.x()), std::max(v2.x(), v3.x()));
    yMin_ = std::min(std::min(v0.y(), v1.y()), std::min(v2.y(), v3.y()));
    yMax_ = std::max(std::max(v0.y(), v1.y()), std::max(v2.y(), v3.y()));
    std::cout << "Создана равномерная (изопараметрическое преобразование) сетка четырехугольных элементов: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

void QuadrilateralMesh2D::triangleDomain(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2)
{
    clear();
    Point2D center = (v0 + v1 + v2) / 3.0; // центр треугольника
    Point2D c01 = (v0 + v1) / 2.0; // центр стороны, соединяющей вершину 0 и 1
    Point2D c12 = (v1 + v2) / 2.0; // центр стороны, соединяющей вершину 1 и 2
    Point2D c20 = (v2 + v0) / 2.0; // центр стороны, соединяющей вершину 2 и 0
    UInteger sideCount = count / 2 + 1;
    double h = 2.0 / (double)(sideCount - 1); // шаг изо-сетки
    Point2D quads [][4] = {
        {c20, v0, c01, center},
        {c01, v1, c12, center},
        {c12, v2, c20, center}
    };
    for (int q = 0; q < 3; q++)
    {
        UInteger nodeNumber[sideCount * sideCount];
        for (UInteger i = 0; i < sideCount; i++)
        {
            double xi = -1.0 + (double) i * h;
            for (UInteger j = 0; j < sideCount; j++)
            {
                double eta = -1.0 + (double) j * h;
                Point2D point = isoFunc(0, xi, eta) * quads[q][0]  +
                        isoFunc(1, xi, eta) * quads[q][1] +
                        isoFunc(2, xi, eta) * quads[q][2] +
                        isoFunc(3, xi, eta) * quads[q][3];
                if (j == 0 && i == sideCount - 1)
                    nodeNumber[i * sideCount + j] = addNode(point, CHARACTER);
                else if (j == 0 || i == sideCount - 1)
                    nodeNumber[i * sideCount + j] = addNode(point, BORDER);
                else
                    nodeNumber[i * sideCount + j] = addNode(point, INNER);
            }
        }
        // формирование массива элементов
        for (UInteger i = 0; i < sideCount - 1; i++)
        {
            for (UInteger j = 0; j < sideCount - 1; j++)
            {
                addElement(nodeNumber[i * sideCount + j],
                        nodeNumber[(i + 1) * sideCount + j],
                        nodeNumber[(i + 1) * sideCount + j + 1],
                        nodeNumber[i * sideCount + j + 1]);
            }
        }
    }
    xMin_ = std::min(std::min(v0.x(), v1.x()), v2.x());
    xMax_ = std::max(std::max(v0.x(), v1.x()), v2.x());
    yMin_ = std::min(std::min(v0.y(), v1.y()), v2.y());
    yMax_ = std::max(std::max(v0.y(), v1.y()), v2.y());
    std::cout << "Создана равномерная сетка четырехугольных элементов для тругольной области: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

void QuadrilateralMesh2D::circleDomain(const UInteger &count, const Point2D &center, const double &radius, unsigned short part)
{
    clear();
    if (part == 1)
    {
        xMin_ = center.x() - radius;
        xMax_ = center.x() + radius;
        yMin_ = center.y() - radius;
        yMax_ = center.y() + radius;

        auto topCircle = [&](const double &xi) {
            return Point2D(center.x() + radius * cos(M_PI_2 - M_PI_2 * xi),
                           center.y() + radius * sin(M_PI_2 - M_PI_2 * xi));
        };
        auto bottomCircle = [&](const double &xi) {
            return Point2D(center.x() + 0.5 * radius * xi,
                           center.y() + 0.5 * radius - 0.5 * radius * xi);
        };
        auto leftCircle = [&](const double &eta) {
            return Point2D(center.x(),
                           center.y() + 0.5 * radius + 0.5 * radius * eta);
        };
        auto rightCircle = [&](const double &eta) {
            return Point2D(center.x() + 0.5 * radius + 0.5 * radius * eta,
                           center.y());
        };
        addTransfiniteMesh(topCircle, bottomCircle, leftCircle, rightCircle, count, count);


        // Отражение по вертикали
        UInteger ec = elementsCount();
        for (UInteger i = 0; i < ec; i++)
        {
            Quadrilateral q = element_[i];
            q[0] = addNode(Point2D(xMin_ + xMax_ - node_[q[0]].point.x(), node_[q[0]].point.y()), node_[q[0]].type);
            q[1] = addNode(Point2D(xMin_ + xMax_ - node_[q[1]].point.x(), node_[q[1]].point.y()), node_[q[1]].type);
            q[2] = addNode(Point2D(xMin_ + xMax_ - node_[q[2]].point.x(), node_[q[2]].point.y()), node_[q[2]].type);
            q[3] = addNode(Point2D(xMin_ + xMax_ - node_[q[3]].point.x(), node_[q[3]].point.y()), node_[q[3]].type);
            std::swap(q[1], q[3]);
            addElement(q);
        }
        // отражение по горизонтали
        UInteger ecc = elementsCount();
        for (UInteger i = 0; i < ecc; i++)
        {
            Quadrilateral q = element_[i];
            q[0] = addNode(Point2D(node_[q[0]].point.x(), yMin_ + yMax_ - node_[q[0]].point.y()), node_[q[0]].type);
            q[1] = addNode(Point2D(node_[q[1]].point.x(), yMin_ + yMax_ - node_[q[1]].point.y()), node_[q[1]].type);
            q[2] = addNode(Point2D(node_[q[2]].point.x(), yMin_ + yMax_ - node_[q[2]].point.y()), node_[q[2]].type);
            q[3] = addNode(Point2D(node_[q[3]].point.x(), yMin_ + yMax_ - node_[q[3]].point.y()), node_[q[3]].type);
            std::swap(q[1], q[3]);
            addElement(q);
        }

        auto topQuad = [&](const double &xi) {
            return Point2D(center.x() + 0.5 * radius * xi,
                           center.y() + 0.5 * radius - 0.5 * radius * xi);
        };
        auto bottomQuad = [&](const double &xi) {
            return Point2D(center.x() - 0.5 * radius + 0.5 * radius * xi,
                           center.y() - 0.5 * radius * xi);
        };
        auto leftQuad = [&](const double &eta) {
            return Point2D(center.x() - 0.5 * radius + 0.5 * radius * eta,
                           center.y() + 0.5 * radius * eta);
        };
        auto rightQuad = [&](const double &eta) {
            return Point2D(center.x() + 0.5 * radius * eta,
                           center.y() - 0.5 * radius + 0.5 * radius * eta);
        };
        addTransfiniteMesh(topQuad, bottomQuad, leftQuad, rightQuad, count, count);
    }
    else
    {
        auto topCircle1 = [&](const double &xi) {
            return Point2D(center.x() + radius * cos(M_PI_2 - M_PI_4 * xi),
                           center.y() + radius * sin(M_PI_2 - M_PI_4 * xi));
        };
        auto bottomCircle1 = [&](const double &xi) {
            return Point2D(center.x() + 0.375 * radius * xi,
                           center.y() + 0.5 * radius - (0.5 - 0.375) * radius * xi);
        };
        auto leftCircle1 = [&](const double &eta) {
            return Point2D(center.x(),
                           center.y() + 0.5 * radius + 0.5 * radius * eta);
        };
        auto rightCircle1 = [&](const double &eta) {
            return Point2D(center.x() + 0.375 * radius + (radius * cos(M_PI_4) - 0.375 * radius) * eta,
                           center.y() + 0.375 * radius + (radius * sin(M_PI_4) - 0.375 * radius) * eta);
        };
        addTransfiniteMesh(topCircle1, bottomCircle1, leftCircle1, rightCircle1, count, count);

        auto topCircle2 = [&](const double &xi) {
            return Point2D(center.x() + radius * cos(M_PI_4 - M_PI_4 * xi),
                           center.y() + radius * sin(M_PI_4 - M_PI_4 * xi));
        };
        auto bottomCircle2 = [&](const double &xi) {
            return Point2D(center.x() + 0.375 * radius  + (0.5 - 0.375) * radius * xi,
                           center.y() + 0.3755 * radius + (0.0 - 0.375) * radius * xi);
        };
        auto leftCircle2 = [&](const double &eta) {
            return Point2D(center.x() + 0.375 * radius + (radius * cos(M_PI_4) - 0.375 * radius) * eta,
                           center.y() + 0.375 * radius + (radius * sin(M_PI_4) - 0.375 * radius) * eta);
        };
        auto rightCircle2 = [&](const double &eta) {
            return Point2D(center.x() + 0.5 * radius + 0.5 * radius * eta,
                           center.y());
        };
        addTransfiniteMesh(topCircle2, bottomCircle2, leftCircle2, rightCircle2, count, count);

        auto topQuad = [&](const double &xi) {
            return Point2D(center.x() + 0.375 * radius * xi,
                           center.y() + 0.5 * radius + (0.5 - 0.375) * radius * xi);
        };
        auto bottomQuad = [&](const double &xi) {
            return Point2D(center.x() + 0.5 * radius * xi,
                           center.y());
        };
        auto leftQuad = [&](const double &eta) {
            return Point2D(center.x(),
                           center.y() + 0.5 * radius * eta);
        };
        auto rightQuad = [&](const double &eta) {
            return Point2D(center.x() + 0.5 * radius + (0.375 - 0.5) * radius * eta,
                           center.y() + 0.375 * radius * eta);
        };
        addTransfiniteMesh(topQuad, bottomQuad, leftQuad, rightQuad, count, count);

        if (part == 4)
        {
            xMin_ = center.x();
            xMax_ = center.x() + radius;
            yMin_ = center.y();
            yMax_ = center.y() + radius;
            // Коррректировка типа для граничных узлов
            for (UInteger i = 0; i < node_.size(); i++)
            {
                if (fabs(node_[i].point.x() - center.x()) < 1.0E-7 ||
                        fabs(node_[i].point.y() - center.y()) < 1.0E-7) node_[i].type = BORDER;
            }
        }
        else
        {
            xMin_ = center.x() - radius;
            xMax_ = center.x() + radius;
            yMin_ = center.y();
            yMax_ = center.y() + radius;
            // Отражение по вертикали
            UInteger ec = elementsCount();
            for (UInteger i = 0; i < ec; i++)
            {
                Quadrilateral q = element_[i];
                q[0] = addNode(Point2D(xMin_ + xMax_ - node_[q[0]].point.x(), node_[q[0]].point.y()), node_[q[0]].type);
                q[1] = addNode(Point2D(xMin_ + xMax_ - node_[q[1]].point.x(), node_[q[1]].point.y()), node_[q[1]].type);
                q[2] = addNode(Point2D(xMin_ + xMax_ - node_[q[2]].point.x(), node_[q[2]].point.y()), node_[q[2]].type);
                q[3] = addNode(Point2D(xMin_ + xMax_ - node_[q[3]].point.x(), node_[q[3]].point.y()), node_[q[3]].type);
                std::swap(q[1], q[3]);
                addElement(q);

            }
            // Коррректировка типа для граничных узлов
            for (UInteger i = 0; i < node_.size(); i++)
            {
                if (fabs(node_[i].point.y() - center.y()) < 1.0E-7) node_[i].type = BORDER;
            }
        }
    }
    // Коррректировка типа для граничных узлов
    for (UInteger i = 0; i < node_.size(); i++)
    {
        auto circle = [&](const double &x, const double &y) { return (x - center.x()) * (x - center.x()) + (y - center.y()) * (y - center.y()) - radius * radius; };
        if (fabs(circle(node_[i].point.x(), node_[i].point.y())) < 1.0E-7) node_[i].type = BORDER;
    }
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const QuadrilateralMesh2D &mesh) :
    Mesh2D(&mesh)
{
    element_ = mesh.element_;
    node_ = mesh.node_;
    xMin_ = mesh.xMin_;
    xMax_ = mesh.xMax_;
    yMin_ = mesh.yMin_;
    yMax_ = mesh.yMax_;
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const QuadrilateralMesh2D *mesh) :
    Mesh2D(mesh)
{
    element_ = mesh->element_;
    node_ = mesh->node_;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
}

void QuadrilateralMesh2D::functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    const double iso_dist = sqrt(hx*hx + hy*hy);
    std::map<UInteger, UInteger> nodesMap;
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

            if (func(x, y) > 0.0)
            {
#ifdef WITH_OPENMP
#pragma omp critical
#endif
                nodesMap[i * yCount + j] = pushNode(point, INNER);
            }
        }
    }
    // формирование начальной сетки
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
#ifdef WITH_OPENMP
#pragma omp critical
#endif
                        addElement(iter0->second, iter1->second, iter2->second, iter3->second);
            } // if
        } // for j
    } // for i

    UInteger baseElementCount = elementsCount();
    std::vector<Point2D> normal(nodesCount()); // нормали к узлам
    // построение нормалей для всех узлов начальной сетки
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for (UInteger i = 0; i < normal.size(); i++)
    {
        Point2D currentPoint = node_[i].point;
        AdjacentSet quads = node_[i].adjacent;
        // обработка смежных элементов
        for (AdjacentSet::iterator it = quads.begin(); it != quads.end(); ++it)
        {
            Quadrilateral quad = element_[*it];
            Point2D prevPoint;
            Point2D nextPoint;
            Point2D prevNormal;
            Point2D nextNormal;
            if (quad[0] == i)
            {
                prevPoint = node_[quad[3]].point;
                nextPoint = node_[quad[1]].point;
            }
            else if (quad[1] == i)
            {
                prevPoint = node_[quad[0]].point;
                nextPoint = node_[quad[2]].point;
            }
            else if (quad[2] == i)
            {
                prevPoint = node_[quad[1]].point;
                nextPoint = node_[quad[3]].point;
            }
            else // if (quad[3] == i)
            {
                prevPoint = node_[quad[2]].point;
                nextPoint = node_[quad[0]].point;
            }
            prevNormal.set((currentPoint.y() - prevPoint.y()), -(currentPoint.x() - prevPoint.x()));
            nextNormal.set((nextPoint.y() - currentPoint.y()), -(nextPoint.x() - currentPoint.x()));
            normal[i] = normal[i] + prevNormal.normalized() + nextNormal.normalized();
        }
    }
    for (short iter = 0; iter < 2; iter++)
    {
        for (UInteger i = 0; i < normal.size(); i++)
        {
            if (normal[i].length() > epsilon_)
            {
                Point2D nn = normal[i].normalized();
                AdjacentSet adjacent = node_[i].adjacent;
                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                {
                    Quadrilateral aquad = element_[*it];
                    for (int nnode = 0; nnode < 4; nnode++)
                    {
                        if (aquad[nnode] == i)
                        {
                            UInteger prev_index = aquad[nnode - 1];
                            UInteger next_index = aquad[nnode + 1];
                            int prev_count = 0;
                            int next_count = 0;
                            for (AdjacentSet::iterator itt = adjacent.begin(); itt != adjacent.end(); ++itt)
                            {
                                Quadrilateral lquad = element_[*itt];
                                if (lquad.in(prev_index))
                                {
                                    ++prev_count;
                                }
                                if (lquad.in(next_index))
                                {
                                    ++next_count;
                                }
                            }
                            if (prev_count == 1)
                            {
                                nn = nn + normal[prev_index].normalized();
                            }
                            if (next_count == 1)
                            {
                                nn = nn + normal[next_index].normalized();
                            }
                            break;
                        }
                    }
                }
                normal[i] = nn.normalized();

            }
        } // for i
    }
    // поиск изо-точек на границе
    std::vector<UInteger> iso(nodesCount()); // изо-точки (номера)
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
            Point2D outer = current + iso_dist * n;
            Point2D borderPoint;
            if (func(outer.x(), outer.y()) > 0.0)
            {
                // внешняя точка перескачила через границу и попала внутрь
                // ищем внешнюю точку пошагово сканированием
                outer = current;
                while (func(outer.x(), outer.y()) >= 0.0)
                    outer = outer + (0.0001 * iso_dist) * n;
            }
            borderPoint = binary(inner, outer, func);
#ifdef WITH_OPENMP
#pragma omp critical
#endif
            iso[i] = pushNode(borderPoint, BORDER);
        } // if
    } // for i
    // для характерных точек, которым не нашлась пара на этапе формирования нормалей, используем метод близжайшего узла
    for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
    {
        std::vector<Node2D>::iterator min_n;
        double min_d = 10.0 * iso_dist;
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
            min_n->type = CHARACTER;
        }
    }
    // Форимрование приграничного слоя элементов
    for (UInteger i = 0; i < baseElementCount; i++)
    {
        Quadrilateral quad = element_[i];
        for (int j = 0; j < 4; j++)
        {
            if (iso[quad[j]] < ULONG_MAX && iso[quad[j + 1]] < ULONG_MAX)
            {
                AdjacentSet adjacent = node_[quad[j]].adjacent;
                int adjacnetCount = 0;
                // обработка смежных элементов
                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                {
                    Quadrilateral aquad = element_[*it];
                    if (aquad.in(quad[j]) && aquad.in(quad[j + 1]))
                        ++adjacnetCount;
                }
                if (adjacnetCount == 1)
                {
                    addElement(quad[j + 1], quad[j], iso[quad[j]], iso[quad[j + 1]]);
                }
            }
        }
    }
    // проект сглаживания путем локальной минимизации функционала
//    for (UInteger i = 0; i < nodesCount(); i++)
//    {
//        if (node_[i].type == INNER)
//        {
//            Point2D point = node_[i].point;
//            double x0[] = {point.x(), point.y()};
//            double x[2];
//            std::cout << functional(x0, i) << std::endl;
//            conjugateGradient(2, x0, x, i, 0.0001 * iso_dist);
//            node_[i].point = Point2D(x[0], x[1]);
//        }
//    } // for i
    // сглаживание Лапласа
    for (int iter = 0; iter < 5; iter++)
    {
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            if (node_[i].type == INNER)
            {
                Point2D nn = node_[i].point;
                AdjacentSet adjacent = node_[i].adjacent;
                int acount = 1;
//                double w_sum = 0.0;
                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                {
                    Quadrilateral aquad = element_[*it];
                    for (int nnode = 0; nnode < 4; nnode++)
                    {
                        if (aquad[nnode] == i)
                        {
                            UInteger prev_index = aquad[nnode - 1];
                            UInteger next_index = aquad[nnode + 1];
                            nn = nn + node_[prev_index].point;
                            ++acount;
                            nn = nn + node_[next_index].point;
                            ++acount;
                        }
                    }
                }
                node_[i].point = (1.0 / (double)acount) * nn;
            }
        } // for i
    }

    // the end.
    std::cout << "Создана сетка четырехугольных элементов для функционального объекта: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

UInteger QuadrilateralMesh2D::elementsCount() const
{
    return element_.size();
}

ElementPointer QuadrilateralMesh2D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

void QuadrilateralMesh2D::minimizeFunctional()
{
    std::vector<double> x0; // первое приближение
    for (UInteger i = 0; i < node_.size(); i++)
    {
            x0.push_back(node_[i].point.x());
            x0.push_back(node_[i].point.y());
    }
    std::cout << functional(x0) << std::endl;
    double h = 1.0E-5 * (xMax_ - xMin_);
    CoordinateFunction func = std::bind(&QuadrilateralMesh2D::functional, this, std::placeholders::_1);
    std::vector<double> x = descentGradient(func, x0, h, epsilon_);
    std::cout << functional(x) << std::endl;
    for (UInteger i = 0; i < node_.size(); i++)
    {
        if (node_[i].type == INNER)
        {
            node_[i].point.set(x[2UL * i], x[2UL * i + 1]);
        }
    }
}

Quadrilateral QuadrilateralMesh2D::quadrilateral(const UInteger &number) const
{
    return element_[number];
}

void QuadrilateralMesh2D::directionChange()
{
    for (UInteger i = 0; i < elementsCount(); i++)
    {
        std::swap(element_[i][1], element_[i][3]);
    }
}

double QuadrilateralMesh2D::area(const UInteger &number) const
{
    Quadrilateral quad = element_[number];
    Point2D p0 = node_[quad[0]].point;
    Point2D p1 = node_[quad[1]].point;
    Point2D p2 = node_[quad[2]].point;
    Point2D p3 = node_[quad[3]].point;
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

void QuadrilateralMesh2D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3)
{
    Quadrilateral quad(node0, node1, node2, node3);
    addElement(quad);
}

void QuadrilateralMesh2D::addElement(const Quadrilateral &quad)
{
    element_.push_back(quad);
    // обновление списка смежных узлов
    node_[quad[0]].adjacent.insert(element_.size() - 1);
    node_[quad[1]].adjacent.insert(element_.size() - 1);
    node_[quad[2]].adjacent.insert(element_.size() - 1);
    node_[quad[3]].adjacent.insert(element_.size() - 1);
}

void QuadrilateralMesh2D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0], nodes_ref[1], nodes_ref[2], nodes_ref[3]);
}

void QuadrilateralMesh2D::clearElements()
{
    element_.clear();
}

double QuadrilateralMesh2D::isoFunc(const UInteger &i, const double &xi, const double &eta)
{
    switch (i)
    {
    case 0:
        return (1.0 - xi) * (1.0 - eta) / 4.0;
    case 1:
        return (1.0 + xi) * (1.0 - eta) / 4.0;
    case 2:
        return (1.0 + xi) * (1.0 + eta) / 4.0;
    case 3:
        return (1.0 - xi) * (1.0 + eta) / 4.0;
    }
    return (1.0 - xi) * (1.0 - eta) / 4.0;
}

double QuadrilateralMesh2D::functional(const std::vector<double> &vars)
{
    double f = 0.0;
    const double alpha = 0.999;
    const double beta = 1.0 - alpha;
    for (UInteger i = 0; i < element_.size(); i++)
    {
        double x[4];
        double y[4];
        for (int j = 0; j < 4; j++)
        {
            if(node_[element_[i][j]].type == INNER)
            {
                x[j] = vars[ 2UL * element_[i][j] ];
                y[j] = vars[ 2UL * element_[i][j] + 1 ];
            }
            else
            {
                // узел граничный
                x[j] = node_[element_[i][j]].point.x();
                y[j] = node_[element_[i][j]].point.y();
            }
        }
        // Рассмотрим 4х-угольник как 4 треугольника, определенных на его углах
        double a [][4] = {{x[0],    x[1],   x[2],   x[3]},
                          {y[0],    y[1],   y[2],   y[3]}};
        double b [][4] = {{x[1],    x[2],   x[3],   x[0]},
                          {y[1],    y[2],   y[3],   y[0]}};
        double c [][4] = {{x[3],    x[0],   x[1],   x[2]},
                          {y[3],    y[0],   y[1],   y[2]}};
        // функция для вычисления квадрата числа (C++0x)
        auto sqr = [](double value) { return value * value; };
        for (int q = 0; q < 4; q++)
        {
//            double length = sqr(c[0][q] - a[0][q]) + sqr(c[1][q] - a[1][q]) + sqr(b[0][q] - a[0][q]) + sqr(b[1][q] - a[1][q]);
//            double ortho = (c[0][q] - a[0][q]) *  (b[0][q] - a[0][q]) + (c[1][q] - a[1][q]) * (b[1][q] - a[1][q]);
            double area = sqr((b[0][q] - a[0][q]) * (c[1][q] - a[1][q]) - (b[1][q] - a[1][q]) * (c[0][q] - a[0][q]));
//            double r = regular(a[0][q], a[1][q], 0.8, 6);
            double r = con(circle(a[0][q], a[1][q], 0.8), -circle(a[0][q], a[1][q], 0.6));
            f += area / (1.0 + (r > 0.0 ? 10.0 * r : -r));
//            f +=  alpha * length + beta *  (1.0 - 1.0 / cosh(r));
//            localValue += sqr(0.5 * alpha);
        }
    }
    return f;
}

template<typename TopFunc, typename BottomFunc, typename LeftFunc, typename RightFunc>
void QuadrilateralMesh2D::addTransfiniteMesh(TopFunc top, BottomFunc bottom, LeftFunc left, RightFunc right, const UInteger &xiCount, const UInteger &etaCount)
{
    const double hXi = 1.0 / (double)(xiCount - 1);
    const double hEta = 1.0 / (double)(etaCount - 1);
    double xi = 0.0;
    double eta = 0.0;
    UInteger nodeNumber[xiCount * etaCount];
    Point2D rb0 = bottom (0.0);
    Point2D rb1 = bottom (1.0);
    Point2D rt0 = top (0.0);
    Point2D rt1 = top (1.0);
    for (UInteger i = 0; i < xiCount; i++)
    {
        eta = 0.0;
        if (i == xiCount - 1) xi = 1.0;
        for (UInteger j = 0; j < etaCount; j++)
        {
            Point2D rt = top (xi);
            Point2D rb = bottom (xi);
            Point2D rl = left (eta);
            Point2D rr = right (eta);
            Point2D p;
            if (j == etaCount - 1) eta = 1.0;
            p = (1.0 - xi) * rl + xi * rr + (1.0 - eta) * rb + eta * rt -
                    (1.0 - xi) * (1.0 - eta) * rb0 - (1.0 - xi) * eta * rt0 -
                    xi * (1.0 - eta) * rb1 - xi * eta * rt1;
            if (i == 0 || j ==0 || i == xiCount - 1 || j == etaCount - 1)
                nodeNumber[i * etaCount + j] = addNode(p, INNER);
            else
                nodeNumber[i * etaCount + j] = pushNode(p, INNER);
            eta += hEta;
        }
        xi += hXi;
    }
    // формирование массива элементов
    for (UInteger i = 0; i < xiCount - 1; i++)
    {
        for (UInteger j = 0; j < etaCount - 1; j++)
        {
            addElement(nodeNumber[i * etaCount + j],
                    nodeNumber[(i + 1) * etaCount + j],
                    nodeNumber[(i + 1) * etaCount + j + 1],
                    nodeNumber[i * etaCount + j + 1]);
        }
    }
}

}
