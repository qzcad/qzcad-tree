#include "quadrilateralmesh2d.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#include <ctime>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "consoleprogress.h"

#include "funcopt.h"
#include "rfunctions.h"
#include "segmentmesh2d.h"

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
    double hx = width / static_cast<double>(xCount - 1);
    double hy = height / static_cast<double>(yCount - 1);
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        double x = xMin + static_cast<double>(i) * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double y = yMin + static_cast<double>(j) * hy;
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
    double hx = 2.0 / static_cast<double>(xCount - 1);
    double hy = 2.0 / static_cast<double>(yCount - 1);
    for (UInteger i = 0; i < xCount; i++)
    {
        double xi = -1.0 + static_cast<double>(i) * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double eta = -1.0 + static_cast<double>(j) * hy;
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
    double h = 2.0 / static_cast<double>(sideCount - 1); // шаг изо-сетки
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
            double xi = -1.0 + static_cast<double>(i) * h;
            for (UInteger j = 0; j < sideCount; j++)
            {
                double eta = -1.0 + static_cast<double>(j) * h;
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
        addTransfiniteMesh(topCircle, bottomCircle, leftCircle, rightCircle, count, count, false);


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
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const QuadrilateralMesh2D *mesh) :
    Mesh2D(mesh)
{
    element_ = mesh->element_;
}

void QuadrilateralMesh2D::functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, int smooth, int optimize)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    std::clock_t start;
    double duration;
    const double hx = width / static_cast<double>(xCount - 1);
    const double hy = height / static_cast<double>(yCount - 1);
    const double iso_dist = sqrt(hx*hx + hy*hy);
    std::cout << "Building initial mesh..." << std::endl;
    ConsoleProgress progress_bar(xCount-1);
    start = std::clock();
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
    for (UInteger i = 0; i < xCount-1; i++)
    {
        double x = xMin + static_cast<double>(i) * hx;
        for (UInteger j = 0; j < yCount-1; j++)
        {
            double y = yMin + static_cast<double>(j) * hy;


            if (func(x, y) > epsilon_ && func(x + hx, y) > epsilon_ && func(x + hx, y + hy) > epsilon_ && func(x, y + hy) > epsilon_)
            {
                Point2D point0(x, y);
                Point2D point1(x + hx, y);
                Point2D point2(x + hx, y + hy);
                Point2D point3(x, y + hy);
#ifdef WITH_OPENMP
#pragma omp critical
#endif
                addElement(addNode(point0, INNER), addNode(point1, INNER), addNode(point2, INNER), addNode(point3, INNER));
            }
        }
        ++progress_bar;
    }
    duration = ( std::clock() - start ) / static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
    UInteger baseElementCount = elementsCount();
    std::vector<Point2D> normal(nodesCount()); // нормали к узлам
    // построение нормалей для всех узлов начальной сетки
    std::cout << "Calculating normals..." << std::endl;
    start = std::clock();
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
    for (short iter = 0; iter < 3; iter++)
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
    duration = ( std::clock() - start ) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
    // поиск изо-точек на границе
    std::vector<UInteger> iso(nodesCount()); // изо-точки (номера)
    std::cout << "Searching iso-points...";
    progress_bar.restart(normal.size());
    start = std::clock();
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
            if (func(outer.x(), outer.y()) > epsilon_)
            {
                // внешняя точка перескачила через границу и попала внутрь
                // ищем внешнюю точку пошагово сканированием
                outer = current;
                while (func(outer.x(), outer.y()) >= epsilon_)
                    outer = outer + (0.0001 * iso_dist) * n;
            }
            borderPoint = binary(inner, outer, func);
#ifdef WITH_OPENMP
#pragma omp critical
#endif
            iso[i] = pushNode(borderPoint, BORDER);
        } // if
        ++progress_bar;
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
    duration = ( std::clock() - start ) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
//    std::cout << "Smoothing border...";
//    for (int iii = 0; iii < 4; iii++)
//    {
//        progress_bar.restart(normal.size());
//#ifdef WITH_OPENMP
//#pragma omp parallel for schedule(static, 1)
//#endif
//        for (UInteger i = 0; i < normal.size(); i++)
//        {
//            if (normal[i].length() > epsilon_ && node_[iso[i]].type != CHARACTER)
//            {
//                AdjacentSet adjacent = node_[i].adjacent;
//                Point2D nn(0.0, 0.0);
//                int acount = 0;
//                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
//                {
//                    Quadrilateral aquad = element_[*it];
//                    for (int nnode = 0; nnode < 4; nnode++)
//                    {
//                        if (aquad[nnode] == i)
//                        {
//                            UInteger prev_index = aquad[nnode - 1];
//                            UInteger next_index = aquad[nnode + 1];
//                            int prev_count = 0;
//                            int next_count = 0;
//                            for (AdjacentSet::iterator itt = adjacent.begin(); itt != adjacent.end(); ++itt)
//                            {
//                                Quadrilateral lquad = element_[*itt];
//                                if (lquad.in(prev_index))
//                                {
//                                    ++prev_count;
//                                }
//                                if (lquad.in(next_index))
//                                {
//                                    ++next_count;
//                                }
//                            }
//                            if (prev_count == 1)
//                            {
//                                nn = nn + node_[iso[prev_index]].point;
//                                acount++;
//                            }
//                            if (next_count == 1)
//                            {
//                                nn = nn + node_[iso[next_index]].point;
//                                acount++;
//                            }
//                            break;
//                        }
//                    }
//                }
//                nn.scale(1.0 / (double)acount);
//                Point2D current = node_[i].point;
//                Point2D n = (nn - current).normalized();
//                // двоичный поиск граничной точки
//                Point2D inner = current;
//                Point2D outer = current + iso_dist * n;
//                Point2D borderPoint = current;
//                if (func(outer.x(), outer.y()) > epsilon_)
//                {
//                    // внешняя точка перескачила через границу и попала внутрь
//                    // ищем внешнюю точку пошагово сканированием
//                    outer = current;
//                    short iic = 0;
//                    while (func(outer.x(), outer.y()) >= epsilon_ && iic < 10000)
//                    {
//                        outer = outer + (0.0001 * iso_dist) * n;
//                        iic++;
//                    }
//                }
//                if (func(outer.x(), outer.y()) < epsilon_) borderPoint = binary(inner, outer, func);
//#ifdef WITH_OPENMP
//#pragma omp critical
//#endif
//                node_[iso[i]].point = borderPoint;
//            } // if
//            ++progress_bar;
//        } // for i
//    }
    std::cout << "Building boundary elements..." << std::endl;
    start = std::clock();
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
    duration = ( std::clock() - start ) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
    SegmentMesh2D smesh;
    std::vector<UInteger> snum(nodesCount()); // изо-точки (номера)
    for (UInteger i = 0; i < nodesCount(); i++)
        if (node_[i].type == BORDER || node_[i].type == CHARACTER)
            snum[i] = smesh.pushNode(node_[i].point, node_[i].type);
        else
            snum[i] = ULONG_MAX;

    for (UInteger i = 0; i < elementsCount(); i++)
    {
        Quadrilateral quad = element_[i];
        for (int j = 0; j < 4; j++)
        {
            if (snum[quad[j]] < ULONG_MAX && snum[quad[j + 1]] < ULONG_MAX)
            {
                AdjacentSet adjacent = node_[quad[j]].adjacent;
                int adjacentCount = 0;
                // обработка смежных элементов
                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                {
                    Quadrilateral aquad = element_[*it];
                    if (aquad.in(quad[j]) && aquad.in(quad[j + 1]))
                        ++adjacentCount;
                }
                if (adjacentCount == 1)
                {
                    smesh.addElement(snum[quad[j+1]], snum[quad[j]]);
                }
            }
        }
    }
    smesh.distlenSmoothing(func, 0.0, optimize);
    for (UInteger i = 0; i < nodesCount(); i++)
        if (snum[i] < ULONG_MAX) node_[i].point = smesh.point2d(snum[i]);

    // сглаживание Лапласа
    laplacianSmoothing(smooth);
    localFubctionalOptimization(optimize);

//    std::list<UInteger> ee;
//    ee.push_back(0); ee.push_back(1); ee.push_back(2);
//    subdivide(ee, func);
//    subdivide(ee, func);
//    minimizeFunctional();
    evalNodalValues(func);
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
    std::clock_t start;
    double h = 0.0;
    for (UInteger i = 0; i < node_.size(); i++)
    {
            x0.push_back(node_[i].point.x());
            x0.push_back(node_[i].point.y());
    }
    for (UInteger i = 0; i < elementsCount(); i++)
    {
        Quadrilateral q = element_[i];
        Point2D p01(node_[q[0]].point, node_[q[1]].point);
        Point2D p12(node_[q[1]].point, node_[q[2]].point);
        Point2D p23(node_[q[2]].point, node_[q[3]].point);
        Point2D p30(node_[q[3]].point, node_[q[0]].point);
        /*if (p01.length() < h)*/ h += p01.length();
        /*if (p12.length() < h)*/ h += p12.length();
        /*if (p23.length() < h)*/ h += p23.length();
        /*if (p30.length() < h)*/ h += p30.length();
    }
    h /= (4.0 *  static_cast<double>(elementsCount()));
    std::cout << "The mesh-defined functional, h = " << 0.001 * h << ", epsilon = " << 0.000001 * h << ", t = " << log10(static_cast<double>(element_.size()) * 4.0) << ", zeroth approximation: " << functional(x0) << std::endl;
    start = std::clock();
    CoordinateFunction func = std::bind(&QuadrilateralMesh2D::functional, this, std::placeholders::_1);
    std::vector<double> x = descentGradient(func, x0, 0.001 * h, 0.00000001 * h, 10000, true);
//    std::vector<double> x = conjugateGradient(func, x0, 0.001 * h, 0.000001 * h, 5000, true);
    double duration = ( std::clock() - start ) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
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

void QuadrilateralMesh2D::subdivide(std::list<UInteger> eNumbers, std::function<double(double, double)> func)
{
//    int quad_table [][2][37] =
//        {
//            /*0*/{{0, 0, 3, 3, -1},
//                  {0, 3, 3, 0, -1}
//            },
//            /*1*/{{0, 0, 1, 1, 0, 0, 3, 1, 1, 3, 3, 1, -1},
//                  {0, 1, 1, 0, 1, 3, 3, 1, 1, 3, 0, 0, -1}
//            },
//            /*2*/{{0, 0, 2, 2, 0, 3, 3, 2, 2, 3, 3, 2, -1},
//                  {0, 3, 1, 0, 3, 3, 1, 1, 1, 1, 0, 0, -1}
//            },
//            /*3*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 3, 2, 1, 3, 3, 2, 2, 3, 3, 2, 2, 2, 1, 1, 2, 1, 2, 2, 1, -1},
//                  {0, 1, 1, 0, 1, 3, 2, 1, 3, 3, 2, 2, 3, 1, 1, 2, 1, 0, 0, 1, 1, 1, 2, 2, 1, 1, 0, 0, -1}
//            },
//            /*4*/{{0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 0, 2, -1},
//                  {0, 3, 3, 2, 2, 3, 3, 2, 2, 0, 0, 2, -1}
//            },
//            /*5*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 1, 1, 2, 1, 2, 3, 1, -1},
//                  {0, 1, 1, 0, 1, 3, 2, 1, 3, 3, 2, 2, 3, 3, 2, 2, 2, 2, 0, 1, 1, 1, 2, 2, 1, 1, 0, 0, -1}
//            },
//            /*6*/{{0, 0, 1, 1, 0, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 0, 1, 2, 2, 1, 1, 2, -1},
//                  {0, 3, 2, 1, 3, 3, 2, 2, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, -1}
//            },
//            /*7*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 1, 1, 2, 2, 1, 1, 2, -1},
//                  {0, 1, 1, 0, 1, 3, 2, 1, 3, 3, 2, 2, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, -1}
//            },
//            /*8*/{{0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 3, 3, -1},
//                  {0, 2, 2, 0, 2, 3, 3, 2, 2, 3, 3, 0, -1}
//            },
//            /*9*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 3, 2, 2, 3, 3, 2, 2, 1, 1, 2, 1, 2, 3, 1, -1},
//                  {0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 0, 1, 1, 1, 2, 2, 1, 1, 0, 0, -1}
//            },
//            /*A*/{{0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 0, 1, 2, 2, 1, 1, 2, -1},
//                  {0, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, -1}
//            },
//            /*B*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 1, 1, 2, 2, 1, 1, 2, -1},
//                  {0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, -1}
//            },
//            /*C*/{{0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 3, 0, 1, 2, 2, 1, 1, 2, -1},
//                  {0, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 0, 1, 2, 0, 0, 1, 1, 1, 1, 2, 2, -1}
//            },
//            /*D*/{{0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 3, 1, 1, 2, 2, 1, 1, 2, -1},
//                  {0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 0, 1, 2, 0, 0, 1, 1, 1, 1, 2, 2, -1}
//            },
//            /*E*/{{0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 3, 3, 2, 2, 2, 2, 0, 1, 1, 1, 2, 2, -1},
//                  {0, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, 2, 1, 1, 2, 1, 0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 1, -1}
//            },
//            /*F*/{{0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, -1},
//                  {0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, -1}
//            },
//        }; // the table of indexes
    int quad_table [][2][37] =
        {
            /*0*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*1*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*2*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*3*/{{0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 1, 0, 3, 2, -1},
                  {0, 3, 2, 0, 0, 2, 2, 0, 0, 2, 3, 0, 2, 3, 3, 2, -1}
            },
            /*4*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*5*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*6*/{{0, 0, 1, 1, 0, 1, 3, 3, 3, 1, 1, 3, 3, 1, 0, 3, -1},
                  {0, 3, 2, 1, 0, 1, 1, 0, 1, 1, 2, 2, 2, 2, 3, 3, -1}
            },
            /*7*/{{0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 2, 1, 3, 3, 1, 0, 3, 3, -1},
                  {0, 3, 2, 0, 0, 2, 1, 0, 0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, -1}
            },
            /*8*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*9*/{{0, 0, 2, 3, 0, 0, 2, 2, 0, 0, 3, 2, 3, 2, 2, 3, -1},
                  {0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 0, 1, 2, 3, -1}
            },
            /*A*/{{0, 0, 3, 3, -1},
                  {0, 3, 3, 0, -1}
            },
            /*B*/{{0, 0, 1, 1, 0, 0, 2, 1, 0, 0, 3, 2, 1, 1, 2, 2, 2, 2, 3, 3, -1},
                  {0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 0, 1, 2, 0, 0, 2, 3, 0, -1}
            },
            /*C*/{{0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 1, 2, 3, -1},
                  {0, 3, 3, 1, 1, 3, 3, 1, 1, 3, 3, 0, 0, 1, 1, 0, -1}
            },
            /*D*/{{0, 0, 2, 3, 0, 0, 1, 2, 0, 0, 1, 1, 2, 1, 1, 2, 3, 2, 2, 3, -1},
                  {0, 1, 1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 1, 2, 3, 3, 0, 1, 3, 3, -1}
            },
            /*E*/{{0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 1, 2, 0, 1, 3, 3, -1},
                  {0, 3, 3, 1, 1, 3, 3, 2, 2, 3, 3, 2, 2, 1, 1, 2, 0, 1, 1, 0, -1}
            },
            /*F*/{{0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, -1},
                  {0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3, 3, 2, -1}
            },
        }; // the table of indexes
    AdjacentSet refined_nodes;
    AdjacentSet refined_quads;
    Node2D local_nodes[4][4];
    const double h = 1.0 / 3.0;
    // build a set of nodes that connect refined elements
    for (UInteger elnum: eNumbers)
    {
        Quadrilateral quad = element_[elnum];
        for (int i = 0; i < 4; i++)
        {
            AdjacentSet a = node_[quad[i]].adjacent;
            refined_quads.insert(a.begin(), a.end());
            refined_nodes.insert(quad[i]);
        }
    }
    for (UInteger elnum: refined_quads)
    {
        Quadrilateral quad = element_[elnum];
        Node2D n0 = node_[quad[0]];
        Node2D n1 = node_[quad[1]];
        Node2D n2 = node_[quad[2]];
        Node2D n3 = node_[quad[3]];
        int code = 0;
        for (int i = 0; i < 4; i++)
        {
            Node2D bottom;
            Node2D top;
            if (i == 0)
            {
                bottom.point = n0.point;
                top.point = n3.point;
            }
            else if (i == 3)
            {
                bottom.point = n1.point;
                top.point = n2.point;
            }
            else
            {
                bottom.point = n0.point + static_cast<double>(i) * h * (n1.point - n0.point);
                top.point = n3.point + static_cast<double>(i) * h * (n2.point - n3.point);
            }

            for (int j = 0; j < 4; j++)
            {
                Node2D current;
                if (i == 0 && j == 0)
                    current.type = n0.type;
                else if (i == 3 && j == 0)
                    current.type = n1.type;
                else if (i == 3 && j == 3)
                    current.type = n2.type;
                else if (i == 0 && j == 3)
                    current.type = n3.type;
                else if(i == 0 && (n3.type == BORDER || n3.type == CHARACTER) && (n0.type == BORDER || n0.type == CHARACTER))
                    current.type = BORDER;
                else if(i == 3 && (n1.type == BORDER || n1.type == CHARACTER) && (n2.type == BORDER || n2.type == CHARACTER))
                    current.type = BORDER;
                else if(j == 0 && (n0.type == BORDER || n0.type == CHARACTER) && (n1.type == BORDER || n1.type == CHARACTER))
                    current.type = BORDER;
                else if(j == 3 && (n2.type == BORDER || n2.type == CHARACTER) && (n3.type == BORDER || n3.type == CHARACTER))
                    current.type = BORDER;
                else
                    current.type = INNER;

                if (j == 0)
                    current.point = bottom.point;
                else if (j == 3)
                    current.point = top.point;
                else
                    current.point = bottom.point + static_cast<double>(j) * h * (top.point - bottom.point);

                if (current.type == BORDER && func != nullptr)
                    current.point = findBorder(current.point, func, 0.5 * (n0.point.distanceTo(n1.point) + n0.point.distanceTo(n2.point)));
//                if (current.type == BORDER) current.type = CHARACTER;
                local_nodes[i][j] = current;

            }
        }
        if (refined_nodes.find(quad[0]) != refined_nodes.end())
            code |= 1;
        if (refined_nodes.find(quad[1]) != refined_nodes.end())
            code |= 2;
        if (refined_nodes.find(quad[2]) != refined_nodes.end())
            code |= 4;
        if (refined_nodes.find(quad[3]) != refined_nodes.end())
            code |= 8;

        if (code != 0)
        {
            node_[quad[0]].adjacent.erase(elnum);
            node_[quad[1]].adjacent.erase(elnum);
            node_[quad[2]].adjacent.erase(elnum);
            node_[quad[3]].adjacent.erase(elnum);
            for(int i = 0; quad_table[code][0][i] != -1; i += 4)
            {
                quad[3] = addNode(local_nodes[quad_table[code][0][i + 0]][quad_table[code][1][i + 0]]);
                quad[2] = addNode(local_nodes[quad_table[code][0][i + 1]][quad_table[code][1][i + 1]]);
                quad[1] = addNode(local_nodes[quad_table[code][0][i + 2]][quad_table[code][1][i + 2]]);
                quad[0] = addNode(local_nodes[quad_table[code][0][i + 3]][quad_table[code][1][i + 3]]);
                if (i == 0)
                {
                    element_[elnum] = quad;
                    node_[quad[0]].adjacent.insert(elnum);
                    node_[quad[1]].adjacent.insert(elnum);
                    node_[quad[2]].adjacent.insert(elnum);
                    node_[quad[3]].adjacent.insert(elnum);
                }
                else
                {
                    addElement(quad);
                }
            }
        }
    }
}

void QuadrilateralMesh2D::localFubctionalOptimization(int maxiter, bool print_values)
{
    std::cout << "Local optimization smoothing..." << std::endl;
    std::clock_t start = std::clock();
    std::vector<double> g(node_.size() * 2UL);
    // проект сглаживания путем локальной минимизации функционала
    for (int iii = 0; iii < maxiter; iii++)
    {
        if (print_values)
        {
            for (UInteger i = 0; i < node_.size(); i++)
            {
                    g[i * 2UL] = node_[i].point.x();
                    g[i * 2UL + 1UL] = node_[i].point.y();
            }
            std::cout << "F = " << functional(g) << std::endl;
        }

        ConsoleProgress progress_bar(nodesCount());

        for (UInteger i = 0; i < nodesCount(); i++)
        {
            if (node_[i].type == INNER)
            {
                Point2D point = node_[i].point;
                std::vector<double> x0(2);
                std::vector<double> x(2);
                AdjacentSet adjacent = node_[i].adjacent;
                double avr_dist = 0.0;

                for (UInteger elnum: adjacent)
                {
                    Quadrilateral el = element_[elnum];
                    int index = el.index(i);
                    Point2D prev = node_[el[index - 1]].point;
                    Point2D next = node_[el[index + 1]].point;
                    avr_dist += point.distanceTo(prev) + point.distanceTo(next);
                }

                avr_dist /= static_cast<double>(2 * adjacent.size());

                auto functor = [&](const std::vector<double> &vars){
                    double f = 0.0;
                    for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                    {
                        Quadrilateral el = element_[*it];
                        double x[4];
                        double y[4];
                        int r = 0;
                        for (int j = 0; j < 4; j++)
                        {
                            UInteger nnode = el[j];
                            if(nnode == i)
                            {
                                x[j] = vars[0];
                                y[j] = vars[1];
                                r = j;
                            }
                            else
                            {
                                x[j] = node_[nnode].point.x();
                                y[j] = node_[nnode].point.y();
                            }
                        }
                        // Рассмотрим 4х-угольник как 4 треугольника, определенных на его углах
                        double a [][4] = {{x[0],    x[1],   x[2],   x[3]},
                                          {y[0],    y[1],   y[2],   y[3]}};
                        double b [][4] = {{x[1],    x[2],   x[3],   x[0]},
                                          {y[1],    y[2],   y[3],   y[0]}};
                        double c [][4] = {{x[3],    x[0],   x[1],   x[2]},
                                          {y[3],    y[0],   y[1],   y[2]}};
                        for (int q = 0; q < 4; q++)
                        {
                            double area = (b[0][q] - a[0][q]) * (c[1][q] - a[1][q]) - (b[1][q] - a[1][q]) * (c[0][q] - a[0][q]);
//                            double area = (b[0][r] - a[0][r]) * (c[1][r] - a[1][r]) - (b[1][r] - a[1][r]) * (c[0][r] - a[0][r]);
                            f += exp(-10.0 * area);
                        }
                    }
                    return f;
                };

                x0[0] = point.x();
                x0[1] = point.y();
                x = conjugateGradient(functor, x0, 0.01 * avr_dist, 0.0001 * avr_dist, 40, false);
                node_[i].point.set(x[0], x[1]);
            }
            ++progress_bar;
        } // for i
    }
    if (print_values)
    {
        for (UInteger i = 0; i < node_.size(); i++)
        {
            g[i * 2UL] = node_[i].point.x();
            g[i * 2UL + 1UL] = node_[i].point.y();
        }
        std::cout << "F = " << functional(g) << std::endl;
    }
    double duration = static_cast<double>(std::clock() - start) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;
}

double QuadrilateralMesh2D::minAngle(const UInteger &elNum) const
{
    Quadrilateral quad = element_[elNum];
    Point2D A = node_[quad[0]].point;
    Point2D B = node_[quad[1]].point;
    Point2D C = node_[quad[2]].point;
    Point2D D = node_[quad[3]].point;
    double a = A.angle(D, B);
    double b = B.angle(A, C);
    double c = C.angle(B, D);
    double d = D.angle(C, A);
    return std::min(std::min(a, b), std::min(c, d));
}

double QuadrilateralMesh2D::maxAngle(const UInteger &elNum) const
{
    Quadrilateral quad = element_[elNum];
    Point2D A = node_[quad[0]].point;
    Point2D B = node_[quad[1]].point;
    Point2D C = node_[quad[2]].point;
    Point2D D = node_[quad[3]].point;
    double a = A.angle(D, B);
    double b = B.angle(A, C);
    double c = C.angle(B, D);
    double d = D.angle(C, A);
    return std::max(std::max(a, b), std::max(c, d));
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
//    double t = log(4.0 * element_.size() + 1.0) + 1.0;
//    double t = log10(static_cast<double>(element_.size()) * 4.0);
    double t = 10.0;
//    const double alpha = 0.5;
//    const double beta = 1.0 - alpha;
    for (auto el: element_)
    {
        double x[4];
        double y[4];
        for (int j = 0; j < 4; j++)
        {
            UInteger nnode = el[j];
            if(node_[nnode].type == INNER)
            {
                x[j] = vars[ 2UL * nnode ];
                y[j] = vars[ 2UL * nnode + 1 ];
            }
            else
            {
                // узел граничный
                x[j] = node_[nnode].point.x();
                y[j] = node_[nnode].point.y();
            }
        }
        // Рассмотрим 4х-угольник как 4 треугольника, определенных на его углах
        double a [][4] = {{x[0],    x[1],   x[2],   x[3]},
                          {y[0],    y[1],   y[2],   y[3]}};
        double b [][4] = {{x[1],    x[2],   x[3],   x[0]},
                          {y[1],    y[2],   y[3],   y[0]}};
        double c [][4] = {{x[3],    x[0],   x[1],   x[2]},
                          {y[3],    y[0],   y[1],   y[2]}};

        for (int q = 0; q < 4; q++)
        {
//            double length = (c[0][q] - a[0][q])*(c[0][q] - a[0][q]) + (c[1][q] - a[1][q])*(c[1][q] - a[1][q]) +
//                    (b[0][q] - a[0][q])*(b[0][q] - a[0][q]) + (b[1][q] - a[1][q])*(b[1][q] - a[1][q]);
//            double ortho = ((c[0][q] - a[0][q]) *  (b[0][q] - a[0][q]) + (c[1][q] - a[1][q]) * (b[1][q] - a[1][q]));
            double area = (b[0][q] - a[0][q]) * (c[1][q] - a[1][q]) - (b[1][q] - a[1][q]) * (c[0][q] - a[0][q]);
//            double r = regular(a[0][q], a[1][q], 0.8, 6);
//            double r = con(circle(a[0][q], a[1][q], 0.8), -circle(a[0][q], a[1][q], 0.6));
//            f += alpha * area * area + beta * ortho * ortho ;// / (1.0 + (r > 0.0 ? 10.0 * r : -r));
//            f += alpha * length + beta * ortho * ortho;
//            f +=  alpha * length + beta *  (1.0 - 1.0 / cosh(r));
//            localValue += sqr(0.5 * alpha);
//            f += alpha * area*area + beta * length;
//            if (area <= 0.0)
//            {
//                area *= 10.;
//                ortho *= 10.;
//                length *= 10.;
//            }
//            area *= 100.; ortho *= 100.; // при малых размерах плохая сходимость поиска
            f += exp(-t*area);//*area + ortho*ortho + length;
        }
    }
    return f;
}

template<typename TopFunc, typename BottomFunc, typename LeftFunc, typename RightFunc>
void QuadrilateralMesh2D::addTransfiniteMesh(TopFunc top, BottomFunc bottom, LeftFunc left, RightFunc right, const UInteger &xiCount, const UInteger &etaCount, bool is_uniform)
{
    const double hXi = 1.0 /  static_cast<double>(xiCount - 1);
    const double hEta = 1.0 /  static_cast<double>(etaCount - 1);
    double xi = 0.0;
    double eta = 0.0;
    UInteger nodeNumber[xiCount * etaCount];
    Point2D rb0 = bottom(0.0);
    Point2D rb1 = bottom(1.0);
    Point2D rt0 = top(0.0);
    Point2D rt1 = top(1.0);
    double pp = 1.9;
    double q = 2.0;
    for (UInteger i = 0; i < xiCount; i++)
    {
        eta = 0.0;
        if (i == xiCount - 1) xi = 1.0;
        Point2D rt = top (xi);
        Point2D rb = bottom (xi);
        for (UInteger j = 0; j < etaCount; j++)
        {
            if (j == etaCount - 1) eta = 1.0;
//            double s = pow(eta, 0.25);
            double s = (is_uniform) ? eta : pp * eta + (1.0 - pp) * (1.0 - tanh(q * (1.0 - eta)) / tanh(q)); // формула растяжения Эйземана, Флетчер, Вычислительные методы в динамике жидкостей, том 2, стр. 123
            Point2D rl = left(s);
            Point2D rr = right(s);
            Point2D p;
            p = (1.0 - xi) * rl + xi * rr + (1.0 - s) * rb + s * rt - (1.0 - xi) * (1.0 - s) * rb0 - (1.0 - xi) * s * rt0 - xi * (1.0 - s) * rb1 - xi * s * rt1;
            if (i == 0 || j == 0 || i == xiCount - 1 || j == etaCount - 1)
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
            addElement(nodeNumber[i * etaCount + j], nodeNumber[(i + 1) * etaCount + j], nodeNumber[(i + 1) * etaCount + j + 1], nodeNumber[i * etaCount + j + 1]);
        }
    }
}

}
