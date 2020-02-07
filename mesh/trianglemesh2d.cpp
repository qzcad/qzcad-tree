#include "trianglemesh2d.h"
#include <iostream>
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#include <algorithm>
#include <float.h>
#include <ctime>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "consoleprogress.h"

#include "segment.h"
#include "funcopt.h"

namespace msh
{

TriangleMesh2D::TriangleMesh2D() : Mesh2D(nullptr)
{
}

void TriangleMesh2D::fromQuadrilateralMesh(const QuadrilateralMesh2D *mesh)
{
    for (UInteger i = 0; i < mesh->nodesCount(); i++)
    {
        pushNode(mesh->point2d(i), mesh->nodeType(i));
    }
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        Quadrilateral q = mesh->quadrilateral(i);
        Point2D p[4];
        for (UInteger j = 0; j < 4; j++)
        {
            p[j] = node_[q[j]].point;
        }
        if (std::min(minAngle(p[0], p[1], p[2]), minAngle(p[0], p[2], p[3])) > std::min(minAngle(p[0], p[1], p[3]), minAngle(p[1], p[2], p[3])))
        {
            addElement(q[0], q[1], q[2]);
            addElement(q[0], q[2], q[3]);
//            std::cout << '+';
        }
        else
        {
            addElement(q[0], q[1], q[3]);
            addElement(q[1], q[2], q[3]);
//            std::cout << '-';
        }
    }
    updateDomain();
}

void TriangleMesh2D::rectangleDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height)
{
    clear();
    double hx = width / static_cast<double>(xCount - 1);
    double hy = height / static_cast<double>(yCount - 1);
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double xCenter = (xMax_ + xMin_) / 2.0;
    const double yCenter = (yMax_ + yMin_) / 2.0;
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

void TriangleMesh2D::functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint, int smooth, int optimize)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / static_cast<double>(xCount - 1);
    const double hy = height / static_cast<double>(yCount - 1);
    const double iso_dist = sqrt(hx*hx + hy*hy);
//    double minDistance = 0.4 * sqrt(hx*hx + hy*hy);
    std::clock_t start;
    double duration;
    std::cout << "Building initial mesh..." << std::endl;
    ConsoleProgress progress_bar(xCount-1);
    start = std::clock();
    const double xCenter = (xMax_ + xMin_) / 2.0;
    const double yCenter = (yMax_ + yMin_) / 2.0;
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
                UInteger p0 = addNode(point0, INNER);
                UInteger p1 = addNode(point1, INNER);
                UInteger p2 = addNode(point2, INNER);
                UInteger p3 = addNode(point3, INNER);
                if ((xCenter <= point0.x() && yCenter <= point0.y() &&
                     xCenter <= point1.x() && yCenter <= point1.y() &&
                     xCenter <= point2.x() && yCenter <= point2.y() &&
                     xCenter <= point3.x() && yCenter <= point3.y())
                        ||
                        (xCenter >= point0.x() && yCenter >= point0.y() &&
                         xCenter >= point1.x() && yCenter >= point1.y() &&
                         xCenter >= point2.x() && yCenter >= point2.y() &&
                         xCenter >= point3.x() && yCenter >= point3.y()))
                {
#ifdef WITH_OPENMP
#pragma omp critical
                    {
#endif
                        addElement(p0, p1, p3);
                        addElement(p1, p2, p3);
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
                    addElement(p0, p1, p2);
                    addElement(p0, p2, p3);
#ifdef WITH_OPENMP
                    } // omp critical
#endif
                } // else
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
            Triangle element = element_[*it];
            Point2D prevPoint = node_[element[element.index(i) - 1]].point;
            Point2D nextPoint = node_[element[element.index(i) + 1]].point;
            Point2D prevNormal;
            Point2D nextNormal;
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
                for (UInteger elnum: adjacent)
                {
                    Triangle element = element_[elnum];
                    UInteger prev_index = element[element.index(i) - 1];
                    UInteger next_index = element[element.index(i) + 1];
                    std::vector<UInteger> prev_intersect;
                    std::vector<UInteger> next_intersect;
                    set_intersection(adjacent.begin(), adjacent.end(), node_[prev_index].adjacent.begin(), node_[prev_index].adjacent.end(), std::back_inserter(prev_intersect));
                    set_intersection(adjacent.begin(), adjacent.end(), node_[next_index].adjacent.begin(), node_[next_index].adjacent.end(), std::back_inserter(next_intersect));
                    if (prev_intersect.size() == 1)
                    {
                        nn = nn + normal[prev_index].normalized();
                    }
                    if (next_intersect.size() == 1)
                    {
                        nn = nn + normal[next_index].normalized();
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
    std::cout << "Building boundary elements..." << std::endl;
    start = std::clock();
    // Форимрование приграничного слоя элементов
    for (UInteger i = 0; i < baseElementCount; i++)
    {
        Triangle element = element_[i];
        for (int j = 0; j < 3; j++)
        {
            if (iso[element[j]] < ULONG_MAX && iso[element[j + 1]] < ULONG_MAX)
            {
                AdjacentSet adjacent0 = node_[element[j]].adjacent;
                AdjacentSet adjacent1 = node_[element[j + 1]].adjacent;
                std::vector<UInteger> intersection;
                set_intersection(adjacent0.begin(), adjacent0.end(), adjacent1.begin(), adjacent1.end(), std::back_inserter(intersection));
                if (intersection.size() == 1)
                {
                    addElement(element[j + 1], element[j], iso[element[j]]);
                    addElement(element[j + 1], iso[element[j]], iso[element[j + 1]]);
                }
            }
        }
    }
    duration = ( std::clock() - start ) /  static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Done in " << duration << " seconds." << std::endl;

    optimizeBorder(func, optimize, 0.0);

    const double t = 1.0 / ((xMax_ - xMin_) * (yMax_ - yMin_) / static_cast<double>(elementsCount()));
    for (int niter = 0; niter < optimize; niter++)
    {
        progress_bar.restart(nodesCount());
        for (UInteger i = 0 ; i < nodesCount(); i++)
        {
            if (node_[i].type == INNER)
            {
                node_[i].point = evalOptimalPosition(i, t);
            }
            ++progress_bar;
        }
        flip(false);
    }

//    laplacianSmoothing(2);
//    flip();
    for (int i = 0; i < smooth; i++)
    {
        laplacianSmoothing(1);
        flip();
    }
    evalNodalValues(func);
    // the end.
    std::cout << "Создана сетка треугольных элементов для функционального объекта: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

void TriangleMesh2D::backgroundGrid(const TriangleMesh2D *mesh, std::function<double (double, double)> primary, std::list<Point2D> charPoint, double level, int smooth, int optimize)
{
    clear();
    const double t = 1.0 / (2.0 * (mesh->xMax_ - mesh->xMin_) * (mesh->yMax_ - mesh->yMin_) / static_cast<double>(mesh->elementsCount()));
    std::cout << "Background Grid Method (Triangles) for Domains from Different Materails" << std::endl;
    std::cout << "Elements: " << mesh->elementsCount() << "; Nodes: " << mesh->nodesCount();
    std::vector<Segment> edges;
//    std::vector<Node2D> nodes = mesh->node_;
    node_ = mesh->node_;
    element_ = mesh->element_;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
    edges = evalEdges();
    std::cout << "; Edges: " << edges.size() << std::endl;
    std::cout << "t = " << t << std::endl;
    bool is_corrected = true;
    std::vector<bool> is_node_corrected(nodesCount(), false);
    while (primary != nullptr && is_corrected)
    {
        is_corrected = false;
        double min_distance = mesh->xMax() - mesh->xMin();
        Segment min_edge(0, 0);
        for (Segment edge: edges)
        {
            Point2D p0 = node_[edge[0]].point;
            Point2D p1 = node_[edge[1]].point;
            double v0 = primary(p0.x(), p0.y());
            double v1 = primary(p1.x(), p1.y());
            double l = p0.distanceTo(p1);
            if (fabs(v0 - level) >= epsilon_ && fabs(v1 - level) >= epsilon_ && signbit(v0) != signbit(v1)/* && !is_node_corrected[edge[0]] && !is_node_corrected[edge[1]]*/)
            {
                Point2D p = p0 + ((-v0) / (v1 - v0)) * (p1 - p0);
                double d = p0.distanceTo(p);
                d = (d < (l - d)) ? d : (l - d);
                if (d < min_distance)
                {
                    min_distance = d;
                    min_edge = edge;
                    is_corrected = true;
                }
            }
        }
        if (is_corrected)
        {
            Point2D p0 = node_[min_edge[0]].point;
            Point2D p1 = node_[min_edge[1]].point;
            Point2D p = binary(p0, p1, primary, level);
            std::cout << min_distance << ": ";
            UInteger optimized_num = 0;
            if (p0.distanceTo(p) < p1.distanceTo(p))
            {
                node_[min_edge[0]].point = p;
                std::cout << p0.distanceTo(p) << std::endl;
                is_node_corrected[min_edge[0]] = true;
                optimized_num = min_edge[0];
            }
            else
            {
                node_[min_edge[1]].point = p;
                std::cout << p1.distanceTo(p) << std::endl;
                is_node_corrected[min_edge[1]] = true;
                optimized_num = min_edge[1];
            }
            for (Segment edge: edges)
            {
                if (edge.in(optimized_num))
                {
                    UInteger i = (edge[0] == optimized_num) ? edge[1] : edge[0];
                    if (node_[i].type == INNER && !is_node_corrected[i])
                    {
                        node_[i].point = evalOptimalPosition(i, t);
                        std::cout << '%';
                    }
                }
            }
            flip(false);
            edges = evalEdges();
        }
    }
//    node_ = nodes;
//    element_ = mesh->element_;

    std::vector<Node2D> nodes = node_;
    std::vector<Triangle> triangles = element_;
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        Triangle triangle = mesh->triangle(i);
        Point2D p0 = nodes[triangle[0]].point;
        Point2D p1 = nodes[triangle[1]].point;
        Point2D p2 = nodes[triangle[2]].point;
        double vp0 = primary(p0.x(), p0.y());
        double vp1 = primary(p1.x(), p1.y());
        double vp2 = primary(p2.x(), p2.y());
        if (fabs(vp0 - level) < epsilon_ && fabs(vp1 - level) < epsilon_ && fabs(vp2 - level) < epsilon_)
        {
            double l0 = p0.distanceTo(p1);
            double l1 = p1.distanceTo(p2);
            double l2 = p2.distanceTo(p0);
            if (l0 < l1 && l0 < l2)
            {
                nodes[triangle[0]].point = nodes[triangle[1]].point;
            }
            else if (l1 < l0 && l1 < l2)
            {
                nodes[triangle[1]].point = nodes[triangle[2]].point;
            }
            else
            {
                nodes[triangle[2]].point = nodes[triangle[0]].point;
            }
        }
    }
    clear();
    std::vector<UInteger> index(nodes.size(), nodes.size()+1);
    for (UInteger i = 0; i < nodes.size(); i++)
    {
        Point2D p = nodes[i].point;
        double v = primary(p.x(), p.y());
        if (fabs(v - level) < epsilon_)
            index[i] = addNode(p, BORDER);
        else if (v - level >= epsilon_)
            index[i] = addNode(p, INNER);
    }
    for (Triangle triangle: triangles)
    {
//        Point2D p0 = nodes[triangle[0]].point;
//        Point2D p1 = nodes[triangle[1]].point;
//        Point2D p2 = nodes[triangle[2]].point;
//        double v0 = primary(p0.x(), p0.y());
//        double v1 = primary(p1.x(), p1.y());
//        double v2 = primary(p2.x(), p2.y());
//        NodeType t0 = OUTER;
//        NodeType t1 = OUTER;
//        NodeType t2 = OUTER;
//        if (fabs(v0 - level) < epsilon_)
//            t0 = BORDER;
//        else if (v0 - level >= epsilon_)
//            t0 = INNER;
//        if (fabs(v1 - level) < epsilon_)
//            t1 = BORDER;
//        else if (v1 - level >= epsilon_)
//            t1 = INNER;
//        if (fabs(v2 - level) < epsilon_)
//            t2 = BORDER;
//        else if (v2 - level >= epsilon_)
//            t2 = INNER;
//        if ((t0 == BORDER || t0 == INNER || t0 == CHARACTER) &&
//                (t1 == BORDER || t1 == INNER || t1 == CHARACTER) &&
//                (t2 == BORDER || t2 == INNER || t2 == CHARACTER))
        if (index[triangle[0]] < nodes.size() && index[triangle[1]] < nodes.size() && index[triangle[2]] < nodes.size() &&
                index[triangle[0]] != index[triangle[1]] && index[triangle[1]] != index[triangle[2]] && index[triangle[2]] != index[triangle[0]])
        {
//            Point2D p0 = nodes[triangle[0]].point;
//            Point2D p1 = nodes[triangle[1]].point;
//            Point2D p2 = nodes[triangle[2]].point;
//            Point2D c = (1.0 / 3.0) * (p0 + p1 + p2);
//            double v = primary(c.x(), c.y());
//            if ((v - level) >= epsilon_)
                addElement(index[triangle[0]], index[triangle[1]], index[triangle[2]]);
        }
    }
//    flip();
    if (optimize > 0)
    {
        optimizeBorder(primary, optimize, level);

        flip(false);
        edges = evalEdges();

        for (int niter = 0; niter < optimize; niter++)
        {
            for (Segment edge: edges)
            {
                if ((node_[edge[0]].type == INNER) ^ (node_[edge[1]].type == INNER))
                {
                    UInteger i = ((node_[edge[0]].type == INNER)) ? edge[0] : edge[1];
                    node_[i].point = evalOptimalPosition(i, t);
                    std::cout << '*';
                }
            }
            flip(false);
            edges = evalEdges();
            std::cout << std::endl;
        }
    }
    if (smooth > 0)
    {
        for (int i = 0; i < smooth; i++)
        {
            laplacianSmoothing();
            flip(false);
        }
    }
    std::cout << "Done. Elements: " << elementsCount() << ", nodes: " << nodesCount() << std::endl;
//    updateDomain();
}

void TriangleMesh2D::backgroundGrid(const TriangleMesh2D *mesh, std::function<double(double, double)> primary, std::function<double(double, double)> secondary, std::list<Point2D> charPoint, double level, int smooth, int optimize)
{
    clear();
    const double t = 1.0 / ((mesh->xMax_ - mesh->xMin_) * (mesh->yMax_ - mesh->yMin_) / static_cast<double>(mesh->elementsCount()));
    std::cout << "Background Grid Method (Triangles) for Domains from Different Materails" << std::endl;
    std::cout << "Elements: " << mesh->elementsCount() << "; Nodes: " << mesh->nodesCount();
    std::vector<Segment> edges;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
    node_ = mesh->node_;
    element_ = mesh->element_;
    edges = evalEdges();
    std::cout << "; Edges: " << edges.size() << std::endl;
    bool is_corrected = true;
    std::vector<bool> is_node_corrected (mesh->nodesCount(), false);
    std::cout << "boundary nodes processing" << std::endl;
    while (primary != nullptr && secondary != nullptr && is_corrected)
    {
        is_corrected = false;
        double min_distance = mesh->xMax() - mesh->xMin();
        Segment min_edge(0, 0);
        for (Segment edge: edges)
        {
            Point2D p0 = node_[edge[0]].point;
            Point2D p1 = node_[edge[1]].point;
            double v0 = secondary(p0.x(), p0.y());
            double v1 = secondary(p1.x(), p1.y());
//            double vp0 = primary(p0.x(), p0.y());
//            double vp1 = primary(p1.x(), p1.y());
            double l = p0.distanceTo(p1);
            AdjacentSet a1 = node_[edge[0]].adjacent;
            AdjacentSet a2 = node_[edge[1]].adjacent;
            std::vector<UInteger> common;
            set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), std::back_inserter(common));
//            int common = 0;
//            for (Triangle triangle: element_)
//            {
//                if (triangle.in(edge[0]) && triangle.in(edge[1])) ++common;
//            }
            if (common.size() == 1 && fabs(v0 - level) >= epsilon_ && fabs(v1 - level) >= epsilon_ && signbit(v0) != signbit(v1) && node_[edge[0]].type == BORDER && node_[edge[1]].type == BORDER /*fabs(vp0 - level) < epsilon_ && fabs(vp1 - level) < epsilon_*/)
            {
                Point2D p = p0 + ((-v0) / (v1 - v0)) * (p1 - p0);
                double d = p0.distanceTo(p);
                d = (d < (l - d)) ? d : (l - d);
                if (d < min_distance)
                {
                    min_distance = d;
                    min_edge = edge;
                    is_corrected = true;
                }
                std::cout <<'*';
            }
            else if (common.size() == 2 && node_[edge[0]].type == BORDER && node_[edge[1]].type == BORDER)
            {
//                UInteger t0 = common[0];
//                UInteger t1 = common[1];
//                UInteger index0 = edge[0];
//                UInteger index1 = edge[1];
//                node_[index0].adjacent.erase(t1);
//                node_[index1].adjacent.erase(t0);
//                element_[t][i + 1] = indexf;
//                node_[indexf].adjacent.insert(t);
//                element_[t1][subindex - 1] = index2;
//                node_[index2].adjacent.insert(t1);
                UInteger n = pushNode(0.5 * (p0 + p1), INNER);
                Triangle t0 = element_[common[0]];
                Triangle t1 = element_[common[1]];
                if (t0[t0.index(edge[0]) + 1] == edge[1])
                {
                    node_[edge[0]].adjacent.erase(common[0]);
                    element_[common[0]][t0.index(edge[0])] = n;
                    node_[n].adjacent.insert(common[0]);
                    addElement(edge[0], n, t0[t0.index(edge[0]) - 1]);
                }
                else
                {
                    node_[edge[1]].adjacent.erase(common[0]);
                    element_[common[0]][t0.index(edge[1])] = n;
                    node_[n].adjacent.insert(common[0]);
                    addElement(edge[1], n, t0[t0.index(edge[1]) - 1]);
                }
                if (t1[t1.index(edge[0]) + 1] == edge[1])
                {
                    node_[edge[0]].adjacent.erase(common[1]);
                    element_[common[1]][t1.index(edge[0])] = n;
                    node_[n].adjacent.insert(common[1]);
                    addElement(edge[0], n, t1[t1.index(edge[0]) - 1]);
                }
                else
                {
                    node_[edge[1]].adjacent.erase(common[1]);
                    element_[common[1]][t1.index(edge[1])] = n;
                    node_[n].adjacent.insert(common[1]);
                    addElement(edge[1], n, t1[t1.index(edge[1]) - 1]);
                }
                std::cout << '+';
                edges = evalEdges();
                is_corrected = true;
                continue;
            }
        }
        if (is_corrected)
        {
            Point2D p0 = node_[min_edge[0]].point;
            Point2D p1 = node_[min_edge[1]].point;
            Point2D p = binary(p0, p1, secondary, level);
            p = findBorder(p, primary, 0.1 * p0.distanceTo(p1), level);
            std::cout << min_distance << "++: ";
            UInteger optimized_num = 0;
            if (p0.distanceTo(p) < p1.distanceTo(p))
            {
                node_[min_edge[0]].point = p;
                std::cout << p0.distanceTo(p) << std::endl;
                is_node_corrected[min_edge[0]] = true;
                optimized_num = min_edge[0];
            }
            else
            {
                node_[min_edge[1]].point = p;
                std::cout << p1.distanceTo(p) << std::endl;
                is_node_corrected[min_edge[1]] = true;
                optimized_num = min_edge[1];
            }
            for (Segment edge: edges)
            {
                if (edge.in(optimized_num))
                {
                    UInteger i = (edge[0] == optimized_num) ? edge[1] : edge[0];
                    if (node_[i].type == INNER && !is_node_corrected[i])
                    {
//                        Point2D nn(0.0, 0.0);
//                        AdjacentSet adjacent = node_[i].adjacent;
//                        int acount = 0;
//                        for (UInteger elnum: adjacent)
//                        {
//                            Triangle triangle = element_[elnum];
//                            int index = triangle.index(i);
//                            nn = nn + node_[triangle[index - 1]].point;
//                            nn = nn + node_[triangle[index + 1]].point;
//                            acount += 2;
//                        }
//                        node_[i].point = (1.0 /  static_cast<double>(acount)) * nn;
                        node_[i].point = evalOptimalPosition(i, t);
                        std::cout << '*';
                    }
                }
            } // for i
            flip(false);
            edges = evalEdges();
        }
    }
    std::cout << "inner nodes processing" << std::endl;
    auto func = [&] (double x, double y) -> double {
        double p = primary(x, y);
        double s = secondary(x, y);
        return p + s - sqrt(p*p + s*s);
    };
    is_corrected = true;
    while (primary != nullptr && secondary != nullptr && is_corrected)
    {
        is_corrected = false;
        double min_distance = mesh->xMax() - mesh->xMin();
        Segment min_edge(0, 0);
        for (Segment edge: edges)
        {
            Point2D p0 = node_[edge[0]].point;
            Point2D p1 = node_[edge[1]].point;
            double v0 = secondary(p0.x(), p0.y());
            double v1 = secondary(p1.x(), p1.y());
            double l = p0.distanceTo(p1);
            if (fabs(v0 - level) >= epsilon_ && fabs(v1 - level) >= epsilon_ && signbit(v0) != signbit(v1)
                    && (node_[edge[0]].type == INNER || node_[edge[1]].type == INNER))
            {
                Point2D p = p0 + ((-v0) / (v1 - v0)) * (p1 - p0);
                double d = p0.distanceTo(p);
                d = (d < (l - d)) ? d : (l - d);
                if (d < min_distance)
                {
                    min_distance = d;
                    min_edge = edge;
                    is_corrected = true;
                }
            }
        }
        if (is_corrected)
        {
            Point2D p0 = node_[min_edge[0]].point;
            Point2D p1 = node_[min_edge[1]].point;
            Point2D p = binary(p0, p1, secondary, level);
            std::cout << min_distance << ": ";
            UInteger optimized_num = 0;
            if (p0.distanceTo(p) < p1.distanceTo(p) || (node_[min_edge[0]].type == INNER && node_[min_edge[1]].type != INNER))
            {
                node_[min_edge[0]].point = p;
                std::cout << p0.distanceTo(p) << std::endl;
                is_node_corrected[min_edge[0]] = true;
                optimized_num = min_edge[0];
            }
            else if (node_[min_edge[1]].type == INNER)
            {
                node_[min_edge[1]].point = p;
                std::cout << p1.distanceTo(p) << std::endl;
                is_node_corrected[min_edge[1]] = true;
                optimized_num = min_edge[1];
            }
            else
            {
                is_corrected = false;
            }
            for (Segment edge: edges)
            {
                if (edge.in(optimized_num))
                {
                    UInteger i = (edge[0] == optimized_num) ? edge[1] : edge[0];
                    if (node_[i].type == INNER && !is_node_corrected[i])
                    {
//                        Point2D nn(0.0, 0.0);
//                        AdjacentSet adjacent = node_[i].adjacent;
//                        int acount = 0;
//                        for (UInteger elnum: adjacent)
//                        {
//                            Triangle triangle = element_[elnum];
//                            int index = triangle.index(i);
//                            nn = nn + node_[triangle[index - 1]].point;
//                            nn = nn + node_[triangle[index + 1]].point;
//                            acount += 2;
//                        }
//                        node_[i].point = (1.0 /  static_cast<double>(acount)) * nn;
                        node_[i].point = evalOptimalPosition(i, t);
                        std::cout << '%';
                    }
                }

            } // for i
            flip(false);
            edges = evalEdges();
        }
    }
    if (optimize > 0)
    {
        SegmentMesh2D smesh;
        std::vector<UInteger> snum(nodesCount()); // изо-точки (номера)
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Point2D p = node_[i].point;
            double v = func(p.x(), p.y());
            if (fabs(v - level) < epsilon_)
                snum[i] = smesh.pushNode(node_[i].point, (node_[i].type == INNER) ? BORDER : node_[i].type);
            else
                snum[i] = ULONG_MAX;
        }

        for (Triangle element: element_)
        {
            Point2D p0 = node_[element[0]].point;
            Point2D p1 = node_[element[1]].point;
            Point2D p2 = node_[element[2]].point;
            double v0 = func(p0.x(), p0.y());
            double v1 = func(p1.x(), p1.y());
            double v2 = func(p2.x(), p2.y());
            if ((fabs(v0 - level) < epsilon_ || v0 - level >= epsilon_) && (fabs(v1 - level) < epsilon_ || v1 - level >= epsilon_) && (fabs(v2 - level) < epsilon_ || v2 - level >= epsilon_))
            for (int j = 0; j < 3; j++)
            {
                if (snum[element[j]] < ULONG_MAX && snum[element[j + 1]] < ULONG_MAX)
                {
                    AdjacentSet adjacent = node_[element[j]].adjacent;
                    int adjacentCount = 0;
                    // обработка смежных элементов
                    for (UInteger adj_num: adjacent)
                    {
                        Triangle adjacent_element = element_[adj_num];
                        Point2D pa0 = node_[adjacent_element[0]].point;
                        Point2D pa1 = node_[adjacent_element[1]].point;
                        Point2D pa2 = node_[adjacent_element[2]].point;
                        double va0 = func(pa0.x(), pa0.y());
                        double va1 = func(pa1.x(), pa1.y());
                        double va2 = func(pa2.x(), pa2.y());
                        if ((fabs(va0 - level) < epsilon_ || va0 - level >= epsilon_) && (fabs(va1 - level) < epsilon_ || va1 - level >= epsilon_) && (fabs(va2 - level) < epsilon_ || va2 - level >= epsilon_) && adjacent_element.in(element[j]) && adjacent_element.in(element[j + 1]))
                            ++adjacentCount;
                    }
                    if (adjacentCount == 1)
                    {
                        smesh.addElement(snum[element[j+1]], snum[element[j]]);
                    }
                }
            }
        }
        smesh.distlenSmoothing(func, level, optimize);
        for (UInteger i = 0; i < nodesCount(); i++)
            if (snum[i] < ULONG_MAX) node_[i].point = smesh.point2d(snum[i]);

        flip(false);
        edges = evalEdges();

        for (int niter = 0; niter < optimize; niter++)
        {
            for (Segment edge: edges)
            {
//                if ((is_node_corrected[edge[0]] ^ is_node_corrected[edge[1]]))
//                {
//                    UInteger i = (is_node_corrected[edge[0]] == true) ? edge[1] : edge[0];
//                    if (node_[i].type == INNER && !is_node_corrected[i])
//                    {
//                        node_[i].point = evalOptimalPosition(i, t);
//                        std::cout << '*';
//                    }
//                }
                if ((is_node_corrected[edge[0]] ^ is_node_corrected[edge[1]])/* || (node_[edge[0]].type == INNER ^ node_[edge[1]].type == INNER)*/)
                {
                    UInteger i = (node_[edge[0]].type == INNER) ? edge[0] : edge[1];
                    if (node_[i].type == INNER && !is_node_corrected[i])
                    {
                        node_[i].point = evalOptimalPosition(i, t);
                        std::cout << '*';
                    }
                }
            }
            flip(false);
            edges = evalEdges();
            std::cout << std::endl;
        }
    }
    if (secondary != nullptr)
    {
        std::vector<double> zones(elementsCount(), 0.0);
        for (UInteger i = 0; i < elementsCount(); i++)
        {
            Triangle t = element_[i];
            Point2D p0 = node_[t[0]].point;
            Point2D p1 = node_[t[1]].point;
            Point2D p2 = node_[t[2]].point;
            double v0 = secondary(p0.x(), p0.y());
            double v1 = secondary(p1.x(), p1.y());
            double v2 = secondary(p2.x(), p2.y());
            if ((fabs(v0 - level) < epsilon_ || v0 - level >= epsilon_) && (fabs(v1 - level) < epsilon_ || v1 - level >= epsilon_) && (fabs(v2 - level) < epsilon_ || v2 - level >= epsilon_))
            {
                zones[i] = 1.0;
            }
        }
        addDataVector("zones", zones);
    }
    std::cout << "Done. Elements: " << elementsCount() << ", nodes: " << nodesCount() << std::endl;
//    updateDomain();
}

TriangleMesh2D::TriangleMesh2D(const TriangleMesh2D &mesh) : Mesh2D(&mesh)
{
    element_ = mesh.element_;
}

TriangleMesh2D::TriangleMesh2D(const TriangleMesh2D *mesh) : Mesh2D(mesh)
{
    element_ = mesh->element_;
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

void TriangleMesh2D::addElement(const Triangle &triangle)
{
    Point2D p0 = node_[triangle[0]].point;
    Point2D p1 = node_[triangle[1]].point;
    Point2D p2 = node_[triangle[2]].point;
    if (p1 < p0 && p1 < p2)
    {
        element_.push_back(Triangle(triangle[1], triangle[2], triangle[0]));
    }
    else if (p2 < p0 && p2 < p1)
    {
        element_.push_back(Triangle(triangle[2], triangle[0], triangle[1]));
    }
    else
    {
        element_.push_back(triangle);
    }
    // обновление списка смежных узлов
    node_[triangle[0]].adjacent.insert(element_.size() - 1);
    node_[triangle[1]].adjacent.insert(element_.size() - 1);
    node_[triangle[2]].adjacent.insert(element_.size() - 1);
}

void TriangleMesh2D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2)
{
    Triangle triangle(node0, node1, node2);
    addElement(triangle);
}

void TriangleMesh2D::addElement(const std::vector<UInteger> &nodes_ref)
{
    Triangle triangle(nodes_ref[0], nodes_ref[1], nodes_ref[2]);
    addElement(triangle);
}

double TriangleMesh2D::jacobian(const UInteger &elementNum) const
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

double TriangleMesh2D::minAngle(const UInteger &elNum) const
{
    const Triangle tri = element_[elNum];
    const Point2D p0 = node_[tri[0]].point;
    const Point2D p1 = node_[tri[1]].point;
    const Point2D p2 = node_[tri[2]].point;
    return minAngle(p0, p1, p2);
}

double TriangleMesh2D::maxAngle(const UInteger &elNum) const
{
    const Triangle tri = element_[elNum];
    const Point2D A = node_[tri[0]].point;
    const Point2D B = node_[tri[1]].point;
    const Point2D C = node_[tri[2]].point;
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    angles(A, B, C, alpha, beta, gamma);
    return std::max(alpha, std::max(beta, gamma));
}

double TriangleMesh2D::angleAspect(const UInteger &elNum) const
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

void TriangleMesh2D::delaunay(const SegmentMesh2D &mesh, std::function<double (double, double)> func)
{
    clear();
    SegmentMesh2D smesh = mesh;
    std::cout << "Delaunay Meshing: " << smesh.nodesCount() << " nodes, " << smesh.elementsCount() << " elements in the initial mesh." << std::endl;
    clock_t start = clock();
    if (func == nullptr) func = std::bind(&SegmentMesh2D::cfunction, &smesh, std::placeholders::_1, std::placeholders::_2);
    Triangulation triangulation = superDelaunay(&smesh, func);
    for (UInteger i = 4; i < triangulation.nodes.size(); i++) pushNode(triangulation.nodes[i], triangulation.types[i]);
    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
         triangle != triangulation.triangles.end(); ++triangle)
    {
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
            Point2D p0 = triangulation.nodes[triangle->vertexNode(0)];
            Point2D p1 = triangulation.nodes[triangle->vertexNode(1)];
            Point2D p2 = triangulation.nodes[triangle->vertexNode(2)];
            Point2D c = (1.0 / 3.0) * (p0 + p1 + p2);
            if (func(c.x(), c.y()) > 0.0)
                addElement(triangle->vertexNode(0) - 4UL, triangle->vertexNode(1) - 4UL, triangle->vertexNode(2) - 4UL);
        }
    }
    std::cout << "Done in " <<  static_cast<double>(clock() - start) /  static_cast<double>(CLOCKS_PER_SEC) << "s: " << nodesCount() << " nodes, " << elementsCount() << " elements." << std::endl;
    evalNodalValues(func);
    xMin_ = smesh.xMin(); xMax_ = smesh.xMax();
    yMin_ = smesh.yMin(); yMax_ = smesh.yMax();
}

void TriangleMesh2D::ruppert(const SegmentMesh2D &mesh, std::function<double (double, double)> func, double alpha, double max_area)
{
    clear();
    SegmentMesh2D smesh = mesh;
    std::cout << "Ruppert Meshing: " << smesh.nodesCount() << " nodes, " << smesh.elementsCount() << " elements in the initial mesh." << std::endl;
    clock_t start = clock();
    if (func == nullptr) func = std::bind(&SegmentMesh2D::cfunction, &smesh, std::placeholders::_1, std::placeholders::_2);
    Triangulation triangulation = superDelaunay(&smesh, func);
    superRuppert(triangulation, &smesh, func, alpha);
    if (max_area > 0.0)
    {
        areaRefinement(max_area, func, triangulation);
        superRuppert(triangulation, &smesh, func);
    }
    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
             triangle != triangulation.triangles.end(); )
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        Point2D center = 1.0 / 3.0 * (A + B + C);
        double val_a = func(A.x(), A.y());
        double val_b = func(B.x(), B.y());
        double val_c = func(C.x(), C.y());
        double val_center = func(center.x(), center.y());
        if (!(val_a > -epsilon_ && val_b > -epsilon_ && val_c > -epsilon_ && val_center > -epsilon_))
        {
            triangle = triangulation.triangles.erase(triangle);
        }
        else
        {
             triangle++;
        }
    }

    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
             triangle != triangulation.triangles.end(); ++triangle)
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        addElement(addNode(A, triangulation.types[triangle->vertexNode(0)]),
                addNode(B, triangulation.types[triangle->vertexNode(1)]),
                addNode(C, triangulation.types[triangle->vertexNode(2)]));
    }
    std::cout << "Done in " <<  static_cast<double>(clock() - start) /  static_cast<double>(CLOCKS_PER_SEC) << "s: " << nodesCount() << " nodes, " << elementsCount() << " elements." << std::endl;
    evalNodalValues(func);
    xMin_ = smesh.xMin(); xMax_ = smesh.xMax();
    yMin_ = smesh.yMin(); yMax_ = smesh.yMax();
}

void TriangleMesh2D::ruppert(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func, std::list<Point2D> charPoint, bool refineArea)
{
    clear();

    SegmentMesh2D mesh;//, contour;
    mesh.MarchingQuads(xCount, yCount, xMin, yMin, width, height, func, charPoint);
//    mesh.frontGraph(xCount, yCount, xMin, yMin, width, height, func, charPoint, 3);

    Triangulation triangulation = superDelaunay(&mesh, func);

    superRuppert(triangulation, &mesh, func);
    if (refineArea)
    {
        areaRefinement(3.0 * (width /  static_cast<double>(xCount - 1) * height /  static_cast<double>(yCount - 1)), func, triangulation);
        superRuppert(triangulation, &mesh, func);
    }

    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
             triangle != triangulation.triangles.end(); )
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        Point2D center = 1.0 / 3.0 * (A + B + C);
        double val_a = func(A.x(), A.y());
        double val_b = func(B.x(), B.y());
        double val_c = func(C.x(), C.y());
        double val_center = func(center.x(), center.y());
        if (!(val_a > -epsilon_ && val_b > -epsilon_ && val_c > -epsilon_ && val_center > -epsilon_))
        {
            triangle = triangulation.triangles.erase(triangle);
        }
        else
        {
             triangle++;
        }
    }

    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
             triangle != triangulation.triangles.end(); ++triangle)
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        addElement(addNode(A, triangulation.types[triangle->vertexNode(0)]),
                addNode(B, triangulation.types[triangle->vertexNode(1)]),
                addNode(C, triangulation.types[triangle->vertexNode(2)]));
    }
    evalNodalValues(func);
    xMin_ = mesh.xMin(); xMax_ = mesh.xMax();
    yMin_ = mesh.yMin(); yMax_ = mesh.yMax();
}

void TriangleMesh2D::ruppert(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double (double, double)> func_a, std::function<double (double, double)> func_b, std::list<Point2D> charPoint, double delta, bool refineArea)
{
    clear();

//    TriangleMesh2D mesh_a, mesh_b;
//    std::list<Point2D> points_a, points_b;
//    xMin_ = xMin;
//    xMax_ = xMin + width;
//    yMin_ = yMin;
//    yMax_ = yMin + height;
//    for (std::list<Point2D>::iterator p = charPoint.begin(); p != charPoint.end(); ++p)
//    {
//        Point2D point = *p;
//        if (fabs(func_a(point.x(), point.y())) < epsilon_) points_a.push_back(point);
//        if (fabs(func_b(point.x(), point.y())) < epsilon_) points_b.push_back(point);
//    }
//    mesh_a.ruppert(xCount, yCount, xMin, yMin, width, height, func_a, points_a, true);
////    for (UInteger i = 0 ; i < mesh_a.nodesCount(); i++)
////    {
////        Point2D point = mesh_a.node_[i].point;
////        if (fabs(func_b(point.x(), point.y())) < epsilon_) points_b.push_back(point);
////    }
//    mesh_b.ruppert(xCount, yCount, xMin, yMin, width, height, func_b, points_b, true);
//    node_ = mesh_a.node_;
//    element_ = mesh_a.element_;

//    for (ElementIterator el_b = mesh_b.element_.begin(); el_b != mesh_b.element_.end(); ++el_b)
//    {
//        Triangle tri = *el_b;
//        Node2D n0 = mesh_b.node_[tri[0]];
//        Node2D n1 = mesh_b.node_[tri[1]];
//        Node2D n2 = mesh_b.node_[tri[2]];
//        tri[0] = addNode(n0.point, n0.type, 100.0 * epsilon_);
//        tri[1] = addNode(n1.point, n1.type, 100.0 * epsilon_);
//        tri[2] = addNode(n2.point, n2.type, 100.0 * epsilon_);
//        addElement(tri);
//    }

    SegmentMesh2D mesh;
    mesh.MarchingQuads(xCount, yCount, xMin, yMin, width, height, func_a, func_b, charPoint, delta);

    Triangulation triangulation = superDelaunay(&mesh, nullptr);

    superRuppert(triangulation, &mesh, nullptr);

    if (refineArea)
    {
        areaRefinement(3.0 * (width /  static_cast<double>(xCount - 1) * height /  static_cast<double>(yCount - 1)), func_a, triangulation);
        areaRefinement(3.0 * (width /  static_cast<double>(xCount - 1) * height /  static_cast<double>(yCount - 1)), func_b, triangulation);
        superRuppert(triangulation, &mesh, nullptr);
    }


//////    superRuppert(triangulation);

//    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
//             triangle != triangulation.triangles.end(); )
//    {
//        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
//        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
//        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
//        Point2D center = 1.0 / 3.0 * (A + B + C);
//        double val_a_a = func_a(A.x(), A.y());
//        double val_b_a = func_a(B.x(), B.y());
//        double val_c_a = func_a(C.x(), C.y());
//        double val_center_a = func_a(center.x(), center.y());
//        double val_a_b = func_b(A.x(), A.y());
//        double val_b_b = func_b(B.x(), B.y());
//        double val_c_b = func_b(C.x(), C.y());
//        double val_center_b = func_b(center.x(), center.y());
//        //            double val_center_b = func_b(center.x(), center.y());
//        if (!((val_a_a >= -epsilon_ || val_a_b >= -epsilon_) && (val_b_a >= -epsilon_ || val_b_b >= -epsilon_) &&
//                (val_c_a >= -epsilon_ || val_c_b >= -epsilon_) && (val_center_a >= -epsilon_ || val_center_b >= -epsilon_)))
//        {
//            triangle = triangulation.triangles.erase(triangle);
//        }
//        else
//        {
//             triangle++;
//        }
//    }

//////    superRuppert(triangulation);

    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
             triangle != triangulation.triangles.end(); ++triangle)
    {
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
            Point2D A = triangulation.nodes[triangle->vertexNode(0)];
            Point2D B = triangulation.nodes[triangle->vertexNode(1)];
            Point2D C = triangulation.nodes[triangle->vertexNode(2)];
            Point2D center = 1.0 / 3.0 * (A + B + C);
            double val_a_a = func_a(A.x(), A.y());
            double val_b_a = func_a(B.x(), B.y());
            double val_c_a = func_a(C.x(), C.y());
            double val_center_a = func_a(center.x(), center.y());
            double val_a_b = func_b(A.x(), A.y());
            double val_b_b = func_b(B.x(), B.y());
            double val_c_b = func_b(C.x(), C.y());
            double val_center_b = func_b(center.x(), center.y());
            //            double val_center_b = func_b(center.x(), center.y());
            if ((val_a_a >= -epsilon_ || val_a_b >= -epsilon_) && (val_b_a >= -epsilon_ || val_b_b >= -epsilon_) &&
                    (val_c_a >= -epsilon_ || val_c_b >= -epsilon_) && (val_center_a >= -epsilon_ || val_center_b >= -epsilon_))
            {
                addElement(addNode(A, triangulation.types[triangle->vertexNode(0)]),
                        addNode(B, triangulation.types[triangle->vertexNode(1)]),
                        addNode(C, triangulation.types[triangle->vertexNode(2)]));
            }
//            else
//            {
//                std::cout << "hmm" << std::endl;
//            }
        }
    }
    xMin_ = mesh.xMin(); xMax_ = mesh.xMax();
    yMin_ = mesh.yMin(); yMax_ = mesh.yMax();
}

Triangle TriangleMesh2D::triangle(const UInteger &number) const
{
    return element_[number];
}

void TriangleMesh2D::clearElements()
{
    element_.clear();
}

void TriangleMesh2D::flip(bool print_messages)
{
    ConsoleProgress *progress=nullptr;
    UInteger fcount = 0;
    if (print_messages)
    {
        std::cout << "Flip" << std::endl;
        progress = new ConsoleProgress(element_.size());
    }
    for (UInteger t = 0; t < element_.size(); t++)
    {
        if (progress != nullptr) progress->inc();
        Triangle triangle = element_[t];
        for (int i = 0; i < 3; i++)
        {
            UInteger index0 = triangle[i];
            UInteger index1 = triangle[i + 1];
            AdjacentSet a1 = node_[index0].adjacent;
            AdjacentSet a2 = node_[index1].adjacent;
            std::vector<UInteger> common;
            set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), std::back_inserter(common));
            if (common.size() == 2)
            {
                Point2D p0 = node_[index0].point;
                Point2D p1 = node_[index1].point;
                UInteger index2 = triangle[i - 1];
                Point2D p2 = node_[index2].point;
                UInteger indexf;
                int subindex;
                Point2D f;
                UInteger t1;
                if (common[0] != t)
                {
                    t1 = common[0];
                }
                else
                {
                    t1 = common[1];
                }

                if (element_[t1][0] != index0 && element_[t1][0] != index1)
                {
                    indexf = element_[t1][0];
                    subindex = 0;
                }
                else if (element_[t1][1] != index0 && element_[t1][1] != index1)
                {
                    indexf = element_[t1][1];
                    subindex = 1;
                }
                else
                {
                    indexf = element_[t1][2];
                    subindex = 2;
                }

                f = node_[indexf].point;
                //                    double min_c = std::min(minAngle(p0, p1, p2), minAngle(p0, f, p1));
                //                    double min_n = std::min(minAngle(p0, f, p2), minAngle(p2, f, p1));
                double tp = 0.0, tq = 0.0, xc, yc, r;
                //                    if (min_n > min_c && isCrossed(p0, p1, p2, f, tp, tq) && tp > epsilon_ && tp < 1.0 - epsilon_ && tq > epsilon_ && tq < 1.0 - epsilon_ )
                if (circumCircle(f.x(), f.y(), p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y(), xc, yc, r) &&
                        isCrossed(p0, p1, p2, f, tp, tq) && tp > epsilon_ && tp < 1.0 - epsilon_ && tq > epsilon_ && tq < 1.0 - epsilon_ )
                {
                    node_[index0].adjacent.erase(t1);
                    node_[index1].adjacent.erase(t);
                    element_[t][i + 1] = indexf;
                    node_[indexf].adjacent.insert(t);
                    element_[t1][subindex - 1] = index2;
                    node_[index2].adjacent.insert(t1);
                    fcount++;
                    break;
                }
            }
        }
    }
    if (progress != nullptr) delete progress;
    if (print_messages) std::cout << "Total flips: " << fcount << std::endl;
}

double TriangleMesh2D::minAngle(const Point2D &A, const Point2D &B, const Point2D &C)
{
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    angles(A, B, C, alpha, beta, gamma);
    return std::min(alpha, std::min(beta, gamma));
}

bool TriangleMesh2D::angles(const Point2D &A, const Point2D &B, const Point2D &C, double &alpha, double &beta, double &gamma)
{
    const double a = B.distanceTo(C); // сторона, противолежащяя вершине A (BC)
    const double b = A.distanceTo(C); // сторона, противолежащяя вершине B (AC)
    const double c = A.distanceTo(B); // сторона, противолежащяя вершине C (AB)
    if (a < epsilon_ || b < epsilon_ || c < epsilon_ || fabs(a + b - c) < epsilon_ || fabs(a + c - b) < epsilon_ || fabs(b + c - a) < epsilon_)
    {
        alpha = beta = gamma = 0.0;
        return false;
    }
    if (a > b && a > c)
    {
        // Теорема косинусов
        alpha = acos((b*b + c*c - a*a) / (2.0 * b * c)); // Угол в вершине A
        // Теорема синусов
        // обеспечение вычислительной устойчисвости: значение может на бесконечно малую отклоняться от единицы для угла 90 градусов
        double sinBeta = sin(alpha) * b / a;
        beta = (sinBeta > 1.0) ? M_PI_2 : asin(sinBeta); // Угол в вершине B
        // Теорема о сумме углов треугольника
        gamma = M_PI - (alpha + beta); // Угол в вершине C
    }
    else if (b > a && b > c)
    {
        // Теорема косинусов
        beta = acos((a*a + c*c - b*b) / (2.0 * a * c)); // Угол в вершине B
        // Теорема синусов
        // обеспечение вычислительной устойчисвости: значение может на бесконечно малую отклоняться от единицы для угла 90 градусов
        double sinAlpha = sin(beta) * a / b; // Угол в вершине A
        alpha = (sinAlpha > 1.0) ? M_PI_2 : asin(sinAlpha);
        // Теорема о сумме углов треугольника
        gamma = M_PI - (alpha + beta); // Угол в вершине C
    }
    else
    {
        // Теорема косинусов
        gamma = acos((b*b + a*a - c*c) / (2.0 * b * a)); // Угол в вершине C
        // Теорема синусов
        // обеспечение вычислительной устойчисвости: значение может на бесконечно малую отклоняться от единицы для угла 90 градусов
        double sinGamma = sin(gamma) * b / c;
        beta = (sinGamma > 1.0) ? M_PI_2 : asin(sinGamma); // Угол в вершине B
        // Теорема о сумме углов треугольника
        alpha = M_PI - (gamma + beta); // Угол в вершине A
    }

    return true;
}

TriangleMesh2D::Triangulation TriangleMesh2D::superDelaunay(SegmentMesh2D *mesh, std::function<double(double, double)> func)
{
    double super_width = mesh->xMax() - mesh->xMin();
    double super_height = mesh->yMax() - mesh->yMin();
    double h = 1.5 * sqrt(super_height*super_height + super_width*super_width);
    Point2D super0(mesh->xMin() - h, mesh->yMin() - h);
    Point2D super1(mesh->xMax() + h, mesh->yMin() - h);
    Point2D super2(mesh->xMax() + h, mesh->yMax() + h);
    Point2D super3(mesh->xMin() - h, mesh->yMax() + h);
    Triangulation triangulation;
    triangulation.nodes.push_back(super0); triangulation.types.push_back(OUTER);
    triangulation.nodes.push_back(super1); triangulation.types.push_back(OUTER);
    triangulation.nodes.push_back(super2); triangulation.types.push_back(OUTER);
    triangulation.nodes.push_back(super3); triangulation.types.push_back(OUTER);
    triangulation.nodes.push_back(mesh->point2d(0)); triangulation.types.push_back(mesh->nodeType(0));
    triangulation.triangles.push_back(Triangle(0, 1, 4));
    triangulation.triangles.push_back(Triangle(1, 2, 4));
    triangulation.triangles.push_back(Triangle(2, 3, 4));
    triangulation.triangles.push_back(Triangle(3, 0, 4));
    ConsoleProgress progress(mesh->nodesCount());
    for (UInteger i = 1; i < mesh->nodesCount(); i++)
    {
        Point2D point = mesh->point2d(i);
        NodeType type = mesh->nodeType(i);
        insertDelaunayNode(point, type, triangulation);
        ++progress;
    }
    std::cout << std::endl << "Delaunay edges refinement...";
    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
    while (triangle != triangulation.triangles.end())
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        UInteger seg_num = 0;
        if (mesh->isCrossedElement(A, B, seg_num))
        {
            Point2D R = mesh->refineMidpoint(seg_num, func);
            if (insertDelaunayNode(R, BORDER, triangulation))
                triangle = triangulation.triangles.begin();
            else
                ++triangle;
            std::cout << '+';
        }
        else if (mesh->isCrossedElement(A, C, seg_num))
        {
            Point2D R = mesh->refineMidpoint(seg_num, func);
            if (insertDelaunayNode(R, BORDER, triangulation))
                triangle = triangulation.triangles.begin();
            else
                ++triangle;
            std::cout << '+';
        }
        else if (mesh->isCrossedElement(B, C, seg_num))
        {
            Point2D R = mesh->refineMidpoint(seg_num, func);
            if (insertDelaunayNode(R, BORDER, triangulation))
                triangle = triangulation.triangles.begin();
            else
                ++triangle;
            std::cout << '+';
        }
        else
        {
            ++triangle;
        }
    }
    std::cout << std::endl;
    return triangulation;
}

void TriangleMesh2D::superRuppert(TriangleMesh2D::Triangulation &triangulation, SegmentMesh2D *mesh, std::function<double(double, double)> func, double alpha)
{
    splitSegments(triangulation);
    ConsoleProgress progress(7000000UL);
    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
    while (triangle != triangulation.triangles.end() && !progress.isExpectedCount())
    {
        ++progress;
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
            Point2D A = triangulation.nodes[triangle->vertexNode(0)];
            Point2D B = triangulation.nodes[triangle->vertexNode(1)];
            Point2D C = triangulation.nodes[triangle->vertexNode(2)];
            UInteger seg_num = 0;
            bool func_val = false;
            if (func == nullptr)
                func_val = true;
            else
                func_val = (func(A.x(), A.y()) > -epsilon_ || func(B.x(), B.y()) > -epsilon_ || func(C.x(), C.y()) > -epsilon_);

            if (mesh->isCrossedElement(A, B, seg_num))
            {
                Point2D R = mesh->refineMidpoint(seg_num, func);
                if (insertDelaunayNode(R, BORDER, triangulation))
                    triangle = triangulation.triangles.begin();
                else
                    ++triangle;
            }
            else if (mesh->isCrossedElement(A, C, seg_num))
            {
                Point2D R = mesh->refineMidpoint(seg_num, func);
                if (insertDelaunayNode(R, BORDER, triangulation))
                    triangle = triangulation.triangles.begin();
                else
                    ++triangle;
            }
            else if (mesh->isCrossedElement(B, C, seg_num))
            {
                Point2D R = mesh->refineMidpoint(seg_num, func);
                if (insertDelaunayNode(R, BORDER, triangulation))
                    triangle = triangulation.triangles.begin();
                else
                    ++triangle;
            }
            else if (minAngle(A, B, C) < alpha && func_val) // 25gr
            {
                Point2D center;

                double xc = 0.0, yc = 0.0, r = 0.0;
                circumCircle(0.0, 0.0, A.x(), A.y(), B.x(), B.y(), C.x(), C.y(), xc, yc, r);
                center.set(xc, yc);
                if (mesh->isEncroached(center, seg_num))
                {
                    Point2D R = mesh->refineMidpoint(seg_num, func);
                    if (insertDelaunayNode(R, BORDER, triangulation))
                        triangle = triangulation.triangles.begin();
                    else
                    {
                        ++triangle;
                    }
                }
                else if (insertDelaunayNode(center, INNER, triangulation))
                {
                    triangle = triangulation.triangles.begin();
                }
                else
                {
                    ++triangle;
                }
            }
            else
                ++triangle;
        }
        else
            ++triangle;
    }
    std::cout << "Ruppert's niter: " << progress.count() << std::endl;
}

void TriangleMesh2D::splitSegments(TriangleMesh2D::Triangulation &triangulation)
{
    ConsoleProgress progress(4294967290UL);
    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
    while (triangle != triangulation.triangles.end() && !progress.isExpectedCount())
    {
        bool divided = false;
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
            Point2D A = triangulation.nodes[triangle->vertexNode(0)];
            Point2D B = triangulation.nodes[triangle->vertexNode(1)];
            Point2D C = triangulation.nodes[triangle->vertexNode(2)];
            Point2D AB = 0.5 * (A + B);
            Point2D BC = 0.5 * (B + C);
            Point2D CA = 0.5 * (C + A);
            double rab = A.distanceTo(B) / 2.0;
            double rbc = B.distanceTo(C) / 2.0;
            double rca = C.distanceTo(A) / 2.0;
            for (UInteger i = 4; i < triangulation.nodes.size(); i++)
            {
                if (!triangle->in(i))
                {
                    Point2D current = triangulation.nodes[i];
                    if (fabs(AB.distanceTo(current) - rab) < -epsilon_*1000.)
                    {
                        insertDelaunayNode(AB, INNER, triangulation);
                        divided =true;
                        break;
                    }
                    if (fabs(BC.distanceTo(current) - rbc) < -epsilon_*1000.)
                    {
                        insertDelaunayNode(BC, INNER, triangulation);
                        divided =true;
                        break;
                    }
                    if (fabs(CA.distanceTo(current) - rca) < -epsilon_*1000.)
                    {
                        insertDelaunayNode(CA, INNER, triangulation);
                        divided =true;
                        break;
                    }
                }
            }
        }
        if (!divided)
        {
            ++triangle;
        }
        else
        {
            triangle = triangulation.triangles.begin();
        }
        ++progress;
    }
    std::cout << "Split segments niter: " << progress.count() << std::endl;
}

void TriangleMesh2D::areaRefinement(double max_area, std::function<double (double, double)> func, TriangleMesh2D::Triangulation &triangulation)
{
    ConsoleProgress progress(4294967290UL);
    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
    while (triangle != triangulation.triangles.end() && !progress.isExpectedCount())
    {
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        Point2D center = (1.0 / 3.0) * (A + B + C);
        if (func(center.x(), center.y()) > 0.0 && fabs(signedArea(A, B, C)) > max_area)
        {
            if (insertDelaunayNode(center, INNER, triangulation))
                triangle = triangulation.triangles.begin();
            else
                ++triangle;
        }
        else
            ++triangle;
        }
        else
            ++triangle;
        ++progress;
    }
    std::cout << "Area niter: " << progress.count() << std::endl;
}

bool TriangleMesh2D::insertDelaunayNode(const Point2D &point, const NodeType &type, Triangulation &triangulation)
{
    UInteger number = triangulation.nodes.size();
    std::vector<int> power;
    std::vector<Segment> edges;

    for (std::vector<Point2D>::iterator p = triangulation.nodes.begin(); p != triangulation.nodes.end(); ++p)
        if (p->isEqualTo(point, epsilon_)) return false;

    triangulation.nodes.push_back(point);
    triangulation.types.push_back(type);

    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin(); triangle != triangulation.triangles.end(); )
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        double xc = 0.0, yc = 0.0, r = 0.0;
        Segment e0(triangle->vertexNode(0), triangle->vertexNode(1));
        Segment e1(triangle->vertexNode(1), triangle->vertexNode(2));
        Segment e2(triangle->vertexNode(2), triangle->vertexNode(0));
        bool flags[] = {false, false, false};
        if (circumCircle(point.x(), point.y(), A.x(), A.y(), B.x(), B.y(), C.x(), C.y(), xc, yc, r))
        {
            for (UInteger j = 0; j < edges.size(); j++)
            {
                if (e0.isSame(edges[j]))
                {
                    power[j] += 1;
                    flags[0] = true;
                }
                else if (e1.isSame(edges[j]))
                {
                    power[j] += 1;
                    flags[1] = true;
                }
                else if (e2.isSame(edges[j]))
                {
                    power[j] += 1;
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
            triangle = triangulation.triangles.erase(triangle);
        }
        else
        {
            ++triangle;
        }
    }

    for (UInteger j = 0; j < edges.size(); j++)
    {
        if ((power[j] % 2) == 1)
        {
            if (signedArea(triangulation.nodes[edges[j][1]], triangulation.nodes[edges[j][0]], triangulation.nodes[number]) > 0.0)
                triangulation.triangles.push_back(Triangle(edges[j][1], edges[j][0], number));
            else
                triangulation.triangles.push_back(Triangle(edges[j][0], edges[j][1], number));
        }
    }
    return true;
}

double TriangleMesh2D::signedArea(const Point2D &A, const Point2D &B, const Point2D &C)
{
    return 0.5 * ( (B.x() - A.x()) * (C.y() - A.y()) - (C.x() - A.x()) * (B.y() - A.y()) );
}

bool TriangleMesh2D::circumCircle(const double &xp, const double &yp, const double &x1, const double &y1, const double &x2, const double &y2, const double &x3, const double &y3, double &xc, double &yc, double &r)
{
    double m1, m2, mx1, mx2, my1, my2;
    double dx, dy, rsqr, drsqr;

    /* Check for coincident points */
    if (fabs(y1 - y2) < epsilon_ && fabs(y2 - y3) < epsilon_)
        return false;
    if (fabs(y2-y1) < epsilon_)
    {
        m2 = - (x3 - x2) / (y3 - y2);
        mx2 = (x2 + x3) / 2.0;
        my2 = (y2 + y3) / 2.0;
        xc = (x2 + x1) / 2.0;
        yc = m2 * (xc - mx2) + my2;
    }
    else if (fabs(y3 - y2) < epsilon_)
    {
        m1 = - (x2 - x1) / (y2 - y1);
        mx1 = (x1 + x2) / 2.0;
        my1 = (y1 + y2) / 2.0;
        xc = (x3 + x2) / 2.0;
        yc = m1 * (xc - mx1) + my1;
    }
    else
    {
        m1 = - (x2 - x1) / (y2 - y1);
        m2 = - (x3 - x2) / (y3 - y2);
        mx1 = (x1 + x2) / 2.0;
        mx2 = (x2 + x3) / 2.0;
        my1 = (y1 + y2) / 2.0;
        my2 = (y2 + y3) / 2.0;
        xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
        yc = m1 * (xc - mx1) + my1;
    }
    dx = x2 - xc;
    dy = y2 - yc;
    rsqr = dx * dx + dy * dy;
    r = sqrt(rsqr);
    dx = xp - xc;
    dy = yp - yc;
    drsqr = dx * dx + dy * dy;
    return (drsqr <= rsqr);
}

void TriangleMesh2D::subdivide(std::list<UInteger> eNumbers, std::function<double(double, double)> func)
{
    int table [][13] = {
        {0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0
        {0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1
        {0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 2
        {0, 3, 2,  3,  1,  2, -1, -1, -1, -1, -1, -1, -1}, // 3
        {0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 4
        {0, 1, 5,  1,  2,  5, -1, -1, -1, -1, -1, -1, -1}, // 5
        {0, 1, 4,  0,  4,  2, -1, -1, -1, -1, -1, -1, -1}, // 6
        {0, 3, 5,  1,  4,  3,  2,  5,  4,  3,  4,  5, -1}  // 7
    };
    AdjacentSet refined_nodes;
    AdjacentSet refined_elements;
    Node2D local[6];
    for (UInteger elnum: eNumbers)
    {
        Triangle triangle = element_[elnum];
        for (int i = 0; i < 3; i++)
        {
            AdjacentSet a = node_[triangle[i]].adjacent;
            refined_elements.insert(a.begin(), a.end());
            refined_nodes.insert(triangle[i]);
        }
    }
    for (UInteger elnum: refined_elements)
    {
        Triangle triangle = element_[elnum];
        Node2D n0 = node_[triangle[0]];
        Node2D n1 = node_[triangle[1]];
        Node2D n2 = node_[triangle[2]];
        int code = 0;
        local[0] = n0;
        local[1] = n1;
        local[2] = n2;
        local[3].point = 0.5 * (n0.point + n1.point); local[3].type = (n0.type == INNER || n1.type == INNER) ? INNER : BORDER;
        local[4].point = 0.5 * (n1.point + n2.point); local[4].type = (n1.type == INNER || n2.type == INNER) ? INNER : BORDER;
        local[5].point = 0.5 * (n2.point + n0.point); local[5].type = (n2.type == INNER || n0.type == INNER) ? INNER : BORDER;
        if (local[3].type == BORDER && func != nullptr)
            local[3].point = findBorder(local[3].point, func, 0.5 * (n0.point.distanceTo(n1.point) + n0.point.distanceTo(n2.point)));
        if (local[4].type == BORDER && func != nullptr)
            local[4].point = findBorder(local[4].point, func, 0.5 * (n0.point.distanceTo(n1.point) + n0.point.distanceTo(n2.point)));
        if (local[5].type == BORDER && func != nullptr)
            local[5].point = findBorder(local[5].point, func, 0.5 * (n0.point.distanceTo(n1.point) + n0.point.distanceTo(n2.point)));
        if (refined_nodes.find(triangle[0]) != refined_nodes.end())
            code |= 1;
        if (refined_nodes.find(triangle[1]) != refined_nodes.end())
            code |= 2;
        if (refined_nodes.find(triangle[2]) != refined_nodes.end())
            code |= 4;
        if (code != 0)
        {
            node_[triangle[0]].adjacent.erase(elnum);
            node_[triangle[1]].adjacent.erase(elnum);
            node_[triangle[2]].adjacent.erase(elnum);
            for(int i = 0; table[code][i] != -1; i += 3)
            {
                triangle[0] = addNode(local[table[code][i]]);
                triangle[1] = addNode(local[table[code][i + 1]]);
                triangle[2] = addNode(local[table[code][i + 2]]);
                if (i == 0)
                {
                    element_[elnum] = triangle;
                    node_[triangle[0]].adjacent.insert(elnum);
                    node_[triangle[1]].adjacent.insert(elnum);
                    node_[triangle[2]].adjacent.insert(elnum);
                }
                else
                {
                    addElement(triangle);
                }
            }
        }
    }
    flip();
}

std::vector<Segment> TriangleMesh2D::evalEdges()
{
    std::vector<Segment> edges;
    for (Triangle triangle: element_)
    {
        for (int j = 0; j < 3; j++)
        {
            Segment edge(triangle[j], triangle[j + 1]);
            bool exist = false;
            for (Segment e: edges)
            {
                if (e.isSame(edge))
                {
                    exist = true;
                    break;
                }
            }
            if (!exist) edges.push_back(edge);
        }
    }
    return edges;
}

Point2D TriangleMesh2D::evalOptimalPosition(const UInteger &i, double t)
{
    Point2D point = node_[i].point;
    std::vector<double> x0(2);
    std::vector<double> x(2);
    AdjacentSet adjacent = node_[i].adjacent;
    double avr_dist = 0.0;
//    double avr_area = 0.0;
//    double min_area = DBL_MAX;
    for (UInteger elnum: adjacent)
    {
        Triangle el = element_[elnum];
        int index = el.index(i);
        Point2D prev = node_[el[index - 1]].point;
        Point2D next = node_[el[index + 1]].point;
//        double tarea = area(elnum);
        avr_dist += (point.distanceTo(prev) + point.distanceTo(next) + prev.distanceTo(next)) / 3.0;
//        avr_area += tarea;
//        if (min_area > tarea) min_area = tarea;
    }

    avr_dist /= static_cast<double>(adjacent.size());

//    t = 1.0;//log(static_cast<double>(adjacent.size() * 3)) / min_area + 0.001;

//    avr_area /= static_cast<double>(adjacent.size());
//    t = log(static_cast<double>(adjacent.size() * 3)) / (2.0 * avr_area);
//    std::cout << t << "\t";
    auto functor = [&](const std::vector<double> &vars)
    {
        double f = 0.0;
        for (UInteger elnum: adjacent)
        {
            Triangle el = element_[elnum];
            Point2D p0 = (el[0] == i) ? Point2D(vars[0], vars[1]) : node_[el[0]].point;
            Point2D p1 = (el[1] == i) ? Point2D(vars[0], vars[1]) : node_[el[1]].point;
            Point2D p2 = (el[2] == i) ? Point2D(vars[0], vars[1]) : node_[el[2]].point;
            double alpha0 = (p1 - p0).product((p2 - p0));
            double alpha1 = (p2 - p1).product((p0 - p1));
            double alpha2 = (p0 - p2).product((p1 - p2));
            double beta0 = 0, beta1 = 0, beta2 = 0;
            angles(p0, p1, p2, beta0, beta1, beta2);
            f += exp(-t * alpha0) + exp(-t * alpha1) + exp(-t * alpha2);
//            f += alpha0*alpha0;
//            if (alpha0 < 0.0) std::cout << "0: " << alpha0 << std::endl; else std::cout << "ok";
//            if (alpha1 < 0.0) std::cout << "1: " << alpha1 << std::endl; else std::cout << "ok";
//            if (alpha2 < 0.0) std::cout << "2: " << alpha2 << std::endl; else std::cout << "ok";
//            double x[3] = {0.0, 0.0, 0.0};
//            double y[3] = {0.0, 0.0, 0.0};
//            for (int j = 0; j < 3; j++)
//            {
//                UInteger nnode = el[j];
//                if(nnode == i)
//                {
//                    x[j] = vars[0];
//                    y[j] = vars[1];
//                }
//                else
//                {
//                    x[j] = node_[nnode].point.x();
//                    y[j] = node_[nnode].point.y();
//                }
//            }
//            // Рассмотрим 3х-угольник как 3 треугольника, определенных на его углах
//            double a [][3] = {{x[0],    x[1],   x[2]},
//                              {y[0],    y[1],   y[2]}};
//            double b [][3] = {{x[1],    x[2],   x[0]},
//                              {y[1],    y[2],   y[0]}};
//            double c [][3] = {{x[2],    x[0],   x[1]},
//                              {y[2],    y[0],   y[1]}};
//            for (int q = 0; q < 3; q++)
//            {
//                double alpha = (b[0][q] - a[0][q]) * (c[1][q] - a[1][q]) - (b[1][q] - a[1][q]) * (c[0][q] - a[0][q]);
////                f += exp(-t * alpha);
//                f += alpha*alpha;
//                if (alpha < 0.0) std::cout << '!';
//            }
        }
        return f;
    };

    x0[0] = node_[i].point.x();
    x0[1] = node_[i].point.y();
    x = descentGradient(functor, x0, 0.001 * avr_dist, 0.00001 * avr_dist, 500, false);
    if (!isnan(x[0]) && !isnan(x[1]))
        return Point2D(x[0], x[1]);
    return node_[i].point;
}

void TriangleMesh2D::optimizeBorder(std::function<double (double, double)> func, int iiter, double level)
{
    SegmentMesh2D smesh;
    std::vector<UInteger> snum(nodesCount()); // изо-точки (номера)
    for (UInteger i = 0; i < nodesCount(); i++)
        if (node_[i].type == BORDER || node_[i].type == CHARACTER)
            snum[i] = smesh.pushNode(node_[i].point, node_[i].type);
        else
            snum[i] = ULONG_MAX;

    for (Triangle element: element_)
    {
        for (int j = 0; j < 3; j++)
        {
            if (snum[element[j]] < ULONG_MAX && snum[element[j + 1]] < ULONG_MAX)
            {
                AdjacentSet adjacent = node_[element[j]].adjacent;
                int adjacentCount = 0;
                // обработка смежных элементов
                for (UInteger adj_num: adjacent)
                {
                    Triangle adjacent_element = element_[adj_num];
                    if (adjacent_element.in(element[j]) && adjacent_element.in(element[j + 1]))
                        ++adjacentCount;
                }
                if (adjacentCount == 1)
                {
                    smesh.addElement(snum[element[j+1]], snum[element[j]]);
                }
            }
        }
    }
    smesh.distlenSmoothing(func, level, iiter);
    for (UInteger i = 0; i < nodesCount(); i++)
        if (snum[i] < ULONG_MAX) node_[i].point = smesh.point2d(snum[i]);
}

}
