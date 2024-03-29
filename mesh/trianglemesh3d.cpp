#include "trianglemesh3d.h"
#define _USE_MATH_DEFINES
#include <cmath>
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <set>
#include <climits>
#include <iostream>
#include <algorithm>
#include <limits>
#include "funcopt.h"
#include "rfunctions.h"
#include "consoleprogress.h"

#include "trianglemesh2d.h"

namespace msh {

TriangleMesh3D::TriangleMesh3D() : Mesh3D(nullptr)
{
    xMin_ = -1.0;
    xMax_ = 1.0;
    yMin_ = -1.0;
    yMax_ = 1.0;
    zMin_ = -1.0;
    zMax_ = 1.0;
}

TriangleMesh3D::TriangleMesh3D(const TriangleMesh3D &mesh) : Mesh3D(&mesh)
{
    element_ = mesh.element_;
}

TriangleMesh3D::TriangleMesh3D(const TriangleMesh3D *mesh) : Mesh3D(mesh)
{
    element_ = mesh->element_;
}

void TriangleMesh3D::cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length)
{
    clear();
    double hphi = 2.0 * M_PI / static_cast<double>(rCount);
    double hl = length / static_cast<double>(lCount);
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
            UInteger p0, p1, p2, p3;
            if (i < rCount - 1)
            {
                p0 = i * (lCount + 1) + j;
                p1 = i * (lCount + 1) + j + 1;
                p2 = (i + 1) * (lCount + 1) + j + 1;
                p3 = (i + 1) * (lCount + 1) + j;
            }
            else
            {
                p0 = i * (lCount + 1) + j;
                p1 = i * (lCount + 1) + j + 1;
                p2 = j + 1;
                p3 = j;
            }
            addElement(p0, p2, p3);
            addElement(p0, p1, p2);
        }
    }
    // размеры области
    xMin_ = zMin_ = -radius;
    xMax_ = zMax_ = radius;
    yMin_ = 0.0;
    yMax_ = length;
    printStats();
}

void TriangleMesh3D::cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, std::function<double (double, double, double)> func)
{
    clear();
    const double xi_max = static_cast<double>(lCount) * 2.0 * M_PI / static_cast<double>(rCount);
    const double xi_max_2 = xi_max / 2.0;

    auto local_func = [&](double xi, double eta)
    {
        double x = radius * cos(eta);
        double y = length * xi / xi_max;
        double z = radius * sin(eta);
        return func(x, y, z);
    };

    auto func2d = [&](double xi, double eta)
    {
        double f = local_func(xi, eta);
        double r = con(xi_max_2 * xi_max_2 - (xi - xi_max_2) * (xi - xi_max_2), M_PI * M_PI - (eta - M_PI) * (eta - M_PI));
        return con(f, r);
    };

    std::list<Point2D> charPoints2d;
    TriangleMesh2D mesh2d;
    std::map<UInteger, UInteger> nodes_map;

    if (func2d(0.0, 0.0) < epsilon_ )
    {
        charPoints2d.push_back(Point2D(0.0, 0.0));
        charPoints2d.push_back(Point2D(0.0, 2.0 * M_PI));
    }
    if (func2d(xi_max, 0.0) < epsilon_ )
    {
        charPoints2d.push_back(Point2D(xi_max, 0.0));
        charPoints2d.push_back(Point2D(xi_max, 2.0 * M_PI));
    }
    double dl = xi_max / static_cast<double>(lCount);
    double dphi = 2.0 * M_PI / static_cast<double>(rCount);
    // двоичный поиск вдоль шва
    for (UInteger i = 0; i < lCount; i++)
    {
        double xi0 = static_cast<double>(i) * dl;
        double xi1 = static_cast<double>(i + 1) * dl;
        double eta = 0.0;
        double val0 = local_func(xi0, eta);
        double val1 = local_func(xi1, eta);
        if ((val0 > epsilon_ && val1 < epsilon_) || (val0 < epsilon_ && val1 > epsilon_))
        {
            Point2D p0(xi0, 0.0);
            Point2D p1(xi1, 0.0);
            Point2D b = Mesh2D::binary(p0, p1, local_func);
            charPoints2d.push_back(Point2D(b.x(), 0.0));
            charPoints2d.push_back(Point2D(b.y(), 2.0 * M_PI));
        }
        if (fabs(val0) < epsilon_ && i > 0)
        {
            charPoints2d.push_back(Point2D(xi0, 0.0));
            charPoints2d.push_back(Point2D(xi0, 2.0 * M_PI));
        }
        if (fabs(val1) < epsilon_ && i < lCount - 1)
        {
            charPoints2d.push_back(Point2D(xi1, 0.0));
            charPoints2d.push_back(Point2D(xi1, 2.0 * M_PI));
        }
    }
    // двоичный поиск вдоль граней
    for (UInteger i = 0; i < rCount; i++)
    {
        double xi0 = 0.0;
        double xi1 = xi_max;
        double eta0 = static_cast<double>(i) * dphi;
        double eta1 = static_cast<double>(i + 1) * dphi;
        for (double xi = xi0; xi <= xi1; xi += xi1)
        {
            double val0 = local_func(xi, eta0);
            double val1 = local_func(xi, eta1);
            if ((val0 > epsilon_ && val1 < epsilon_) || (val0 < epsilon_ && val1 > epsilon_))
            {
                Point2D p0(xi, eta0);
                Point2D p1(xi, eta1);
                Point2D b = Mesh2D::binary(p0, p1, local_func);
                charPoints2d.push_back(Point2D(xi, b.y()));
            }
            if (fabs(val0) < epsilon_ && i > 0)
            {
                charPoints2d.push_back(Point2D(xi, eta0));
            }
            if (fabs(val1) < epsilon_ && i < rCount - 1)
            {
                charPoints2d.push_back(Point2D(xi, eta1));
            }
        }
    }
    mesh2d.ruppert(lCount, rCount, -0.001, -0.001, xi_max + 0.002, 2.0 * M_PI + 0.002, func2d, charPoints2d, true);

    // сглаживание Лапласа, взвешенное по растояниям до соседних узлов
    for (int iter = 0; iter < 5; iter++)
    {
        for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
        {
            if (mesh2d.nodeType(i) == INNER)
            {
                Point2D s(0.0, 0.0);
                Point2D c = mesh2d.point2d(i);
                Point3D c3 (radius * cos(c.y()), c.x() * length / xi_max, radius * sin(c.y()));
                AdjacentSet adjacent = mesh2d.adjacent(i);
                double w_sum = 0.0;
                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                {
                    Triangle t = mesh2d.triangle(*it);
                    for (int nnode = 0; nnode < 3; nnode++)
                    {
                        if (t[nnode] != i)
                        {
                            Point2D p = mesh2d.point2d(t[nnode]);
                            Point3D p3(radius * cos(p.y()), p.x() * length / xi_max, radius * sin(p.y()));
                            double w = c3.distanceTo(p3);
                            w_sum += w;
                            s = s + (w * p);
                        }
                    }
                }
                mesh2d.setPoint(i, (1.0 / w_sum) * s);
            }
        } // for i
    }
    // сглаживание Лапласа, взвешенное по мерам углов
    //    for (int iter = 0; iter < 10; iter++)
    //    {
    //        for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
    //        {
    //            if (mesh2d.nodeType(i) == INNER)
    //            {
    //                Point2D s(0.0, 0.0);
    //                Point2D c = mesh2d.point2d(i);
    //                Point3D c3 (radius * cos(c.y()), c.x(), radius * sin(c.y()));
    //                AdjacentSet adjacent = mesh2d.adjacent(i);
    //                double w_sum = 0.0;
    //                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
    //                {
    //                    Triangle t = mesh2d.triangle(*it);
    //                    for (int nnode = 0; nnode < 3; nnode++)
    //                    {
    //                        if (t[nnode] == i)
    //                        {
    //                            Point2D a = mesh2d.point2d(t[nnode + 1]);
    //                            Point2D b = mesh2d.point2d(t[nnode - 1]);
    //                            Point3D a3(radius * cos(a.y()), a.x(), radius * sin(a.y()));
    //                            Point3D b3(radius * cos(b.y()), b.x(), radius * sin(b.y()));
    //                            double wa = Point3D(a3, c3).normalized() * Point3D(a3, b3).normalized();
    //                            double wb = Point3D(b3, a3).normalized() * Point3D(b3, c3).normalized();
    //                            double wc = Point3D(c3, a3).normalized() * Point3D(c3, b3).normalized();
    //                            w_sum += fabs(wa) + fabs(wb) + 2. * fabs(wc);
    //                            s = s + ((fabs(wa) + fabs(wc)) * a + (fabs(wb) + fabs(wc)) * b);
    //                        }
    //                    }
    //                }
    //                mesh2d.setPoint(i, (1.0 / w_sum) * s);
    //            }
    //        } // for i
    //    }

    //    for (int iter = 0; iter < 5; iter++)
    //    {
    //        for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
    //        {
    //            if (mesh2d.nodeType(i) == INNER)
    //            {
    ////                Point2D s(0.0, 0.0);
    //                Point2D current = mesh2d.point2d(i);
    //                std::vector<double> x0(2);
    //                std::vector<double> result;
    ////                Point3D c3 (radius * cos(c.y()), c.x(), radius * sin(c.y()));
    //                AdjacentSet adjacent = mesh2d.adjacent(i);
    //                auto functor = [&](const std::vector<double> &xxxx)
    //                {
    //                    double f = 0.0;
    //                    Point2D c(xxxx[0], xxxx[1]);
    //                    Point3D c3 (radius * cos(c.y()), c.x() * length / xi_max, radius * sin(c.y()));
    //                    for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
    //                    {
    //                        Triangle t = mesh2d.triangle(*it);
    //                        for (int nnode = 0; nnode < 3; nnode++)
    //                        {
    //                            if (t[nnode] == i)
    //                            {
    //                                Point2D a = mesh2d.point2d(t[nnode - 1]);
    //                                Point2D b = mesh2d.point2d(t[nnode + 1]);
    //                                Point3D a3(radius * cos(a.y()), a.x() * length / xi_max, radius * sin(a.y()));
    //                                Point3D b3(radius * cos(b.y()), b.x() * length / xi_max, radius * sin(b.y()));
    ////                                double alpha = 0.0, beta = 0.0, gamma = 0.0;
    //                                double area = Point3D(c3, a3).product(Point3D(c3, b3)).length();
    ////                                angles(a3, b3, c3, alpha, beta, gamma);
    //                                f +=  area * area;
    //                                break;
    //                            }
    //                        }
    //                    }
    //                    return f;
    //                };
    //                x0[0] = current.x();
    //                x0[1] = current.y();
    ////                std::cout << "f = " << functor(x0) << std::endl;
    //                result = descentGradient(functor, x0, epsilon_ * 1000.0, epsilon_ * 100., 1000, false);
    ////                std::cout << "f = " << functor(result) << std::endl;
    //                mesh2d.setPoint(i, Point2D(result[0], result[1]));
    ////                break;
    //            }
    //        } // for i
    //    }

    for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
    {
        Point2D p = mesh2d.point2d(i);
        NodeType nodeType = (fabs(func2d(p.x(), p.y())) < epsilon_ || mesh2d.nodeType(i) == CHARACTER) ? CHARACTER : BORDER;
        if (fabs(p.y() - 2.0 * M_PI) < epsilon_)
            nodes_map[i] = addNode(Point3D(radius * cos(p.y()), p.x() * length / xi_max, radius * sin(p.y())), CHARACTER);
        else
            nodes_map[i] = pushNode( Point3D(radius * cos(p.y()), p.x() * length / xi_max, radius * sin(p.y())), nodeType);
    }
    for (UInteger i = 0; i < mesh2d.elementsCount(); i++)
    {
        Triangle t = mesh2d.triangle(i);
        addElement(nodes_map[t[0]], nodes_map[t[1]], nodes_map[t[2]]);
        //        addElement(t[1], t[0], t[2]);
    }

    // размеры области
    xMin_ = zMin_ = -radius;
    xMax_ = zMax_ = radius;
    yMin_ = 0.0;
    yMax_ = length;
    printStats();
}

void TriangleMesh3D::coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length)
{
    clear();
    double hphi = 2.0 * M_PI / static_cast<double>(rCount);
    double hl = length / static_cast<double>(lCount);
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
            UInteger p0, p1, p2, p3;
            if (i < rCount - 1)
            {
                p0 = i * (lCount + 1) + j;
                p1 = i * (lCount + 1) + j + 1;
                p2 = (i + 1) * (lCount + 1) + j + 1;
                p3 = (i + 1) * (lCount + 1) + j;
            }
            else
            {
                p0 = i * (lCount + 1) + j;
                p1 = i * (lCount + 1) + j + 1;
                p2 = j + 1;
                p3 = j;
            }
            addElement(p0, p2, p3);
            addElement(p0, p1, p2);
        }
    }
    // размеры области
    xMin_ = zMin_ = -(bottom_radius > top_radius ? bottom_radius : top_radius);
    xMax_ = zMax_ = (bottom_radius > top_radius ? bottom_radius : top_radius);
    yMin_ = 0.0;
    yMax_ = length;
    printStats();
}

void TriangleMesh3D::coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, std::function<double (double, double, double)> func)
{
    clear();
    const double xi_max = static_cast<double>(lCount) * 2.0 * M_PI / static_cast<double>(rCount);
    const double xi_max_2 = xi_max / 2.0;

    auto conePoint3d = [&](const double &xi, const double &eta)
    {
        double radius = bottom_radius + xi / xi_max * (top_radius - bottom_radius);
        double x = radius * cos(eta);
        double y = length * xi / xi_max;
        double z = radius * sin(eta);
        return Point3D(x, y, z);
    };

    auto local_func = [&](double xi, double eta)
    {
        double radius = bottom_radius + xi / xi_max * (top_radius - bottom_radius);
        double x = radius * cos(eta);
        double y = length * xi / xi_max;
        double z = radius * sin(eta);
        return func(x, y, z);
    };

    auto func2d = [&](double xi, double eta)
    {
        double f = local_func(xi, eta);
        double r = con(xi_max_2 * xi_max_2 - (xi - xi_max_2) * (xi - xi_max_2), M_PI * M_PI - (eta - M_PI) * (eta - M_PI));
        return con(f, r);
    };

    std::list<Point2D> charPoints2d;
    TriangleMesh2D mesh2d;
    std::map<UInteger, UInteger> nodes_map;

    if (local_func(0.0, 0.0) > -epsilon_ )
    {
        charPoints2d.push_back(Point2D(0.0, 0.0));
        charPoints2d.push_back(Point2D(0.0, 2.0 * M_PI));
    }
    if (local_func(xi_max, 0.0) > -epsilon_ )
    {
        charPoints2d.push_back(Point2D(xi_max, 0.0));
        charPoints2d.push_back(Point2D(xi_max, 2.0 * M_PI));
    }
    double dxi = xi_max / static_cast<double>(lCount);
    double dphi = 2.0 * M_PI / static_cast<double>(rCount);
    // двоичный поиск вдоль шва
    for (UInteger i = 0; i < lCount; i++)
    {
        double xi0 = static_cast<double>(i) * dxi;
        double xi1 = static_cast<double>(i + 1) * dxi;
        double eta = 0.0;
        double val0 = local_func(xi0, eta);
        double val1 = local_func(xi1, eta);
        if ((val0 > epsilon_ && val1 < epsilon_) || (val0 < epsilon_ && val1 > epsilon_))
        {
            Point2D p0(xi0, 0.0);
            Point2D p1(xi1, 0.0);
            Point2D b = Mesh2D::binary(p0, p1, local_func);
            charPoints2d.push_back(Point2D(b.x(), 0.0));
            charPoints2d.push_back(Point2D(b.y(), 2.0 * M_PI));
        }
        if (fabs(val0) < epsilon_ && i > 0)
        {
            charPoints2d.push_back(Point2D(xi0, 0.0));
            charPoints2d.push_back(Point2D(xi0, 2.0 * M_PI));
        }
        if (fabs(val1) < epsilon_ && i < lCount - 1)
        {
            charPoints2d.push_back(Point2D(xi1, 0.0));
            charPoints2d.push_back(Point2D(xi1, 2.0 * M_PI));
        }
    }
    // двоичный поиск вдоль граней
    for (UInteger i = 0; i < rCount; i++)
    {
        double xi0 = 0.0;
        double xi1 = length;
        double eta0 = static_cast<double>(i) * dphi;
        double eta1 = static_cast<double>(i + 1) * dphi;
        for (double xi = xi0; xi <= xi1; xi += xi1)
        {
            double val0 = local_func(xi, eta0);
            double val1 = local_func(xi, eta1);
            if ((val0 > epsilon_ && val1 < epsilon_) || (val0 < epsilon_ && val1 > epsilon_))
            {
                Point2D p0(xi, eta0);
                Point2D p1(xi, eta1);
                Point2D b = Mesh2D::binary(p0, p1, local_func);
                charPoints2d.push_back(Point2D(xi, b.y()));
            }
            if (fabs(val0) < epsilon_ && i > 0)
            {
                charPoints2d.push_back(Point2D(xi, eta0));
            }
            if (fabs(val1) < epsilon_ && i < rCount - 1)
            {
                charPoints2d.push_back(Point2D(xi, eta1));
            }
        }
    }
    mesh2d.ruppert(lCount, rCount, -0.001, -0.001, xi_max + 0.002, 2.0 * M_PI + 0.002, func2d, charPoints2d, true);

    //    // сглаживание Лапласа, взвешенное по растояниям до соседних узлов
    //    for (int iter = 0; iter < 5; iter++)
    //    {
    //        for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
    //        {
    //            if (mesh2d.nodeType(i) == INNER)
    //            {
    //                Point2D s(0.0, 0.0);
    //                Point2D c = mesh2d.point2d(i);
    //                Point3D c3 = conePoint3d(c.x(), c.y());
    //                AdjacentSet adjacent = mesh2d.adjacent(i);
    //                double w_sum = 0.0;
    //                for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
    //                {
    //                    Triangle t = mesh2d.triangle(*it);
    //                    for (int nnode = 0; nnode < 3; nnode++)
    //                    {
    //                        if (t[nnode] != i)
    //                        {
    //                            Point2D p = mesh2d.point2d(t[nnode]);
    //                            Point3D p3 = conePoint3d(p.x(), p.y());
    //                            double w = c3.distanceTo(p3);
    //                            w_sum += w;
    //                            s = s + (w * p);
    //                        }
    //                    }
    //                }
    //                mesh2d.setPoint(i, (1.0 / w_sum) * s);
    //            }
    //        } // for i
    //    }

    for (int iter = 0; iter < 1; iter++)
    {
        for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
        {
            if (mesh2d.nodeType(i) == INNER)
            {
                //                Point2D s(0.0, 0.0);
                Point2D current = mesh2d.point2d(i);
                std::vector<double> x0(2);
                std::vector<double> result;
                //                Point3D c3 (radius * cos(c.y()), c.x(), radius * sin(c.y()));
                AdjacentSet adjacent = mesh2d.adjacent(i);
                auto functor = [&](const std::vector<double> &xxxx)
                {
                    double f = 0.0;
                    Point2D c(xxxx[0], xxxx[1]);
                    Point3D c3 = conePoint3d(xxxx[0], xxxx[1]);
                    for (AdjacentSet::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
                    {
                        Triangle t = mesh2d.triangle(*it);
                        for (int nnode = 0; nnode < 3; nnode++)
                        {
                            if (t[nnode] == i)
                            {
                                Point2D a = mesh2d.point2d(t[nnode - 1]);
                                Point2D b = mesh2d.point2d(t[nnode + 1]);
                                Point3D a3 = conePoint3d(a.x(), a.y());
                                Point3D b3 = conePoint3d(b.x(), b.y());
                                //                                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                                double alpha = minAngle(a3, c3, b3);
                                double tarea = -0.5 * ( (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y()) );
                                f += (tarea <= 0.0) ? 1000. : (1.0 / (alpha * alpha));
                                //                                break;
                            }
                        }
                    }
                    return f;
                };
                x0[0] = current.x();
                x0[1] = current.y();
                //                                std::cout << "f = " << functor(x0) << std::endl;
                result = descentGradient(functor, x0, epsilon_ * 100.0, epsilon_, 1000, false);
                //                                std::cout << "f = " << functor(result) << std::endl;
                if (functor(x0) > functor(result))mesh2d.setPoint(i, Point2D(result[0], result[1]));
                //                break;
            }
            if (i % 200 == 0) std::cout << i << " of " << mesh2d.nodesCount() << std::endl;
        } // for i
    }
    for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
    {
        Point2D p = mesh2d.point2d(i);
        Point3D p3 = conePoint3d(p.x(), p.y());
        NodeType nodeType = (fabs(func2d(p.x(), p.y())) < epsilon_ || mesh2d.nodeType(i) == CHARACTER) ? CHARACTER : BORDER;
        if (fabs(p.y() - 2.0 * M_PI) < epsilon_)
            nodes_map[i] = addNode(p3, CHARACTER);
        else
            nodes_map[i] = pushNode(p3, nodeType);
    }
    for (UInteger i = 0; i < mesh2d.elementsCount(); i++)
    {
        Triangle t = mesh2d.triangle(i);
        addElement(nodes_map[t[0]], nodes_map[t[1]], nodes_map[t[2]]);
    }
    flip(); // !
    // размеры области
    xMin_ = zMin_ = -(bottom_radius > top_radius ? bottom_radius : top_radius);
    xMax_ = zMax_ = (bottom_radius > top_radius ? bottom_radius : top_radius);
    yMin_ = 0.0;
    yMax_ = length;
    printStats();
}

void TriangleMesh3D::parametricDomain(const UInteger &uCount, const UInteger &vCount, std::function<Point3D (double, double)> domainFunction, std::function<double (double, double, double)> rfunc)
{
    clear();

    if (rfunc == nullptr)
    {
        typedef struct{
            Point2D p0;
            Point2D p1;
            Point2D p2;
        } ParametricTriangle;
        std::list<ParametricTriangle> triangles;
        double du = 1.0 / static_cast<double>(uCount - 1);
        double dv = 1.0 / static_cast<double>(vCount - 1);
        double u = 0.0;
        for (UInteger i = 0; i < uCount - 1; i++)
        {
            double v = 0.0;
            for (UInteger j = 0; j < vCount - 1; j++)
            {
                double un = (i == uCount - 2) ? 1.0 : u + du;
                double vn = (j == vCount - 2) ? 1.0 : v + dv;
                ParametricTriangle t0 = {Point2D(u, v), Point2D(un, v), Point2D(un, vn)};
                ParametricTriangle t1 = {Point2D(un, vn), Point2D(u, vn), Point2D(u, v)};
                triangles.push_back(t0);
                triangles.push_back(t1);
//                if ((0.5 <= un && 0.5 <= vn && 0.5 <= u && 0.5 <= v) || (0.5 >= un && 0.5 >= vn && 0.5 >= u && 0.5 >= v))
//                {
//                    ParametricTriangle t0 = {Point2D(u, v), Point2D(un, v), Point2D(un, vn)};
//                    ParametricTriangle t1 = {Point2D(un, vn), Point2D(u, vn), Point2D(u, v)};
//                    triangles.push_back(t0);
//                    triangles.push_back(t1);
//                } // if
//                else
//                {
//                    ParametricTriangle t0 = {Point2D(u, v), Point2D(un, v), Point2D(u, vn)};
//                    ParametricTriangle t1 = {Point2D(un, vn), Point2D(u, vn), Point2D(un, v)};
//                    triangles.push_back(t0);
//                    triangles.push_back(t1);
//                } // else
                v += dv;
            }
            u += du;
        }
        /*int count = 0;
        std::list<ParametricTriangle>::iterator t = triangles.begin();
        while ( t != triangles.end() && count < 500000)
        {
            ParametricTriangle tri = *t;
            Point3D v0 = domainFunction(tri.p0.x(), tri.p0.y());
            Point3D v1 = domainFunction(tri.p1.x(), tri.p1.y());
            Point3D v2 = domainFunction(tri.p2.x(), tri.p2.y());
            if (!v0.isEqualTo(v1) && !v1.isEqualTo(v2) && !v2.isEqualTo(v0))
            {
                Point3D c01 = 0.5 * (v0 + v1);
                Point3D b01 = domainFunction(0.5 * (tri.p0.x() + tri.p1.x()), 0.5 * (tri.p0.y() + tri.p1.y()));
                Point3D c12 = 0.5 * (v1 + v2);
                Point3D b12 = domainFunction(0.5 * (tri.p1.x() + tri.p2.x()), 0.5 * (tri.p1.y() + tri.p2.y()));
                Point3D c20 = 0.5 * (v2 + v0);
                Point3D b20 = domainFunction(0.5 * (tri.p2.x() + tri.p0.x()), 0.5 * (tri.p2.y() + tri.p0.y()));
                bool e01 = (c01.distanceTo(b01) / v0.distanceTo(v1)) >= 0.05;
                bool e12 = (c12.distanceTo(b12) / v1.distanceTo(v2)) >= 0.05;
                bool e20 = (c20.distanceTo(b20) / v2.distanceTo(v0)) >= 0.05;
                double l01 = v0.distanceTo(v1);
                double l12 = v1.distanceTo(v2);
                double l20 = v2.distanceTo(v0);
                if ((e01 && !e12 && !e20) ||
                        (e01 && e12 && !e20 && l01 >= l12) ||
                        (e01 && !e12 && e20 && l01 >= l20) ||
                        (e01 && e12 && e20 && l01 >= l12 && l01 >= l20))
                {
                    Point2D c = 0.5 * (tri.p0 + tri.p1);
                    ParametricTriangle t0 = {tri.p0, c, tri.p2}; //
                    ParametricTriangle t1 = {c, tri.p1, tri.p2};
                    triangles.erase(t);
                    triangles.push_back(t0); //
                    triangles.push_back(t1);
                    t = triangles.begin();
                    continue;
                }
                if ((!e01 && e12 && !e20) ||
                        (e01 && e12 && !e20 && l12 >= l01) ||
                        (!e01 && e12 && e20 && l12 >= l20) ||
                        (e01 && e12 && e20 && l12 >= l01 && l12 >= l20))
                {
                    Point2D c = 0.5 * (tri.p1 + tri.p2);
                    ParametricTriangle t0 = {tri.p0, tri.p1, c};
                    ParametricTriangle t1 = {tri.p0, c, tri.p2};  //
                    triangles.erase(t);
                    triangles.push_back(t0);
                    triangles.push_back(t1);//
                    t = triangles.begin();
                    continue;
                }
                if ((!e01 && !e12 && e20) ||
                        (e01 && !e12 && e20 && l20 >= l01) ||
                        (!e01 && e12 && e20 && l20 >= l12) ||
                        (e01 && e12 && e20 && l20 >= l01 && l20 >= l12))
                {
                    Point2D c = 0.5 * (tri.p2 + tri.p0);
                    ParametricTriangle t0 = {tri.p2, c, tri.p1};
                    ParametricTriangle t1 = {tri.p1, c, tri.p0};
                    triangles.erase(t);
                    triangles.push_back(t0);
                    triangles.push_back(t1);
                    t = triangles.begin();
                    continue;
                }
            }
            ++t;
            count++;
        }
        std::cout << count << std::endl;*/
        for (std::list<ParametricTriangle>::iterator t = triangles.begin(); t != triangles.end(); ++t)
        {
            ParametricTriangle tri = *t;
            UInteger p0 = addNode(domainFunction(tri.p0.x(), tri.p0.y()), BORDER);
            UInteger p1 = addNode(domainFunction(tri.p1.x(), tri.p1.y()), BORDER);
            UInteger p2 = addNode(domainFunction(tri.p2.x(), tri.p2.y()), BORDER);
            if (p0 != p1 && p1 != p2 && p0 != p2) addElement(p0, p1, p2);
        }
        updateDomain();
        printStats();
        return;
    }

    auto distance = [&](Point2D p0, Point2D p1)
    {
        Point3D pp0 = domainFunction(p0.x(), p0.y());
        Point3D pp1 = domainFunction(p1.x(), p1.y());
        return pp0.distanceTo(pp1);
    };

    auto local_func = [&](double u, double v)
    {
        Point3D p3 = domainFunction(u, v);
        return (rfunc != nullptr) ? rfunc(p3.x(), p3.y(), p3.z()) : 1.0;
    };

    auto func2d = [&](double u, double v)
    {
        double f = local_func(u, v);
        double r = con(0.5 * 0.5 - (u - 0.5) * (u - 0.5), 0.5 * 0.5 - (v - 0.5) * (v - 0.5));
        return con(f, r);
    };
    std::list<Point2D> charPoints2d;
    TriangleMesh2D mesh2d;

    if (local_func(0.0, 0.0) > -epsilon_ )
    {
        charPoints2d.push_back(Point2D(0.0, 0.0));
    }
    if (local_func(1.0, 0.0) > -epsilon_ )
    {
        charPoints2d.push_back(Point2D(1.0, 0.0));
    }
    if (local_func(1.0, 1.0) > -epsilon_ )
    {
        charPoints2d.push_back(Point2D(1.0, 1.0));
    }
    if (local_func(0.0, 1.0) > -epsilon_ )
    {
        charPoints2d.push_back(Point2D(0.0, 1.0));
    }
    double du = 1.0 / static_cast<double>(uCount);
    double dv = 1.0 / static_cast<double>(vCount);
    // двоичный поиск по первому направлению
    for (UInteger i = 0; i < uCount; i++)
    {
        double u0 = static_cast<double>(i) * du;
        double u1 = static_cast<double>(i + 1) * du;
        double v[] = {0.0, 1.0};
        for (int j = 0; j < 2; j++)
        {
            double val0 = local_func(u0, v[j]);
            double val1 = local_func(u1, v[j]);
            if (fabs(val0) > epsilon_ && fabs(val1) > epsilon_ && val0 * val1 < 0.0)
            {
                Point2D p0(u0, v[j]);
                Point2D p1(u1, v[j]);
                Point2D b = Mesh2D::binary(p0, p1, local_func);
                charPoints2d.push_back(Point2D(b.x(), b.y()));
            }
            if (fabs(val0) < epsilon_ && i > 0)
            {
                charPoints2d.push_back(Point2D(u0, v[j]));
            }
            if (fabs(val1) < epsilon_ && i < uCount - 1)
            {
                charPoints2d.push_back(Point2D(u1, v[j]));
            }
        }
    }
    // двоичный поиск по второму направлению
    for (UInteger i = 0; i < vCount; i++)
    {
        double u[] = {0.0, 1.0};
        double v0 = static_cast<double>(i) * dv;
        double v1 = static_cast<double>(i + 1) * dv;
        for (int j = 0; j < 2; j++)
        {
            double val0 = local_func(u[j], v0);
            double val1 = local_func(u[j], v1);
            if (fabs(val0) > epsilon_ && fabs(val1) > epsilon_ && val0 * val1 < 0.0)
            {
                Point2D p0(u[j], v0);
                Point2D p1(u[j], v1);
                Point2D b = Mesh2D::binary(p0, p1, local_func);
                charPoints2d.push_back(Point2D(u[j], b.y()));
            }
            if (fabs(val0) < epsilon_ && i > 0)
            {
                charPoints2d.push_back(Point2D(u[j], v0));
            }
            if (fabs(val1) < epsilon_ && i < vCount - 1)
            {
                charPoints2d.push_back(Point2D(u[j], v1));
            }
        }
    }
    SegmentMesh2D smesh;
    smesh.MarchingQuads(uCount, vCount, 0.0, 0.0, 1.0, 1.0, func2d, charPoints2d, 0.0, 0, 0, distance);
    TriangleMesh2D::Triangulation triangulation = TriangleMesh2D::superDelaunay(&smesh, func2d);
    TriangleMesh2D::superRuppert(triangulation, &smesh, func2d);
    std::list<Triangle>::iterator triangle = triangulation.triangles.begin();
    UInteger niter = 0;
    while (triangle != triangulation.triangles.end() && niter < 4294967290UL)
    {
        if (triangle->vertexNode(0) > 3 && triangle->vertexNode(1) > 3 && triangle->vertexNode(2) > 3)
        {
            Point2D A = triangulation.nodes[triangle->vertexNode(0)];
            Point2D B = triangulation.nodes[triangle->vertexNode(1)];
            Point2D C = triangulation.nodes[triangle->vertexNode(2)];
            Point3D A3 = domainFunction(A.x(), A.y());
            Point3D B3 = domainFunction(B.x(), B.y());
            Point3D C3 = domainFunction(C.x(), C.y());
            //            // 3d центры
            Point3D AB3 = domainFunction(0.5 * (A.x() + B.x()), 0.5 * (A.y() + B.y()));
            Point3D AC3 = domainFunction(0.5 * (A.x() + C.x()), 0.5 * (A.y() + C.y()));
            Point3D BC3 = domainFunction(0.5 * (B.x() + C.x()), 0.5 * (B.y() + C.y()));
            double bc = distance(B, C);
            double ac = distance(A, C);
            double ab = distance(A, B);
            //            if (minAngle(A3, B3, C3) < 0.3 && (0.0 <= A.x() && A.x() <= 1.0 && 0.0 <= B.x() && B.x() <= 1.0 && 0.0 <= C.x() && C.x() <= 1.0)
            //                    && (0.0 <= A.y() && A.y() <= 1.0 && 0.0 <= B.y() && B.y() <= 1.0 && 0.0 <= C.y() && C.y() <= 1.0))
            if ((AB3.distanceTo(0.5 * (A3 + B3)) / ab > 0.04 || AC3.distanceTo(0.5 * (A3 + C3)) / ac > 0.04 || BC3.distanceTo(0.5 * (B3 + C3)) / bc > 0.04) && (func2d(A.x(), A.y()) >= -epsilon_ || func2d(B.x(), B.y()) >= -epsilon_ || func2d(C.x(), C.y()) >= -epsilon_))
            {
                Point2D center;
                double xc = 0.0, yc = 0.0, r = 0.0;
                TriangleMesh2D::circumCircle(0.0, 0.0, A.x(), A.y(), B.x(), B.y(), C.x(), C.y(), xc, yc, r);
                center.set(xc, yc);
                //                center = circumCylinderCenter(A, B, C, domainFunction);
                UInteger seg_num = 0;
                if (smesh.isEncroached(center, seg_num))
                {
                    Point2D R = smesh.refineMidpoint(seg_num, func2d);
                    if (TriangleMesh2D::insertDelaunayNode(R, INNER, triangulation))
                        triangle = triangulation.triangles.begin();
                    else
                    {
                        std::cout << R.x() << " " << R.y() << std::endl;
                        ++triangle;
                    }
                }
                else if (TriangleMesh2D::insertDelaunayNode(center, INNER, triangulation))
                {
                    triangle = triangulation.triangles.begin();
                    std::cout << niter << ' ';
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

        ++niter;
    }
    std::cout << "Ruppert's niter: " << niter << std::endl;

    for (std::list<Triangle>::iterator triangle = triangulation.triangles.begin(); triangle != triangulation.triangles.end(); ++triangle)
    {
        Point2D A = triangulation.nodes[triangle->vertexNode(0)];
        Point2D B = triangulation.nodes[triangle->vertexNode(1)];
        Point2D C = triangulation.nodes[triangle->vertexNode(2)];
        Point2D center = 1.0 / 3.0 * (A + B + C);
        double val_a = func2d(A.x(), A.y());
        double val_b = func2d(B.x(), B.y());
        double val_c = func2d(C.x(), C.y());
        double val_center = func2d(center.x(), center.y());
        if (val_a > -epsilon_ && val_b > -epsilon_ && val_c > -epsilon_ && val_center > -epsilon_)
            addElement(addNode(domainFunction(A.x(), A.y()), BORDER), addNode(domainFunction(B.x(), B.y()), BORDER), addNode(domainFunction(C.x(), C.y()), BORDER));
    }
    flip();
    printStats();
    updateDomain();
}

void TriangleMesh3D::marchingCubes(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth, std::function<double (double, double, double)> func, double level, bool slice_x, bool slice_y, bool slice_z, int smooth, int optimize, bool useFlip)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    zMin_ = zMin;
    zMax_ = zMin + depth;
    int edges[256] = { 0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
                       0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
                       0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
                       0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
                       0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
                       0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
                       0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
                       0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
                       0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
                       0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
                       0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
                       0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
                       0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
                       0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
                       0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
                       0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
                       0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
                       0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
                       0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
                       0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
                       0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
                       0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
                       0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
                       0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
                       0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
                       0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
                       0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
                       0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
                       0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
                       0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
                       0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
                       0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0 };
    int triangles[256][16] = {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
                              {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
                              {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
                              {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
                              {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
                              {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
                              {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
                              {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
                              {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
                              {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
                              {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
                              {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
                              {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
                              {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
                              {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
                              {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
                              {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
                              {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
                              {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
                              {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
                              {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
                              {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
                              {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
                              {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
                              {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
                              {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
                              {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
                              {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
                              {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
                              {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
                              {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
                              {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
                              {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
                              {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
                              {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
                              {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
                              {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
                              {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
                              {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
                              {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
                              {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
                              {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
                              {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
                              {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
                              {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
                              {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
                              {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
                              {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
                              {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
                              {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
                              {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
                              {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
                              {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
                              {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
                              {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
                              {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
                              {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
                              {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
                              {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
                              {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
                              {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
                              {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
                              {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
                              {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
                              {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
                              {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
                              {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
                              {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
                              {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
                              {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
                              {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
                              {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
                              {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
                              {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
                              {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
                              {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
                              {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
                              {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
                              {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
                              {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
                              {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
                              {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
                              {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
                              {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
                              {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
                              {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
                              {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
                              {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
                              {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
                              {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
                              {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
                              {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
                              {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
                              {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
                              {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
                              {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
                              {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
                              {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
                              {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
                              {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
                              {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
                              {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
                              {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
                              {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
                              {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
                              {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
                              {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
                              {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
                              {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
                              {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
                              {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
                              {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
                              {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
                              {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
                              {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
                              {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
                              {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
                              {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
                              {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
                              {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
                              {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
                              {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
                              {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
                              {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
                              {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
                              {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
                              {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
                              {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
                              {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
                              {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
                              {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
                              {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
                              {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
                              {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
                              {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
                              {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
                              {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
                              {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
                              {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
                              {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
                              {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
                              {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
                              {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
                              {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
                              {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
                              {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
                              {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
                              {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
                              {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
                              {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
                              {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
                              {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
                              {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
                              {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
                              {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
                              {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
                              {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
                              {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
                              {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
                              {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
                              {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
                              {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
                              {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
                              {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
                              {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
                              {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
                              {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
                              {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
                              {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
                              {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
                              {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
                              {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
                              {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
                              {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
                              {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
                              {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
                              {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
                              {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                              {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1} };
    double hx = width / static_cast<double>(xCount - 1);
    double hy = height / static_cast<double>(yCount - 1);
    double hz = depth / static_cast<double>(zCount - 1);
    double xc = xMin + width / 2.0;
    double yc = yMin + height / 2.0;
    double zc = zMin + depth / 2.0;

    std::list<CooTriangle> cootriangles;
    //
    auto func_slice = [&](double x, double y, double z)
    {
        double r = func(x, y, z);
        if (slice_x) r = con(r, xc - x);
        if (slice_y) r = con(r, yc - y);
        if (slice_z) r = con(r, zc - z);
        return r;
    };
    //
    ConsoleProgress progress(xCount - 1);
    double x = xMin;
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        double y = yMin;
        ++progress;
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            double z = zMin;
            for (UInteger k = 0; k < zCount - 1; k++)
            {
                int index = 0;
                Point3D border[12];
                Point3D p[8];
                p[0].set(x, y, z);
                p[1].set(x + hx, y, z);
                p[2].set(x + hx, y + hy, z);
                p[3].set(x, y + hy, z);
                p[4].set(x, y, z + hz);
                p[5].set(x + hx, y, z + hz);
                p[6].set(x + hx, y + hy, z + hz);
                p[7].set(x, y + hy, z + hz);
                if (func_slice(p[0].x(), p[0].y(), p[0].z()) - level < epsilon_) index |= 1;
                if (func_slice(p[1].x(), p[1].y(), p[1].z()) - level < epsilon_) index |= 2;
                if (func_slice(p[2].x(), p[2].y(), p[2].z()) - level < epsilon_) index |= 4;
                if (func_slice(p[3].x(), p[3].y(), p[3].z()) - level < epsilon_) index |= 8;
                if (func_slice(p[4].x(), p[4].y(), p[4].z()) - level < epsilon_) index |= 16;
                if (func_slice(p[5].x(), p[5].y(), p[5].z()) - level < epsilon_) index |= 32;
                if (func_slice(p[6].x(), p[6].y(), p[6].z()) - level < epsilon_) index |= 64;
                if (func_slice(p[7].x(), p[7].y(), p[7].z()) - level < epsilon_) index |= 128;
                if (index != 0)
                {
                    if (edges[index] & 1) border[0] = binary(p[0], p[1], func_slice, level);
                    if (edges[index] & 2) border[1] = binary(p[1], p[2], func_slice, level);
                    if (edges[index] & 4) border[2] = binary(p[2], p[3], func_slice, level);
                    if (edges[index] & 8) border[3] = binary(p[3], p[0], func_slice, level);
                    if (edges[index] & 16) border[4] = binary(p[4], p[5], func_slice, level);
                    if (edges[index] & 32) border[5] = binary(p[5], p[6], func_slice, level);
                    if (edges[index] & 64) border[6] = binary(p[6], p[7], func_slice, level);
                    if (edges[index] & 128) border[7] = binary(p[7], p[4], func_slice, level);
                    if (edges[index] & 256) border[8] = binary(p[0], p[4], func_slice, level);
                    if (edges[index] & 512) border[9] = binary(p[1], p[5], func_slice, level);
                    if (edges[index] & 1024) border[10] = binary(p[2], p[6], func_slice, level);
                    if (edges[index] & 2048) border[11] = binary(p[3], p[7], func_slice, level);
                    for (int t = 0; triangles[index][t] != -1; t += 3)
                    {
                        CooTriangle coot;
                        coot.a = border[triangles[index][t]];
                        coot.b = border[triangles[index][t + 1]];
                        coot.c = border[triangles[index][t + 2]];
                        cootriangles.push_back(coot);
                    }
                }
                z += hz;
            }
            y += hy;
        }
        x += hx;
    }
    double delta = 0.2 * sqrt(hx*hx + hy*hy + hz*hz);
//    std::cout << cootriangles.size() << std::endl;
    clearCooTriangles(cootriangles, func_slice, level, delta);
    for (CooTriangle t: cootriangles)
    {
        addElement(addNode(t.a, BORDER, delta), addNode(t.b, BORDER, delta), addNode(t.c, BORDER, delta));
    }
    if (useFlip)
        flip();

    laplacianSmoothing(func_slice, level, smooth, useFlip);

    distlenSmoothing(func_slice, level, optimize, useFlip);

    evalNodalValues(func);
    printStats();
}

void TriangleMesh3D::marchingTetrahedrons(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth, std::function<double (double, double, double)> func, double level, bool slice_x, bool slice_y, bool slice_z, int smooth, int optimize, bool useFlip)
{
    clear();
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    zMin_ = zMin;
    zMax_ = zMin + depth;
    double hx = width / static_cast<double>(xCount - 1);
    double hy = height / static_cast<double>(yCount - 1);
    double hz = depth / static_cast<double>(zCount - 1);
    double xc = xMin + width / 2.0;
    double yc = yMin + height / 2.0;
    double zc = zMin + depth / 2.0;
    // номера вершин куба, образующих шесть тетраэдров
    const int tetrahedronsInACube[6][4] =
    {
        {0,5,1,6},
        {0,1,2,6},
        {0,2,3,6},
        {0,3,7,6},
        {0,7,4,6},
        {0,4,5,6}
    };
    // коды ребер тетрыэдра, которые пересечены границей
    const int tetrahedronEdgeFlags[16]=
    {
        0x00, 0x0d, 0x13, 0x1e, 0x26, 0x2b, 0x35, 0x38, 0x38, 0x35, 0x2b, 0x26, 0x1e, 0x13, 0x0d, 0x00
    };
    // номера пар вершин, образующих ребра тетраэра
    const int tetrahedronEdgeConnection[6][2] =
    {
        {0,1},  {1,2},  {2,0},  {0,3},  {1,3},  {2,3}
    };
    // список треугольников
    const int tetrahedronTriangles[16][7] =
    {
            {-1, -1, -1, -1, -1, -1, -1},
            { 0,  3,  2, -1, -1, -1, -1},
            { 0,  1,  4, -1, -1, -1, -1},
            { 1,  4,  2,  2,  4,  3, -1},
            { 1,  2,  5, -1, -1, -1, -1},
            { 0,  3,  5,  0,  5,  1, -1},
            { 0,  2,  5,  0,  5,  4, -1},
            { 5,  4,  3, -1, -1, -1, -1},
            { 3,  4,  5, -1, -1, -1, -1},
            { 4,  5,  0,  5,  2,  0, -1},
            { 1,  5,  0,  5,  3,  0, -1},
            { 5,  2,  1, -1, -1, -1, -1},
            { 3,  4,  2,  2,  4,  1, -1},
            { 4,  1,  0, -1, -1, -1, -1},
            { 2,  3,  0, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1},
    };
    std::list<CooTriangle> cootriangles;
    // lambda-функция для построенгия сечений области
    auto func_slice = [&](double x, double y, double z)
    {
        double r = func(x, y, z);
        if (slice_x) r = con(r, xc - x);
        if (slice_y) r = con(r, yc - y);
        if (slice_z) r = con(r, zc - z);
        return r;
    };
    //
    ConsoleProgress progress(xCount - 1);
    double x = xMin;
    for (UInteger i = 0; i < xCount - 1; i++)
    {
        double y = yMin;
        ++progress;
        for (UInteger j = 0; j < yCount - 1; j++)
        {
            double z = zMin;
            for (UInteger k = 0; k < zCount - 1; k++)
            {
                Point3D p[8];
                double cubeValue[8];
                double tetrahedronValue[4];
                Point3D tetrahedronPosition[4];
                p[0].set(x, y, z);
                p[1].set(x + hx, y, z);
                p[2].set(x + hx, y + hy, z);
                p[3].set(x, y + hy, z);
                p[4].set(x, y, z + hz);
                p[5].set(x + hx, y, z + hz);
                p[6].set(x + hx, y + hy, z + hz);
                p[7].set(x, y + hy, z + hz);
                for (int icv = 0; icv < 8; icv++)
                    cubeValue[icv] = func_slice(p[icv].x(), p[icv].y(), p[icv].z()) - level;
                for (int iTetrahedron = 0; iTetrahedron < 6; iTetrahedron++)
                {
                    for (int iVertex = 0; iVertex < 4; iVertex++)
                    {
                        int iVertexInACube = tetrahedronsInACube[iTetrahedron][iVertex];
                        tetrahedronPosition[iVertex] = p[iVertexInACube];
                        tetrahedronValue[iVertex] = cubeValue[iVertexInACube];
                    }
                    int iFlagIndex = 0;
                    Point3D edgeVertex[6];
                    for (int iVertex = 0; iVertex < 4; iVertex++)
                    {
                        if(tetrahedronValue[iVertex] < epsilon_)
                            iFlagIndex |= 1<<iVertex; // если внешний или граничный
                    }
                    int iEdgeFlags = tetrahedronEdgeFlags[iFlagIndex];
                    if (iEdgeFlags == 0)
                        continue;
                    for (int iEdge = 0; iEdge < 6; iEdge++)
                    {
                        //if there is an intersection on this edge
                        if(iEdgeFlags & (1<<iEdge))
                        {
                            int iVert0 = tetrahedronEdgeConnection[iEdge][0];
                            int iVert1 = tetrahedronEdgeConnection[iEdge][1];
                            edgeVertex[iEdge] = binary(tetrahedronPosition[iVert0], tetrahedronPosition[iVert1], func_slice, level);
                        }
                    }
                    for (int iTriangle = 0; iTriangle < 2 && tetrahedronTriangles[iFlagIndex][3*iTriangle] >= 0; iTriangle++)
                    {
                        Point3D t0 = edgeVertex[tetrahedronTriangles[iFlagIndex][3*iTriangle]];
                        Point3D t1 = edgeVertex[tetrahedronTriangles[iFlagIndex][3*iTriangle+1]];
                        Point3D t2 = edgeVertex[tetrahedronTriangles[iFlagIndex][3*iTriangle+2]];
                        CooTriangle coot;
                        coot.a = t0;
                        coot.b = t1;
                        coot.c = t2;
                        cootriangles.push_back(coot);
                        /*if (!t0.isEqualTo(t1, epsilon_) && !t0.isEqualTo(t2, epsilon_) && !t1.isEqualTo(t2, epsilon_))
                        {
                            addElementOrdered(t2, t1, t0);
                        }*/
                    }
                }
                z += hz;
            }
            y += hy;
        }
        x += hx;
    }
    double delta = 0.2 * sqrt(hx*hx + hy*hy + hz*hz);
    clearCooTriangles(cootriangles, func_slice, level, delta);
    for (CooTriangle t: cootriangles)
    {
        addElement(addNode(t.a, BORDER, delta), addNode(t.b, BORDER, delta), addNode(t.c, BORDER, delta));
    }
    if (useFlip)
        flip();

    laplacianSmoothing(func_slice, level, smooth, useFlip);

    distlenSmoothing(func_slice, level, optimize, useFlip);

    evalNodalValues(func);
    printStats();
}

std::list<ElementPointer> TriangleMesh3D::backgroundGrid(const TetrahedralMesh3D *mesh, std::function<double (double, double, double)> func, double level, int smooth, int optimize, bool useFlip)
{
    clear();
    std::list<ElementPointer> inner;
    double h = mesh->xMax() - mesh->xMin();
    ConsoleProgress progress(mesh->elementsCount());
    // make a list of background elements and a set of border points in this list
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        ElementPointer el = mesh->element(i);
        int m = el->verticesCount();
        std::vector<double> values(m);
        int code = 0;
        for (int j = 0; j < m; j++)
        {
            Point3D point = mesh->point3d(el->vertexNode(j));
            double value = func(point.x(), point.y(), point.z());
            values[j] = value;
            if (value - level < epsilon_)
                code |= (1 << j);
        }
        if (code == 0)
        {
            inner.push_back(el);
        }
        ++progress;
    }
    // delete all isolated elements and add border faces into the mesh
    progress.restart(inner.size());
    for (std::list<ElementPointer>::iterator it = inner.begin(); it != inner.end();)
    {
        bool isWeak = false;
        UInteger weak0, weak1;
        ElementPointer el = *it;
        for (int j = 0; j < el->facesCount(); j++)
        {
            UIntegerVector f = el->face(j);
            Point3D p[3];
            p[0] = mesh->point3d(f[0]);
            p[1] = mesh->point3d(f[1]);
            p[2] = mesh->point3d(f[2]);
            double d = p[0].distanceTo(p[1]);
            bool isInner = false;
            // check an other element with these vertices
            for (ElementPointer el1: inner)
            {
                if (el1 != el && el1->in(f[0]) && el1->in(f[1]) && el1->in(f[2]))
                {
                    isInner = true;
                    break;
                }
            }
            if (!isInner)
            {
                UInteger index0 = addNode(p[0], BORDER);
                UInteger index1 = addNode(p[1], BORDER);
                UInteger index2 = addNode(p[2], BORDER);
                addElement(index0, index1, index2);
                AdjacentSet a0 = node_[index0].adjacent;
                AdjacentSet a1 = node_[index1].adjacent;
                AdjacentSet a2 = node_[index2].adjacent;
                std::vector<UInteger> common01;
                std::vector<UInteger> common12;
                std::vector<UInteger> common20;
                set_intersection(a0.begin(), a0.end(), a1.begin(), a1.end(), std::back_inserter(common01));
                set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), std::back_inserter(common12));
                set_intersection(a2.begin(), a2.end(), a0.begin(), a0.end(), std::back_inserter(common20));
//                if (common01.size() > 2 || common12.size() > 2 || common20.size() > 2)
//                {
//                    isWeak = true;
//                    break;
//                }
                if (common01.size() > 2)
                {
                    isWeak = true;
                    weak0 = f[0];
                    weak1 = f[1];
                    break;
                }
                if (common12.size() > 2)
                {
                    isWeak = true;
                    weak0 = f[1];
                    weak1 = f[2];
                    break;
                }
                if (common20.size() > 2)
                {
                    isWeak = true;
                    weak0 = f[2];
                    weak1 = f[0];
                    break;
                }
            }
            if (d < h)
                h = d;
        }
        if (isWeak)
        {
            UInteger power = 1000000;
            std::list<ElementPointer>::iterator min_ptr = inner.end();
            for (std::list<ElementPointer>::iterator iit = inner.begin(); iit != inner.end(); iit++)
            {
                ElementPointer ell = *iit;
                if (ell->in(weak0) || ell->in(weak1))
                {
                    UInteger p = 0;
                    for (int j = 0; j < ell->verticesCount(); j++)
                    {
                        for (ElementPointer el1: inner)
                        {
                            if (el1->in(ell->vertexNode(j))) p++;
                        }
                    }
//                        p += mesh->adjacentCount(ell->vertexNode(j));
                    if (p < power)
                    {
                        power = p;
                        min_ptr = iit;

                    }
                }
            }
            if (min_ptr != inner.end())
            {
                inner.erase(min_ptr);
                it = inner.begin();
                progress.restart(inner.size());
                clear();
            }
            else it++;
        }
        else
        {
            it++;
            ++progress;
        }
    }
    std::vector<Point3D> normals(node_.size());
    for (UInteger i = 0; i != node_.size(); ++i)
    {
        Point3D point = node_[i].point;
        AdjacentSet adjacent = node_[i].adjacent;
        Point3D n(0.0, 0.0, 0.0);
        for (UInteger elnum: adjacent)
        {
            Triangle element = element_[elnum];
            int index = element.index(i);
            Point3D prev = node_[element[index - 1]].point;
            Point3D next = node_[element[index + 1]].point;
            n = n + normal3(point, prev, next);
        }
        normals[i] = n.normalized();
    }
    for (short j = 0; j < 10; j++)
    {
        for (UInteger i = 0; i != node_.size(); ++i)
        {
            AdjacentSet adjacent = node_[i].adjacent;
            Point3D n = normals[i];
            AdjacentSet neighbours;
            for (UInteger elnum: adjacent)
            {
                Triangle element = element_[elnum];
                int index = element.index(i);
                neighbours.insert(element[index + 1]);
                neighbours.insert(element[index - 1]);
            }
            for (UInteger npointer: neighbours)
            {
                n = n + normals[npointer];
            }
            normals[i] = n.normalized();
        }
    }
    // loop nodes of the mesh and find correct position
    progress.restart(node_.size());
    std::vector<Point3D> surface(node_.size());
    for (UInteger i = 0; i != node_.size(); ++i)
    {
        surface[i] = findBorder(node_[i].point, normals[i], func, 0.5 * h, level);
        ++progress;
    }
    for (UInteger i = 0; i != node_.size(); ++i) node_[i].point = surface[i];
    laplacianSmoothing(func, level, smooth, useFlip);
    distlenSmoothing(func, level, optimize, useFlip);
//    std::list<UInteger> ee;
//    ee.push_back(0); ee.push_back(1); ee.push_back(2);ee.push_back(3);ee.push_back(4);
//    subdivide(ee, func);
//    subdivide(ee, func);
    printStats();
    updateDomain();
    return  inner;
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

void TriangleMesh3D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0], nodes_ref[1], nodes_ref[2]);
}

void TriangleMesh3D::addElementOrdered(const Point3D &t0, const Point3D &t1, const Point3D &t2, double epsilon)
{
    UInteger p0;
    UInteger p1;
    UInteger p2;
    if (t0 <= t1 && t0 <= t2)
    {
        p0 = addNode(t0, BORDER, epsilon);
        if (t1 <= t2)
        {
            p1 = addNode(t1, BORDER, epsilon);
            p2 = addNode(t2, BORDER, epsilon);
        }
        else
        {
            p2 = addNode(t2, BORDER, epsilon);
            p1 = addNode(t1, BORDER, epsilon);
        }
    }
    else if (t1 <= t0 && t1 <= t2)
    {
        p1 = addNode(t1, BORDER, epsilon);
        if (t0 <= t2)
        {
            p0 = addNode(t0, BORDER, epsilon);
            p2 = addNode(t2, BORDER, epsilon);
        }
        else
        {
            p2 = addNode(t2, BORDER, epsilon);
            p0 = addNode(t0, BORDER, epsilon);
        }
    }
    else
    {
        p2 = addNode(t2, BORDER, epsilon);
        if (t1 <= t0)
        {
            p1 = addNode(t1, BORDER, epsilon);
            p0 = addNode(t0, BORDER, epsilon);
        }
        else
        {
            p0 = addNode(t0, BORDER, epsilon);
            p1 = addNode(t1, BORDER, epsilon);
        }
    }
    addElement(p0, p1, p2);
}

double TriangleMesh3D::area(const UInteger &number) const
{
    Triangle triangle = element_[number];
    return area(node_[triangle[0]].point, node_[triangle[1]].point, node_[triangle[2]].point);
}

double TriangleMesh3D::area(const Point3D &p0, const Point3D &p1, const Point3D &p2) const
{
    double a = p0.distanceTo(p1);
    double b = p1.distanceTo(p2);
    double c = p2.distanceTo(p0);
    double p = (a + b + c) / 2.0;
    return sqrt(p * (p - a) * (p - b) * (p - c)); // формула Герона
}

double TriangleMesh3D::minAngle(const UInteger &elnum) const
{
    const Triangle tri = element_[elnum];
    const Point3D p0 = node_[tri[0]].point;
    const Point3D p1 = node_[tri[1]].point;
    const Point3D p2 = node_[tri[2]].point;
    return minAngle(p0, p1, p2);
}

double TriangleMesh3D::maxAngle(const UInteger &elnum) const
{
    const Triangle tri = element_[elnum];
    const Point3D A = node_[tri[0]].point;
    const Point3D B = node_[tri[1]].point;
    const Point3D C = node_[tri[2]].point;
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    angles(A, B, C, alpha, beta, gamma);
    return std::max(alpha, std::max(beta, gamma));
}

void TriangleMesh3D::clearElements()
{
    element_.clear();
}

void TriangleMesh3D::add(const TriangleMesh3D *mesh)
{
    UInteger count = 0;
    std::vector<UInteger> nodesPointers(mesh->nodesCount());
    ConsoleProgress progress_bar(mesh->nodesCount());
    for (UInteger i = 0; i < mesh->nodesCount(); i++)
    {
        UInteger nc = nodesCount();
        nodesPointers[i] = addNode(mesh->node_[i]);
        if (nc <= nodesPointers[i]) count++;
        ++progress_bar;
    }
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        Triangle triangle = mesh->element_[i];
        addElement(nodesPointers[triangle[0]], nodesPointers[triangle[1]], nodesPointers[triangle[2]]);
    }
    std::cout << "Added " << count << " new nodes." << std::endl;
    updateDomain();
}

double TriangleMesh3D::cfunction(const double &x, const double &y, const double &z)
{
    double min_distance = std::numeric_limits<double>::max();
    double sign = +1.0;
    UInteger count = 0;
    Point3D point(x, y, z);
    for (std::vector<Triangle>::iterator t = element_.begin(); t != element_.end(); ++t)
    {
        Triangle triangle = *t;
        Point3D p0 = node_[triangle[0]].point;
        Point3D p1 = node_[triangle[1]].point;
        Point3D p2 = node_[triangle[2]].point;
        double d = point.distanceTo(p0, p1, p2);
        if (d < min_distance) min_distance = d;

        if ((p0.z() < z && (z <= p1.z() || z <= p2.z())) || (p1.z() < z && (z <= p2.z() || z <= p0.z())) || (p2.z() < z && (z <= p0.z() || z <= p1.z())))
        {
            std::vector<Point3D> plane;
            if (fabs(z - p0.z()) < epsilon_) plane.push_back(p0);
            if (fabs(z - p1.z()) < epsilon_) plane.push_back(p1);
            if (fabs(z - p2.z()) < epsilon_) plane.push_back(p2);
            if (plane.size() < 2)
            {
                // треугольник не лежит в плоскости Z = z
                if (signbit(z - p0.z()) != signbit(z - p1.z()) && !(fabs(z - p0.z()) < epsilon_) && !(fabs(z - p1.z()) < epsilon_))
                {
                    double t = (z - p0.z()) / (p1.z() - p0.z());
                    Point3D cross = p0 + t * (p1 - p0);
                    plane.push_back(cross);
                }
                if (signbit(z - p1.z()) != signbit(z - p2.z()) && !(fabs(z - p1.z()) < epsilon_) && !(fabs(z - p2.z()) < epsilon_))
                {
                    double t = (z - p1.z()) / (p2.z() - p1.z());
                    Point3D cross = p1 + t * (p2 - p1);
                    plane.push_back(cross);
                }
                if (signbit(z - p2.z()) != signbit(z - p0.z()) && !(fabs(z - p2.z()) < epsilon_) && !(fabs(z - p0.z()) < epsilon_))
                {
                    double t = (z - p2.z()) / (p0.z() - p2.z());
                    Point3D cross = p2 + t * (p0 - p2);
                    plane.push_back(cross);
                }
            }
            if (plane.size() == 2)
            {
                Point2D q1(plane[0].x(), plane[0].y());
                Point2D q2(plane[1].x(), plane[1].y());
                double tp = -1.0, tq = -1.0;
                if ( isCrossed(Point2D(x, y), Point2D(x + 1.0, y), q1, q2, tp, tq) && tp > 0.0 && ((q1.y() < y && y <= q2.y()) || (q2.y() < y && y <= q1.y())) )
                    count++;
            }
        }
    }
    sign = (count % 2 == 0 ) ? -1.0 : 1.0;
    return sign * min_distance;
}

void TriangleMesh3D::laplacianSmoothing(std::function<double (double, double, double)> func, double level, int iter_num, bool useFlip)
{
//    auto functor = [&](const AdjacentSet &adjasentset, const UInteger &nnode)
//    {
//        const double alpha = 0.001;
//        const double beta = 1.0 - alpha;
//        double F = 0.0; // функционал
//        for (UInteger adj: adjasentset)
//        {
//            Triangle t = element_[adj];
//            int index = t.index(nnode);
//            Point3D A = node_[t[index]].point;
//            Point3D B = node_[t[index + 1]].point;
//            Point3D C = node_[t[index + 2]].point;
//            Point3D AB(A, B);
//            Point3D AC(A, C);
//            double ab = AB.length();
//            double ac = AC.length();
//            Point3D center = (A + B + C) / 3.0;
//            Point3D p = findBorder(center, func, min(ab, ac) * 0.25, level);
//            double d2b = center.distanceTo(p);
//            F += alpha * (ab*ab + ac*ac) + beta * d2b*d2b;
//        }
//        return F;
//    };
    std::cout << "Laplacian Smoothing: nodes - " << nodesCount() << ", elements - " << elementsCount() << "." << std::endl;
    for (short iit = 0; iit < iter_num; iit++)
    {
        ConsoleProgress progress(nodesCount());
        for (UInteger nnode = 0; nnode < nodesCount(); nnode++)
        {
            AdjacentSet adjacent = node_[nnode].adjacent;
//            Point3D prev = node_[nnode].point;
            Point3D point(0.0, 0.0, 0.0);
            AdjacentSet neighbours;
            double avr_len = 0.0; // средняя длина ребра
//            double f = functor(adjacent, nnode);
            for (UInteger elnum: adjacent)
            {
                Triangle triangle = element_[elnum];
                int index = triangle.index(nnode);
                neighbours.insert(triangle[index + 1]);
                neighbours.insert(triangle[index + 2]);
            }
            for (UInteger npointer: neighbours)
            {
                point = point + node_[npointer].point;
            }
            point.scale(1.0 / static_cast<double>(neighbours.size()));
            for (UInteger elnum: adjacent)
            {
                Triangle triangle = element_[elnum];
                int index = triangle.index(nnode);
                Point3D a(point, node_[triangle[index - 1]].point);
                Point3D b(point, node_[triangle[index + 1]].point);
                avr_len += 0.5 * (a.length() + b.length());
            }
            avr_len /= static_cast<double>(adjacent.size());
            node_[nnode].point = findBorder(point, func, 0.1 * avr_len, level);
//            if (f < functor(adjacent, nnode)) node_[nnode].point = prev;
            ++progress;
        }
        if (useFlip)
            flip();
    }
}

void TriangleMesh3D::distlenSmoothing(std::function<double (double, double, double)> func, double level, int iter_num, bool useFlip)
{
    auto functor = [&](const AdjacentSet &adjasentset, const UInteger &nnode)
    {
        // alpha = 0.001 is optimal for gradient searching
        const double alpha = 0.001;
        const double beta = 1.0 - alpha;
        double F = 0.0; // функционал
        for (UInteger adj: adjasentset)
        {
            Triangle t = element_[adj];
            int index = t.index(nnode);
            Point3D A = node_[t[index]].point;
            Point3D B = node_[t[index + 1]].point;
            Point3D C = node_[t[index + 2]].point;
//            double d2b = distToBorder(A, B, C, func, 0.33333333, 0.33333333, level);
            Point3D AB(A, B);
            Point3D AC(A, C);
            double ab = AB.length();
            double ac = AC.length();
            double h = std::min(ab, ac) * 0.25;
            Point3D center = (A + B + C) / 3.0;
            Point3D p = findBorder(center, func, h, level);
            double d2b = center.distanceTo(p);
            F += alpha * (ab*ab + ac*ac) + beta * d2b*d2b;
        }
        return F;
    };

    std::cout << "Сглаживание функционала расстояния-длины: " << nodesCount() << " узлов, " << elementsCount() << " элементов." << std::endl;
    for (short iit = 0; iit < iter_num; iit++)
    {
        if (useFlip/* && iit < 4*/) flip();
        std::cout << iit+1 << " from " << iter_num << "...";
        ConsoleProgress progress(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Node3D node = node_[i];
            AdjacentSet adjasent = node.adjacent;
            Point3D point = node.point;
            double f_current = functor(adjasent, i);

//            for (std::set<UInteger>::iterator it = adjasent.begin(); it != adjasent.end(); ++it)
//            {
//                Triangle t = element_[*it];
//                int index = t.index(i);
//                Point3D B = node_[t[index + 1]].point;
//                Point3D C = node_[t[index + 2]].point;
//                double step = 0.1;
//                int iic = 0;
//                do
//                {
//                    node_[i].point = findBorder(node_[i].point, B, C, func, step, 0.0, level);
//                    double f_dir = functor(adjasent, i);
//                    if (f_dir >= f_current)
//                    {
//                        node_[i].point = point;
//                        f_dir = f_current;
//                        step /= 10.0;
//                    }
//                    else
//                    {
//                        point = node_[i].point;
//                        f_current = f_dir;
//                    }
//                    iic++;
//                } while (step >= 0.001 && iic < 20);
//                step = 0.1;
//                iic = 0;
//                do
//                {
//                    node_[i].point = findBorder(node_[i].point, B, C, func, 0.0, step, level);
//                    double f_dir = functor(adjasent, i);
//                    if (f_dir >= f_current)
//                    {
//                        node_[i].point = point;
//                        f_dir = f_current;
//                        step /= 10.0;
//                    }
//                    else
//                    {
//                        point = node_[i].point;
//                        f_current = f_dir;
//                    }
//                    iic++;
//                } while (step >= 0.001 && iic < 20);
//                step = 0.1;
//                iic = 0;
//                do
//                {
//                    node_[i].point = findBorder(node_[i].point, B, C, func, step, step, level);
//                    double f_dir = functor(adjasent, i);
//                    if (f_dir >= f_current)
//                    {
//                        node_[i].point = point;
//                        f_dir = f_current;
//                        step /= 10.0;
//                    }
//                    else
//                    {
//                        point = node_[i].point;
//                        f_current = f_dir;
//                    }
//                    iic++;
//                } while (step >= 0.001 && iic < 10);

////                Triangle t = element_[*it];
////                int index = t.index(i);
////                Point3D A = node_[t[index]].point;
////                Point3D B = node_[t[index + 1]].point;
////                Point3D C = node_[t[index + 2]].point;
////                double step = 0.01;
////                bool isOptimizedTriangle = false;
////                double f_dir = f_current;
////                double l = step;
////                double f_dir = 0.0;
////                node_[i].point = findBorder(A, B, C, func_slice, l, l, level);
////                f_dir = functor(adjasent);
////                if (f_dir > f_current)
////                {
////                    node_[i].point = point; // revert changes
////                }
////                else
////                {
////                    while (f_dir < f_current && l < 0.3)
////                    {
////                        point = node_[i].point;
////                        f_current = f_dir;
////                        l += step;
////                        node_[i].point = findBorder(A, B, C, func_slice, l, l, level);
////                        f_dir = functor(adjasent);
////                    }
////                    node_[i].point = point;
////                    optimized = isOptimizedTriangle = true;
////                }
////                if (!isOptimizedTriangle)
////                {
////                    l = step;
////                    node_[i].point = findBorder(A, B, C, func_slice, l, 0.0, level);
////                    f_dir = functor(adjasent);
////                    if (f_dir > f_current)
////                    {
////                        node_[i].point = point; // revert changes
////                    }
////                    else
////                    {
////                        while (f_dir < f_current && l < 0.6)
////                        {
////                            point = node_[i].point;
////                            f_current = f_dir;
////                            l += step;
////                            node_[i].point = findBorder(A, B, C, func_slice, l, 0.0, level);
////                            f_dir = functor(adjasent);
////                        }
////                        node_[i].point = point;
////                        optimized = isOptimizedTriangle = true;
////                    }
////                }
////                if (!isOptimizedTriangle)
////                {
////                    l = step;
////                    node_[i].point = findBorder(A, B, C, func_slice, 0, l, level);
////                    f_dir = functor(adjasent);
////                    if (f_dir > f_current)
////                    {
////                        node_[i].point = point; // revert changes
////                    }
////                    else
////                    {
////                        while (f_dir < f_current && l < 0.6)
////                        {
////                            point = node_[i].point;
////                            f_current = f_dir;
////                            l += step;
////                            node_[i].point = findBorder(A, B, C, func_slice, 0.0, l, level);
////                            f_dir = functor(adjasent);
////                        }
////                        node_[i].point = point;
////                        optimized =true;
////                    }
////                }
//            }

            std::set<UInteger> neigbours;
            std::list<Point3D> points;
            for (UInteger elnum: adjasent)
            {
                Triangle t = element_[elnum];
                int index = t.index(i);
                UInteger code1 = t[index + 1];
                UInteger code2 = t[index + 2];
                neigbours.insert(code1);
                neigbours.insert(code2);
                points.push_back((node_[code1].point + node_[code2].point) / 2.0);
            }
            for (UInteger nnode: neigbours)
                points.push_back(node_[nnode].point);
            for (Point3D p: points)
            {
                double step = 0.1;
                double h = (p - node_[i].point).length() * 0.25;
                short local_iter = 0;
                do
                {
                    Point3D A = node_[i].point;
                    Point3D pp = (step * (p - A)) + A;
                    node_[i].point = findBorder(pp, func, h, level);
                    double f_dir = functor(adjasent, i);
                    if (f_dir >= f_current)
                    {
                        node_[i].point = point;
                        f_dir = f_current;
                        step /= 10.0;
                    }
                    else
                    {
                        point = node_[i].point;
                        f_current = f_dir;
                    }
                    ++local_iter;
                } while (step >= 0.001 && local_iter < 10);
            }
            ++progress;
        }
//        flip();
    }
}

void TriangleMesh3D::printStats() const
{
    std::cout << "Triangle mesh (surface): " << nodesCount() << " node(s), " << elementsCount() << " element(s)." << std::endl;
}

void TriangleMesh3D::subdivide(std::list<UInteger> eNumbers, std::function<double (double, double, double)> func)
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
    Node3D local[6];
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
        Node3D n0 = node_[triangle[0]];
        Node3D n1 = node_[triangle[1]];
        Node3D n2 = node_[triangle[2]];
        int code = 0;
        local[0] = n0;
        local[1] = n1;
        local[2] = n2;
        local[3].point = 0.5 * (n0.point + n1.point); local[3].type = BORDER;
        local[4].point = 0.5 * (n1.point + n2.point); local[4].type = BORDER;
        local[5].point = 0.5 * (n2.point + n0.point); local[5].type = BORDER;
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

    if (func != nullptr)
    {
        double h = sqrt((xMax_ - xMin_) * (xMax_ - xMin_) + (yMax_ - yMin_) * (yMax_ - yMin_) + (zMax_ - zMin_) * (zMax_ - zMin_)) / static_cast<double>(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            if (node_[i].type == BORDER)
                node_[i].point = findBorder(node_[i].point, func, h);
        }
    }

    flip();
}

void TriangleMesh3D::transformGrid(const TriangleMesh2D *mesh, std::function<Point3D (Point2D)> func)
{
    clear();

    for (UInteger i = 0; i < mesh->nodesCount(); i++)
        pushNode(func(mesh->point2d(i)), BORDER);

    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        Triangle t = mesh->triangle(i);
        std::swap(t[0], t[2]);
        addElement(t);
    }

    updateDomain();
}

bool TriangleMesh3D::angles(const Point3D &A, const Point3D &B, const Point3D &C, double &alpha, double &beta, double &gamma) const
{
    Point3D AB = B - A;
    Point3D AC = C - A;
    Point3D N = AB.product(AC);
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vz.product(Vx);

    Point3D a = (A - A).inCoordSystem(Vx, Vy, Vz);
    Point3D b = (B - A).inCoordSystem(Vx, Vy, Vz);
    Point3D c = (C - A).inCoordSystem(Vx, Vy, Vz);

    return TriangleMesh2D::angles(a, b, c, alpha, beta, gamma);
}

double TriangleMesh3D::minAngle(const Point3D &A, const Point3D &B, const Point3D &C) const
{
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    angles(A, B, C, alpha, beta, gamma);
//    double alpha = A.angle(C, B);
//    double beta = B.angle(A, C);
//    double gamma = M_PI - alpha - beta;
    return std::min(alpha, std::min(beta, gamma));
}

bool TriangleMesh3D::inCircumSphere(const Point3D &P, const Point3D &A, const Point3D &B, const Point3D &C)
{
    Point3D AB = B - A;
    Point3D AC = C - A;
    Point3D N = AB.product(AC);
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vz.product(Vx);

    Point3D p = (P - A).inCoordSystem(Vx, Vy, Vz);
    Point3D a = (A - A).inCoordSystem(Vx, Vy, Vz);
    Point3D b = (B - A).inCoordSystem(Vx, Vy, Vz);
    Point3D c = (C - A).inCoordSystem(Vx, Vy, Vz);

    double xc = 0.0, yc = 0.0, r = 0.0;

//    return TriangleMesh2D::circumCircle(p.x(), p.y(), a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), xc, yc, r);

    return epsilon_ < r - p.distanceTo(Point3D(xc, yc, 0.0));
}

bool TriangleMesh3D::inCircumCylinder(const Point2D &P, const Point2D &A, const Point2D &B, const Point2D &C, std::function<Point3D (double, double)> domainFunction)
{
    Point3D A3 = domainFunction(A.x(), A.y());
    Point3D B3 = domainFunction(B.x(), B.y());
    Point3D C3 = domainFunction(C.x(), C.y());
    Point3D P3 = domainFunction(P.x(), P.y());

    Point3D AB = B3 - A3;
    Point3D AC = C3 - A3;
    Point3D N = AB.product(AC);
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vz.product(Vx);

    Point3D a = (A3 - A3).inCoordSystem(Vx, Vy, Vz);
    Point3D b = (B3 - A3).inCoordSystem(Vx, Vy, Vz);
    Point3D c = (C3 - A3).inCoordSystem(Vx, Vy, Vz);
    Point3D p = (P3 - A3).inCoordSystem(Vx, Vy, Vz);

    double xc = 0.0, yc = 0.0, r = 0.0;

    TriangleMesh2D::circumCircle(0, 0, a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), xc, yc, r);

    return p.distanceTo(Point3D(xc, yc, 0.0)) < r - epsilon_;
}

Point2D TriangleMesh3D::circumCylinderCenter(const Point2D &A, const Point2D &B, const Point2D &C, std::function<Point3D (double, double)> domainFunction)
{
    Point3D A3 = domainFunction(A.x(), A.y());
    Point3D B3 = domainFunction(B.x(), B.y());
    Point3D C3 = domainFunction(C.x(), C.y());

    Point3D AB = B3 - A3;
    Point3D AC = C3 - A3;
    Point3D N = AB.product(AC);
    Point3D Vx = AB.normalized();
    Point3D Vz = N.normalized();
    Point3D Vy = Vz.product(Vx);

    Point3D a = (A3 - A3).inCoordSystem(Vx, Vy, Vz);
    Point3D b = (B3 - A3).inCoordSystem(Vx, Vy, Vz);
    Point3D c = (C3 - A3).inCoordSystem(Vx, Vy, Vz);
    double xc = 0.0, yc = 0.0, r = 0.0;

    TriangleMesh2D::circumCircle(0.0, 0.0, a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), xc, yc, r);

    Point3D circum(xc, yc, 0.0);
    double u_min = std::min(A.x(), std::min(B.x(), C.x()));
    double v_min = std::min(A.y(), std::min(B.y(), C.y()));
    double u_max = std::max(A.x(), std::max(B.x(), C.x()));
    double v_max = std::max(A.y(), std::max(B.y(), C.y()));
    double du = u_max - u_min;
    double dv = v_max - v_min;
    double u0 = (u_min + u_max) / 2.0;
    double v0 = (v_min + v_max) / 2.0;
    Point3D p0 = domainFunction(u0, v0);
    double d0 = circum.distanceTo((p0 - A3).inCoordSystem(Vx, Vy, Vz));
    u_min -= du;
    u_max += du;
    v_min -= dv;
    v_max += dv;
    du *= 3.0;
    dv *= 3.0;
    for (int i = 0; i < 21; i++)
    {
        double u = u_min + static_cast<double>(i) * du / 20.0;

        for (int j = 0; j < 21; j++)
        {
            double v = v_min + static_cast<double>(j) * dv / 20.0;
            Point3D p = domainFunction(u, v);
            double d = circum.distanceTo((p - A3).inCoordSystem(Vx, Vy, Vz));
            if (d < d0)
            {
                d0 = d;
                p0 = p;
                u0 = u;
                v0 = v;
            }
        }
    }
    return Point2D(u0, v0);
}

bool TriangleMesh3D::insertDelaunayNode(const Point2D &point, const NodeType &type, TriangleMesh2D::Triangulation &triangulation, std::function<Point3D(double, double)> domainFunction)
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
        if ((0.0 <= A.x() && A.x() <= 1.0 && 0.0 <= B.x() && B.x() <= 1.0 && 0.0 <= C.x() && C.x() <= 1.0)
                && (0.0 <= A.y() && A.y() <= 1.0 && 0.0 <= B.y() && B.y() <= 1.0 && 0.0 <= C.y() && C.y() <= 1.0))
        {
            Segment e0(triangle->vertexNode(0), triangle->vertexNode(1));
            Segment e1(triangle->vertexNode(1), triangle->vertexNode(2));
            Segment e2(triangle->vertexNode(2), triangle->vertexNode(0));
            bool flags[] = {false, false, false};
            if (inCircumCylinder(point, A, B, C, domainFunction))
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
                std::cout << '-';
            }
            else
            {
                ++triangle;
            }
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
            if (TriangleMesh2D::signedArea(triangulation.nodes[edges[j][1]], triangulation.nodes[edges[j][0]], triangulation.nodes[number]) > 0.0)
                triangulation.triangles.push_back(Triangle(edges[j][1], edges[j][0], number));
            else
                triangulation.triangles.push_back(Triangle(edges[j][0], edges[j][1], number));
        }
    }
    return true;
}

void TriangleMesh3D::flip()
{
    std::cout << "Flip routine ";
    bool were_flips = true;
    int iic = 0;
    while (were_flips && iic < 40)
    {
        were_flips = false;
        ++iic;
        std::cout << ' ';
        for (UInteger triangle_index = 0; triangle_index < element_.size(); triangle_index++)
        {
            Triangle triangle = element_[triangle_index];
            for (int i = 0; i < 3; i++)
            {
                UInteger index0 = triangle[i];
                UInteger index1 = triangle[i + 1];
                UInteger index2 = triangle[i - 1];
                AdjacentSet a0 = node_[index0].adjacent;
                AdjacentSet a1 = node_[index1].adjacent;
                AdjacentSet a2 = node_[index2].adjacent;
                std::vector<UInteger> common;
                set_intersection(a0.begin(), a0.end(), a1.begin(), a1.end(), std::back_inserter(common));
                if (common.size() > 2)
                {
                    std::cout << "Bad mesh" << std::endl;
                    were_flips = false;
                    return;
                }
                if (common.size() == 2)
                {
                    Point3D p0 = node_[index0].point;
                    Point3D p1 = node_[index1].point;
                    Point3D p2 = node_[index2].point;
                    UInteger index_of_flip_node;
                    int subindex;
                    Point3D f;
                    UInteger flip_triangle_index;
                    AdjacentSet af;
                    if (common[0] != triangle_index)
                    {
                        flip_triangle_index = common[0];
                    }
                    else
                    {
                        flip_triangle_index = common[1];
                    }

                    if (element_[flip_triangle_index][0] != index0 && element_[flip_triangle_index][0] != index1)
                    {
                        index_of_flip_node = element_[flip_triangle_index][0];
                        subindex = 0;
                    }
                    else if (element_[flip_triangle_index][1] != index0 && element_[flip_triangle_index][1] != index1)
                    {
                        index_of_flip_node = element_[flip_triangle_index][1];
                        subindex = 1;
                    }
                    else
                    {
                        index_of_flip_node = element_[flip_triangle_index][2];
                        subindex = 2;
                    }
                    f = node_[index_of_flip_node].point;
//                    node_[index_of_flip_node].adjacent.sort();
                    af = node_[index_of_flip_node].adjacent;
                    std::vector<UInteger> a2_af_common;
                    set_intersection(a2.begin(), a2.end(), af.begin(), af.end(), std::back_inserter(a2_af_common));
                    double tp = 0.0, tq = 0.0;
                    if (a2_af_common.size() == 0 && isSkew(p0, p1, p2, f, tp, tq) && tp > epsilon_ && tp < (1.0 - epsilon_) && tq > epsilon_ && tq < (1.0 - epsilon_))
                    {
//                        double min_c = std::min(minAngle(p0, p1, p2), minAngle(p0, f, p1));
//                        double min_n = std::min(minAngle(p0, f, p2), minAngle(p2, f, p1));
                        if (/*(min_n - min_c) > epsilon_*/inCircumSphere(f, p0, p1, p2))
                        {
//                            std::cout << min_c << " " << min_n << std::endl;
                            std::cout << '.';
                            node_[index0].adjacent.erase(flip_triangle_index);
                            node_[index1].adjacent.erase(triangle_index);
                            element_[triangle_index][i + 1] = index_of_flip_node;
                            node_[index_of_flip_node].adjacent.insert(triangle_index);
                            element_[flip_triangle_index][subindex - 1] = index2;
                            node_[index2].adjacent.insert(flip_triangle_index);
                            were_flips = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    std::cout << iic << std::endl;
}

void TriangleMesh3D::clearCooTriangles(std::list<CooTriangle> &cootriangles, std::function<double(double, double, double)> func, double level, double delta)
{
    if (delta < epsilon_) delta = epsilon_;
    double h = (0.5 * delta > epsilon_) ? 0.5 * delta : 3.0 * epsilon_;
    std::list<CooTriangle>::iterator it = cootriangles.begin();
    while (it != cootriangles.end())
    {
        Point3D a = (*it).a;
        Point3D b = (*it).b;
        Point3D c = (*it).c;
        if ((a.isEqualTo(b, delta) && a.isEqualTo(c, delta)) ||
                (b.isEqualTo(c, delta) && b.isEqualTo(a, delta)) ||
                (c.isEqualTo(a, delta) && c.isEqualTo(b, delta)))
        {
            Point3D mc = findBorder((1.0 / 3.0) * (a + b + c), func, h, level);
            it = cootriangles.erase(it);
            it = cootriangles.begin();
            for (std::list<CooTriangle>::iterator itt = cootriangles.begin(); itt != cootriangles.end(); ++itt)
            {
                if ((*itt).a.isEqualTo(a, delta) || (*itt).a.isEqualTo(b, delta) || (*itt).a.isEqualTo(c, delta))
                    (*itt).a = mc;
                if ((*itt).b.isEqualTo(a, delta) || (*itt).b.isEqualTo(b, delta) || (*itt).b.isEqualTo(c, delta))
                    (*itt).b = mc;
                if ((*itt).c.isEqualTo(a, delta) || (*itt).c.isEqualTo(b, delta) || (*itt).c.isEqualTo(c, delta))
                    (*itt).c = mc;
            }
        }
        else if (a.isEqualTo(b, delta))
        {
            Point3D mc = findBorder(0.5 * (a + b), func, h, level);
            it = cootriangles.erase(it);
            it = cootriangles.begin();
            for (std::list<CooTriangle>::iterator itt = cootriangles.begin(); itt != cootriangles.end(); ++itt)
            {
                if ((*itt).a.isEqualTo(a, delta) || (*itt).a.isEqualTo(b, delta))
                    (*itt).a = mc;
                if ((*itt).b.isEqualTo(a, delta) || (*itt).b.isEqualTo(b, delta))
                    (*itt).b = mc;
                if ((*itt).c.isEqualTo(a, delta) || (*itt).c.isEqualTo(b, delta))
                    (*itt).c = mc;
            }
        }
        else if (a.isEqualTo(c, delta))
        {
            Point3D mc = findBorder(0.5 * (a + c), func, h, level);
            it = cootriangles.erase(it);
            it = cootriangles.begin();
            for (std::list<CooTriangle>::iterator itt = cootriangles.begin(); itt != cootriangles.end(); ++itt)
            {
                if ((*itt).a.isEqualTo(a, delta) || (*itt).a.isEqualTo(c, delta))
                    (*itt).a = mc;
                if ((*itt).b.isEqualTo(a, delta) || (*itt).b.isEqualTo(c, delta))
                    (*itt).b = mc;
                if ((*itt).c.isEqualTo(a, delta) || (*itt).c.isEqualTo(c, delta))
                    (*itt).c = mc;
            }
        }
        else if (c.isEqualTo(b, delta))
        {
            Point3D mc = findBorder(0.5 * (c + b), func, h, level);
            it = cootriangles.erase(it);
            it = cootriangles.begin();
            for (std::list<CooTriangle>::iterator itt = cootriangles.begin(); itt != cootriangles.end(); ++itt)
            {
                if ((*itt).a.isEqualTo(c, delta) || (*itt).a.isEqualTo(b, delta))
                    (*itt).a = mc;
                if ((*itt).b.isEqualTo(c, delta) || (*itt).b.isEqualTo(b, delta))
                    (*itt).b = mc;
                if ((*itt).c.isEqualTo(c, delta) || (*itt).c.isEqualTo(b, delta))
                    (*itt).c = mc;
            }

        }
        else
        {
            ++it;
        }
    }
    ConsoleProgress pb(cootriangles.size());
    it = cootriangles.begin();
    while (it != cootriangles.end())
    {
        Point3D a = (*it).a;
        Point3D b = (*it).b;
        Point3D c = (*it).c;
        bool isFind = false;
        std::list<CooTriangle>::iterator itt = std::next(it);
        while (itt != cootriangles.end())
        {
            Point3D ai = (*itt).a;
            Point3D bi = (*itt).b;
            Point3D ci = (*itt).c;
            if ((a.isEqualTo(ai, delta) || a.isEqualTo(bi, delta) || a.isEqualTo(ci, delta)) &&
                    (b.isEqualTo(ai, delta) || b.isEqualTo(bi, delta) || b.isEqualTo(ci, delta)) &&
                    (c.isEqualTo(ai, delta) || c.isEqualTo(bi, delta) || c.isEqualTo(ci, delta)))
            {
                itt = cootriangles.erase(itt);
                isFind = true;
            }
            else ++itt;
        }

        if (isFind) it = cootriangles.erase(it);
        else ++it;

        ++pb;
    }
}

}
