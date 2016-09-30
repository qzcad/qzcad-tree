#include "trianglemesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#include <iostream>
#include <algorithm>
#include "funcopt.h"
#include "rfunctions.h"

#include "trianglemesh2d.h"

namespace msh {

TriangleMesh3D::TriangleMesh3D() : Mesh3D(NULL)
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
    node_ = mesh.node_;
    element_ = mesh.element_;
}

TriangleMesh3D::TriangleMesh3D(const TriangleMesh3D *mesh) : Mesh3D(mesh)
{
    node_ = mesh->node_;
    element_ = mesh->element_;
}

void TriangleMesh3D::cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length)
{
    clear();
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
}

void TriangleMesh3D::cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, std::function<double (double, double, double)> func)
{
    clear();
    const double xi_max = (double)lCount * 2.0 * M_PI / (double)rCount;
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
    double dl = xi_max / (double)lCount;
    double dphi = 2.0 * M_PI / (double)rCount;
    // двоичный поиск вдоль шва
    for (UInteger i = 0; i < lCount; i++)
    {
        double xi0 = (double)i * dl;
        double xi1 = (double)(i + 1) * dl;
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
        double eta0 = (double)i * dphi;
        double eta1 = (double)(i + 1) * dphi;
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
}

void TriangleMesh3D::coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length)
{
    clear();
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
}

void TriangleMesh3D::coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, std::function<double (double, double, double)> func)
{
    clear();
    const double xi_max = (double)lCount * 2.0 * M_PI / (double)rCount;
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
    double dxi = xi_max / (double)lCount;
    double dphi = 2.0 * M_PI / (double)rCount;
    // двоичный поиск вдоль шва
    for (UInteger i = 0; i < lCount; i++)
    {
        double xi0 = (double)i * dxi;
        double xi1 = (double)(i + 1) * dxi;
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
        double eta0 = (double)i * dphi;
        double eta1 = (double)(i + 1) * dphi;
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
}

void TriangleMesh3D::parametricDomain(const UInteger &uCount, const UInteger &vCount, std::function<Point3D (double, double)> domainFunction, std::function<double (double, double, double)> rfunc)
{
    clear();

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
    double du = 1.0 / (double)uCount;
    double dv = 1.0 / (double)vCount;
    // двоичный поиск по первому направлению
    for (UInteger i = 0; i < uCount; i++)
    {
        double u0 = (double)i * du;
        double u1 = (double)(i + 1) * du;
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
        double v0 = (double)i * dv;
        double v1 = (double)(i + 1) * dv;
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
    smesh.functionalDomain(uCount, vCount, 0.0, 0.0, 1.0, 1.0, func2d, charPoints2d, 0.0, false, distance);
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

void TriangleMesh3D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0], nodes_ref[1], nodes_ref[2]);
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

double TriangleMesh3D::minAngle(const UInteger &elNum)
{
    const Triangle tri = element_[elNum];
    const Point3D p0 = node_[tri[0]].point;
    const Point3D p1 = node_[tri[1]].point;
    const Point3D p2 = node_[tri[2]].point;
//    if (minAngle(p0, p1, p2) < 0.1) std::cout << elNum << " ";
    return minAngle(p0, p1, p2);
}

void TriangleMesh3D::clearElements()
{
    element_.clear();
}

bool TriangleMesh3D::angles(const Point3D &A, const Point3D &B, const Point3D &C, double &alpha, double &beta, double &gamma)
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

double TriangleMesh3D::minAngle(const Point3D &A, const Point3D &B, const Point3D &C)
{
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    angles(A, B, C, alpha, beta, gamma);
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

    TriangleMesh2D::circumCircle(0, 0, a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), xc, yc, r);

    return p.distanceTo(Point3D(xc, yc, 0.0)) < r - epsilon_;
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
        double u = u_min + (double)i * du / 20.0;

        for (int j = 0; j < 21; j++)
        {
            double v = v_min + (double)j * dv / 20.0;
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
    bool were_flips = true;
    while (were_flips)
    {
        were_flips = false;
        for (UInteger t = 0; t < element_.size(); t++)
        {
            Triangle triangle = element_[t];
            for (int i = 0; i < 3; i++)
            {
                UInteger index0 = triangle[i];
                UInteger index1 = triangle[i + 1];
                AdjacentSet a1 = node_[index0].adjacent;
                AdjacentSet a2 = node_[index1].adjacent;
                std::vector<UInteger> common;
                set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), std::back_inserter(common));
                if (common.size() > 2)std::cout << ">";
                if (common.size() == 2)
                {
                    //                std::cout << common[0] << " " << common[1] << std::endl;
                    Point3D p0 = node_[index0].point;
                    Point3D p1 = node_[index1].point;
                    UInteger index2 = triangle[i - 1];
                    Point3D p2 = node_[index2].point;
                    UInteger indexf;
                    int subindex;
                    Point3D f;
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
                    double tp = 0.0, tq = 0.0;
//                    if (min_n > min_c && isSkew(p0, p1, p2, f, tp, tq) && tp > epsilon_ && tp < (1.0 - epsilon_) && tq > epsilon_ && tq < (1.0 - epsilon_))
//                    {
//                        std::cout << " flip " << std::endl;
//                        node_[index0].adjacent.erase(t1);
//                        node_[index1].adjacent.erase(t);
//                        element_[t][i + 1] = indexf;
//                        node_[indexf].adjacent.insert(t);
//                        element_[t1][subindex - 1] = index2;
//                        node_[index2].adjacent.insert(t1);
//                        were_flips = true;
//                        break;
//                    }
                    if (isSkew(p0, p1, p2, f, tp, tq) && tp > epsilon_ && tp < (1.0 - epsilon_) && tq > epsilon_ && tq < (1.0 - epsilon_))
                    {
                        if (inCircumSphere(f, p0, p1, p2))
                        {
                            std::cout << " flip " << std::endl;
                            node_[index0].adjacent.erase(t1);
                            node_[index1].adjacent.erase(t);
                            element_[t][i + 1] = indexf;
                            node_[indexf].adjacent.insert(t);
                            element_[t1][subindex - 1] = index2;
                            node_[index2].adjacent.insert(t1);
                            were_flips = true;
                            break;
                        }
                    }
                }
            }
        }
    }
}

}
