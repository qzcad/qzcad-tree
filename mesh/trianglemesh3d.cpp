#include "trianglemesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <map>
#include <climits>
#include <iostream>

#include "trianglemesh2d.h"

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
    double length_2 = length / 2.0;
    auto con = [](const double &x, const double &y)
    {
        return x + y - sqrt(x*x + y*y);
    };

    auto func2d = [&](double xi, double eta)
    {
        double x = radius * cos(eta);
        double y = xi;
        double z = radius * sin(eta);
        double f = func(x, y, z);
        double r = con(length_2 * length_2 - (xi - length_2) * (xi - length_2), M_PI * M_PI - (eta - M_PI) * (eta - M_PI));
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
    if (func2d(length, 0.0) < epsilon_ )
    {
        charPoints2d.push_back(Point2D(length, 0.0));
        charPoints2d.push_back(Point2D(length, 2.0 * M_PI));
    }
    double dl = length / (double)lCount;
    double dphi = 2.0 * M_PI / (double)rCount;
    // двоичный поиск вдоль шва
    for (UInteger i = 0; i < lCount; i++)
    {
        double xi0 = (double)i * dl;
        double xi1 = (double)(i + 1) * dl;
        double eta = 0.0;
        double x = radius * cos(eta);
        double y = radius * sin(eta);
        double val0 = func(x, xi0, y);
        double val1 = func(x, xi1, y);
        if ((val0 > epsilon_ && val1 < epsilon_) || (val0 < epsilon_ && val1 > epsilon_))
        {
            double c, val;
            do
            {
                c = 0.5 * (xi0 + xi1);
                val = func(x, c, y);
                if (val0 * val < 0.0)
                {
                    xi1 = c;
                    val1 = val;
                }
                else
                {
                    xi0 = c;
                    val0 = val;
                }
            } while (fabs(xi0 - xi1) > epsilon_ && fabs(val) > epsilon_);
            charPoints2d.push_back(Point2D(c, 0.0));
            charPoints2d.push_back(Point2D(c, 2.0 * M_PI));
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
        double val0 = func(radius * cos(eta0), xi0, radius * sin(eta0));
        double val1 = func(radius * cos(eta1), xi0, radius * sin(eta1));
        if ((val0 > epsilon_ && val1 < epsilon_) || (val0 < epsilon_ && val1 > epsilon_))
        {
            double c, val;
            do
            {
                c = 0.5 * (eta0 + eta1);
                val = func(radius * cos(c), xi0, radius * cos(c));
                if (val0 * val < 0.0)
                {
                    eta1 = c;
                    val1 = val;
                }
                else
                {
                    eta0 = c;
                    val0 = val;
                }
            } while (fabs(eta0 - eta1) > epsilon_ && fabs(val) > epsilon_);
            charPoints2d.push_back(Point2D(xi0, c));
            std::cout << "...b...";
        }
        if (fabs(val0) < epsilon_ && i > 0)
        {
            charPoints2d.push_back(Point2D(xi0, eta0));
            std::cout << "...0...";
        }
        if (fabs(val1) < epsilon_ && i < rCount - 1)
        {
            charPoints2d.push_back(Point2D(xi0, eta1));
            std::cout << "...1...";
        }
    }
    mesh2d.ruppert(lCount, rCount, -0.001, -0.001, length + 0.002, 2.0 * M_PI + 0.002, func2d, charPoints2d, true);

    for (UInteger i = 0; i < mesh2d.nodesCount(); i++)
    {
        Point2D p = mesh2d.point2d(i);
        NodeType nodeType = (fabs(p.x()) < epsilon_ || fabs(p.x() - length) < epsilon()) ? CHARACTER : BORDER;
        if (fabs(p.y() - 2.0 * M_PI) < epsilon_)
            nodes_map[i] = addNode(Point3D(radius * cos(p.y()), p.x(), radius * sin(p.y())), CHARACTER);
        else
            nodes_map[i] = pushNode( Point3D(radius * cos(p.y()), p.x(), radius * sin(p.y())), nodeType);
    }
    for (UInteger i = 0; i < mesh2d.elementsCount(); i++)
    {
        Triangle t = mesh2d.triangle(i);
        addElement(nodes_map[t[1]], nodes_map[t[0]], nodes_map[t[2]]);
//        addElement(t[1], t[0], t[2]);
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

double TriangleMesh3D::minAngle(const UInteger &elNum)
{
    const Triangle tri = element_[elNum];
    const Point3D p0 = node_[tri[0]].point;
    const Point3D p1 = node_[tri[1]].point;
    const Point3D p2 = node_[tri[2]].point;
    return minAngle(p0, p1, p2);
}

bool TriangleMesh3D::angles(const Point3D &A, const Point3D &B, const Point3D &C, double &alpha, double &beta, double &gamma)
{
    const double a = B.distanceTo(C); // сторона, противолежащяя вершине A (BC)
    const double b = A.distanceTo(C); // сторона, противолежащяя вершине B (AC)
    const double c = A.distanceTo(B); // сторона, противолежащяя вершине C (AB)
    if (a < epsilon_ || b < epsilon_ || c < epsilon_)
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

double TriangleMesh3D::minAngle(const Point3D &A, const Point3D &B, const Point3D &C)
{
    double alpha = 0.0, beta = 0.0, gamma = 0.0;
    angles(A, B, C, alpha, beta, gamma);
    return std::min(alpha, std::min(beta, gamma));
}

}
