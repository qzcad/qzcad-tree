#include "quadrilateralmesh2d.h"
#include <iostream>
#include <math.h>

namespace msh
{
QuadrilateralMesh2D::QuadrilateralMesh2D()
{
    xMin_ = -1.0;
    xMax_ = 1.0;
    yMin_ = -1.0;
    yMax_ = 1.0;
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount,
                                         const double &xMin,
                                         const double &yMin, const double &width, const double &height)
{
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

QuadrilateralMesh2D::QuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount,
                                         const Point2D &v0, const Point2D &v1,
                                         const Point2D &v2, const Point2D &v3)
{
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

QuadrilateralMesh2D::QuadrilateralMesh2D(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2)
{
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

QuadrilateralMesh2D::QuadrilateralMesh2D(const UInteger &count, const Point2D &center, const double &radius, unsigned short part)
{
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

QuadrilateralMesh2D::QuadrilateralMesh2D(const QuadrilateralMesh2D &mesh)
{
    element_ = mesh.element_;
    node_ = mesh.node_;
    xMin_ = mesh.xMin_;
    xMax_ = mesh.xMax_;
    yMin_ = mesh.yMin_;
    yMax_ = mesh.yMax_;
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const QuadrilateralMesh2D *mesh)
{
    element_ = mesh->element_;
    node_ = mesh->node_;
    xMin_ = mesh->xMin_;
    xMax_ = mesh->xMax_;
    yMin_ = mesh->yMin_;
    yMax_ = mesh->yMax_;
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
//    UInteger currNumber = 0; // количество внутренних узлов сетки
//    UInteger variablesCount; // количество переменных в функционале
    std::vector<double> x0; // первое приближение
    std::vector<double> x; // минимизация
    std::vector<UInteger> nodeVariable; // номер переменной, соответвствующей узлу
    for (UInteger i = 0; i < node_.size(); i++)
    {
        if (node_[i].type == INNER)
        {
            nodeVariable.push_back(x0.size());
            x0.push_back(node_[i].point.x());
            x0.push_back(node_[i].point.y());
        }
        else
        {
            nodeVariable.push_back(2 * node_.size() + 1);
        }
    }
    x.resize(x0.size());
    std::cout << functional(x0.data(), nodeVariable) << std::endl;
    conjugateGradient(x0.size(), x0.data(), x.data(), nodeVariable);
    std::cout << functional(x.data(), nodeVariable) << std::endl;
    for (UInteger i = 0; i < node_.size(); i++)
    {
        if (node_[i].type == INNER)
        {
            node_[i].point.set(x[nodeVariable[i]], x[nodeVariable[i] + 1]);
        }
    }
}

bool QuadrilateralMesh2D::isBorderElement(const UInteger &number) const
{
    for (int i = 0; i < 4; i++)
    {
        if (node_[element_[number].vertexNode(i)].type == BORDER || node_[element_[number].vertexNode(i)].type == CHARACTER)
            return true;
    }
    return false;
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
    element_.push_back(quad);
    // обновление списка смежных узлов
    node_[node0].adjacent.insert(element_.size() - 1);
    node_[node1].adjacent.insert(element_.size() - 1);
    node_[node2].adjacent.insert(element_.size() - 1);
    node_[node3].adjacent.insert(element_.size() - 1);
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

double QuadrilateralMesh2D::functional(double *vars, const std::vector<UInteger> &nodeVariable)
{
    double f = 0.0;
    for (UInteger i = 0; i < element_.size(); i++)
    {
        double x[4];
        double y[4];
        for (int j = 0; j < 4; j++)
        {
            if(node_[element_[i][j]].type == INNER)
            {
                x[j] = vars[nodeVariable[element_[i][j]]];
                y[j] = vars[nodeVariable[element_[i][j]] + 1];
            }
            else
            {
                // узел граничный
                x[j] = node_[element_[i][j]].point.x();
                y[j] = node_[element_[i][j]].point.y();
            }
        }
        // Рассмотрим 4х-угольник как 4 треугольника, определенных на его углах
        double xc = (x[0] + x[1] + x[2] + x[3]) / 4.0;
        double yc = (y[0] + y[1] + y[2] + y[3]) / 4.0;
        double a [][4] = {{x[0],    x[1],   x[2],   x[3]},
                            {y[0],    y[1],   y[2],   y[3]}};
        double b [][4] = {{x[1],    x[2],   x[3],   x[0]},
                            {y[1],    y[2],   y[3],   y[0]}};
        double c [][4] = {{xc,    xc,   xc,   xc},
                            {yc,    yc,   yc,   yc}};
        double localValue = 0.0;
        // функция для вычисления квадрата числа (C++0x)
        auto sqr = [](double value) { return value * value; };
        for (int q = 0; q < 4; q++)
        {
            double l = sqr(c[0][q] - a[0][q]) + sqr(c[1][q] - a[1][q]) + sqr(b[0][q] - a[0][q]) + sqr(b[1][q] - a[1][q]);
            double o = (c[0][q] - a[0][q]) *  (b[0][q] - a[0][q]) + (c[1][q] - a[1][q]) * (b[1][q] - a[1][q]);
            double alpha = (b[0][q] - a[0][q]) * (c[1][q] - a[1][q]) - (b[1][q] - a[1][q]) * (c[0][q] - a[0][q]);
            //            f += sqr(alpha / 2.0); // площадь
            //            f += 0.3 * sqr(alpha / 2.0) + 0.7 * sqr(o);
            localValue += l;

        }
        f += 0.95 * localValue + 0.05 * (sqr(0.49 - sqr(xc) - sqr(yc)));
    }
    return f;
}

void QuadrilateralMesh2D::nabla(const UInteger &size, double *x, const std::vector<UInteger> &nodeVariable, double *gradient, double h)
{
    double x_plus_h; // x + h
    double x_minus_h; // x - h
    double x_plus_h_h; // x + 2 * h
    double x_minus_h_h; // x - 2 * h
    double f_x_plus_h; // f(x + h)
    double f_x_minus_h; // f(x - h)
    double f_x_plus_h_h; // f(x + 2 * h)
    double f_x_minus_h_h; // f(x - 2 * h)
    double x_i;
    for (UInteger i = 0; i < size; i++)
    {
        x_i = x[i];
        x_plus_h = x_i + h;
        x_minus_h = x_i - h;
        x_plus_h_h = x_plus_h + h;
        x_minus_h_h = x_minus_h - h;
        x[i] = x_plus_h; f_x_plus_h = functional(x, nodeVariable);
        x[i] = x_minus_h; f_x_minus_h = functional(x, nodeVariable);
        x[i] = x_plus_h_h; f_x_plus_h_h = functional(x, nodeVariable);
        x[i] = x_minus_h_h; f_x_minus_h_h = functional(x, nodeVariable);
        gradient[i] = (f_x_minus_h_h - 8.0 * f_x_minus_h + 8.0 * f_x_plus_h - f_x_plus_h_h) / (12.0 * h);
        x[i] = x_i;
    }
}

double QuadrilateralMesh2D::lambda(const UInteger &size, double *x, double *s, const double &lambda_val, const std::vector<UInteger> &nodeVariable)
{
    double x_lambda[size];
    for (UInteger i = 0; i < size; i++)
    {
        x_lambda[i] = x[i] + lambda_val * s[i];
    }
    return functional(x_lambda, nodeVariable);
}

double QuadrilateralMesh2D::goldenRatio(const UInteger &size, const double &a, const double &b, double *x0, double *s, const std::vector<UInteger> &nodeVariable, double epsilon, UInteger maxIter)
{
    double left = a;
    double right = b;
    const double phi = (sqrt(5.0) + 1.0) / 2.0; // golden ratio
    double x1 = b - (b - a) / phi;
    double x2 = a + (b - a) / phi;
    double y1 = lambda(size, x0, s, x1, nodeVariable);
    double y2 = lambda(size, x0, s, x2, nodeVariable);

    for(UInteger count = 0; count < maxIter; count++)
    {
        if(y1 <= y2)
        {
            right = x2;
            x2 = x1;
            x1 = right - (right - left) / phi;
            y2 = y1;
            y1 = lambda(size, x0, s, x1, nodeVariable);
        }
        else
        {
            left = x1;
            x1 = x2;
            x2 = left + (right - left) / phi;
            y1 = y2;
            y2 = lambda(size, x0, s, x2, nodeVariable);
        }
        if((right - left) < epsilon) break;
    }
    if(y1 <= y2) return x1;
    return x2;
}

double QuadrilateralMesh2D::norm2(const UInteger &size, double *x)
{
    double norm = 0.0;
    for (UInteger i = 0; i < size; i++)
    {
        norm += x[i] * x[i];
    }
    return norm;
}

void QuadrilateralMesh2D::conjugateGradient(const UInteger &size, double *x0, double *xMin, const std::vector<UInteger> &nodeVariable, double epsilon, UInteger maxIter)
{
    double xk[size];
    double xkj[size];
    double xkj1[size];
    double dxkj[size];
    double nablaxkj[size];
    double nablaxkj1[size];
    double Skj[size];
    double lambda;
    double omega;
    double nkj, nkj1;

    for(UInteger i = 0; i < size; i++) xk[i] = x0[i];

    //        if(dialog != NULL)
    //        {
    //            dialog->SetPosition(0);
    //            dialog->SetMaxPosition(maxSearchIterations);
    //            dialog->SetMessage("Оптимизация. Минимизация функционала методом сопряженных градиентов");
    //        }

    for(UInteger k = 0; k < maxIter; k++)
    {
        //        if(dialog != NULL)
        //        {
        //            dialog->IncPosition();
        //        }

        for(UInteger i = 0; i < size; i++)
            xkj[i] = xk[i];
        for(UInteger j = 0; j < maxIter; j++)
        {
            nabla(size, xkj, nodeVariable, Skj);
            lambda = goldenRatio(size, -1.0, 1.0, xkj, Skj, nodeVariable, epsilon, maxIter);
            for(UInteger i = 0; i < size; i++) xkj1[i] = xkj[i] + lambda * Skj[i];
            nabla(size, xkj, nodeVariable, nablaxkj);
            nabla(size, xkj1, nodeVariable, nablaxkj1);
            nkj = norm2(size, xkj);
            nkj1 = norm2(size, xkj1);
            omega = (nkj1 * nkj1) / (nkj * nkj);

            for(UInteger i = 0; i < size; i++)
                Skj[i] = -nablaxkj[i] + omega * Skj[i];

            for(UInteger i = 0; i < size; i++)
                dxkj[i] = xkj1[i] - xkj[i];

            if(norm2(size, dxkj) < epsilon)    // norm (varCount, Skj) < epsilon ||
            {
                for(UInteger i = 0; i < size; i++) xMin[i] = xkj1[i];
                return;
            }
            for(UInteger i = 0; i < size; i++) xkj[i] = xkj1[i];
        }
        for(UInteger i = 0; i < size; i++) xk[i] = xkj1[i];
    }
    for(UInteger i = 0; i < size; i++) xMin[i] = xk[i];
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
