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
                                         const Floating &xMin,
                                         const Floating &yMin, const Floating &width, const Floating &height)
{
    Floating hx = width / (Floating)(xCount - 1);
    Floating hy = height / (Floating)(yCount - 1);
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        Floating x = xMin + (Floating) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            Floating y = yMin + (Floating) j * hy;
            Point2D point(x, y);

            if ((i == 0 && j == 0) || (i == 0 & j == yCount - 1) || (i == xCount - 1 && j == 0) || (i == xCount - 1 & j == yCount - 1))
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
    minimizeFunctional();
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount,
                                         const Point2D &v0, const Point2D &v1,
                                         const Point2D &v2, const Point2D &v3)
{
    Floating hx = 2.0 / (Floating)(xCount - 1);
    Floating hy = 2.0 / (Floating)(yCount - 1);
    for (UInteger i = 0; i < xCount; i++)
    {
        Floating xi = -1.0 + (Floating) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            Floating eta = -1.0 + (Floating) j * hy;
            Point2D point = isoFunc(0, xi, eta) * v0  + isoFunc(1, xi, eta) * v1 + isoFunc(2, xi, eta) * v2 + isoFunc(3, xi, eta) * v3;

            if ((i == 0 && j == 0) || (i == 0 & j == yCount - 1) || (i == xCount - 1 && j == 0) || (i == xCount - 1 & j == yCount - 1))
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
}

QuadrilateralMesh2D::QuadrilateralMesh2D(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2)
{
    Point2D center = (v0 + v1 + v2) / 3.0; // центр треугольника
    Point2D c01 = (v0 + v1) / 2.0; // центр стороны, соединяющей вершину 0 и 1
    Point2D c12 = (v1 + v2) / 2.0; // центр стороны, соединяющей вершину 1 и 2
    Point2D c20 = (v2 + v0) / 2.0; // центр стороны, соединяющей вершину 2 и 0
    UInteger sideCount = count / 2 + 1;
    Floating h = 2.0 / (Floating)(sideCount - 1); // шаг изо-сетки
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
            Floating xi = -1.0 + (Floating) i * h;
            for (UInteger j = 0; j < sideCount; j++)
            {
                Floating eta = -1.0 + (Floating) j * h;
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
    std::vector<Floating> x0; // первое приближение
    std::vector<Floating> x; // минимизация
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

Floating QuadrilateralMesh2D::nodeValue(const UInteger &number) const
{
    return node_[number].type;
}

bool QuadrilateralMesh2D::isBorderElement(const UInteger &number) const
{
    if (node_[number].type == BORDER || node_[number].type == CHARACTER)
        return true;
    return false;
}

Quadrilateral QuadrilateralMesh2D::quadrilateral(const UInteger &number) const
{
    return element_[number];
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

Floating QuadrilateralMesh2D::isoFunc(const UInteger &i, const Floating &xi, const Floating &eta)
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

Floating QuadrilateralMesh2D::functional(Floating *vars, const std::vector<UInteger> &nodeVariable)
{
    Floating f = 0.0;
    for (UInteger i = 0; i < element_.size(); i++)
    {
        Floating x[4];
        Floating y[4];
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
        Floating xc = (x[0] + x[1] + x[2] + x[3]) / 4.0;
        Floating yc = (y[0] + y[1] + y[2] + y[3]) / 4.0;
        Floating a [][4] = {{x[0],    x[1],   x[2],   x[3]},
                            {y[0],    y[1],   y[2],   y[3]}};
        Floating b [][4] = {{x[1],    x[2],   x[3],   x[0]},
                            {y[1],    y[2],   y[3],   y[0]}};
        Floating c [][4] = {{xc,    xc,   xc,   xc},
                            {yc,    yc,   yc,   yc}};
        Floating localValue = 0.0;
        auto sqr = [](Floating v) { return v * v; };
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

void QuadrilateralMesh2D::nabla(const UInteger &size, Floating *x, const std::vector<UInteger> &nodeVariable, Floating *gradient, Floating h)
{
    Floating x_plus_h; // x + h
    Floating x_minus_h; // x - h
    Floating x_plus_h_h; // x + 2 * h
    Floating x_minus_h_h; // x - 2 * h
    Floating f_x_plus_h; // f(x + h)
    Floating f_x_minus_h; // f(x - h)
    Floating f_x_plus_h_h; // f(x + 2 * h)
    Floating f_x_minus_h_h; // f(x - 2 * h)
    Floating x_i;
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

Floating QuadrilateralMesh2D::lambda(const UInteger &size, Floating *x, Floating *s, const Floating &lambda_val, const std::vector<UInteger> &nodeVariable)
{
    Floating x_lambda[size];
    for (UInteger i = 0; i < size; i++)
    {
        x_lambda[i] = x[i] + lambda_val * s[i];
    }
    return functional(x_lambda, nodeVariable);
}

Floating QuadrilateralMesh2D::goldenRatio(const UInteger &size, const Floating &a, const Floating &b, Floating *x0, Floating *s, const std::vector<UInteger> &nodeVariable, Floating epsilon, UInteger maxIter)
{
    Floating left = a;
    Floating right = b;
    const Floating phi = (sqrt(5.0) + 1.0) / 2.0; // golden ratio
    Floating x1 = b - (b - a) / phi;
    Floating x2 = a + (b - a) / phi;
    Floating y1 = lambda(size, x0, s, x1, nodeVariable);
    Floating y2 = lambda(size, x0, s, x2, nodeVariable);

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

Floating QuadrilateralMesh2D::norm2(const UInteger &size, Floating *x)
{
    Floating norm = 0.0;
    for (UInteger i = 0; i < size; i++)
    {
        norm += x[i] * x[i];
    }
    return norm;
}

void QuadrilateralMesh2D::conjugateGradient(const UInteger &size, Floating *x0, Floating *xMin, const std::vector<UInteger> &nodeVariable, Floating epsilon, UInteger maxIter)
{
    Floating xk[size];
    Floating xkj[size];
    Floating xkj1[size];
    Floating dxkj[size];
    Floating nablaxkj[size];
    Floating nablaxkj1[size];
    Floating Skj[size];
    Floating lambda;
    Floating omega;
    Floating nkj, nkj1;

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
}
