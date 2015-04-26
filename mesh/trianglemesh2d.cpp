#include "trianglemesh2d.h"
#include <iostream>
#include <math.h>
#include <map>
#include <climits>
#include <float.h>

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
            addElement(i * yCount + j, (i + 1) * yCount + j, i * yCount + j + 1);
            addElement((i + 1) * yCount + j, (i + 1) * yCount + j + 1, i * yCount + j + 1);
        }
    }
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
//    minimizeFunctional();
    std::cout << "Создана равномерная сетка треугольных элементов: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
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


double TriangleMesh2D::minAngle(const Point2D &A, const Point2D &B, const Point2D &C)
{
    double alpha, beta, gamma;
    angles(A, B, C, alpha, beta, gamma);
    return std::min(alpha, std::min(beta, gamma));
}

bool TriangleMesh2D::angles(const Point2D &A, const Point2D &B, const Point2D &C, double &alpha, double &beta, double &gamma)
{
    const double epsilon = 1.0E-12;
    const double a = B.distanceTo(C); // сторона, противолежащяя вершине A (BC)
    const double b = A.distanceTo(C); // сторона, противолежащяя вершине B (AC)
    const double c = A.distanceTo(B); // сторона, противолежащяя вершине C (AB)
    if (a < epsilon || b < epsilon || c < epsilon)
    {
        alpha = beta = gamma = 0.0;
        return false;
    }
    // Теорема косинусов
    alpha = acos((b*b + c*c - a*a) / (2.0 * b * c)); // Угол в вершине A
    // Теорема синусов
    beta = asin(sin(alpha) * b / a); // Угол в вершине B
    // Теорема о сумме углов треугольника
    gamma = M_PI - (alpha + beta); // Угол в вершине C
    return true;
}

TriangleMesh2D::TriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint)
{
    const double epsilon = 1.0E-6;
    xMin_ = xMin;
    xMax_ = xMin + width;
    yMin_ = yMin;
    yMax_ = yMin + height;
    const double hx = width / (double)(xCount - 1);
    const double hy = height / (double)(yCount - 1);
    const double minDistance = 0.25 * sqrt(hx*hx + hy*hy);
    std::map<UInteger, UInteger> nodesMap;
    // формирование массива узлов
    for (UInteger i = 0; i < xCount; i++)
    {
        double x = xMin + (double) i * hx;
        for (UInteger j = 0; j < yCount; j++)
        {
            double y = yMin + (double) j * hy;
            Point2D point(x, y);

            if (func(x, y) > 0.0)
            {
                nodesMap[i * yCount + j] = pushNode(point, INNER);
            }
        }
    }
    // формирование начальной сетки
    const double xCenter = (xMax_ + xMin_) / 2.0;
    const double yCenter = (yMax_ + yMin_) / 2.0;
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
                    addElement(iter0->second, iter1->second, iter3->second);
                    addElement(iter1->second, iter2->second, iter3->second);
                }
                else
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
    for (UInteger i = 0; i < normal.size(); i++)
    {
        iso[i] = ULONG_MAX;
        if (normal[i].length() > epsilon)
        {
            Point2D current = node_[i].point;
            // двоичный поиск граничной точки
            Point2D inner = current;
            Point2D outer = current + sqrt(hx*hx + hy*hy) * normal[i];
            Point2D mid;
            do
            {
                mid = 0.5 * (inner + outer);
                double val = func(mid.x(), mid.y());
                if (val < 0.0)
                    outer = mid;
                else if (val > 0.0)
                    inner = mid;
                else
                    break;
            } while(inner.distanceTo(outer) > epsilon);
            if (current.distanceTo(mid) < minDistance)
            {
                node_[i].point = mid;
                node_[i].type = BORDER;
                iso[i] = i;
            }
            else
            {
                // сравнение с существующими изо-точками перед всатвкой
                bool isExist = false;
                for (UInteger j = normal.size(); j < nodesCount(); j++)
                {
                    Point2D border = node_[j].point;
                    if (border.distanceTo(mid) < minDistance)
                    {
                        iso[i] = j;
                        isExist = true;
                        break;
                    }
                }
                if (!isExist)
                {
                    // поиск соответствующей характерной точки
                    for (std::list<Point2D>::iterator cPoint = charPoint.begin(); cPoint != charPoint.end(); ++cPoint)
                    {
                        if (mid.distanceTo(*cPoint) < minDistance)
                        {
                            mid = *cPoint;
                            charPoint.erase(cPoint);
                            break;
                        }
                    }
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
            double d = (n->point).distanceTo(*cPoint);
            if (d < min_d)
            {
                min_d = d;
                min_n = n;
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
                        if (triangle[j + 1] != iso[triangle[j + 1]])
                            addElement(triangle[j + 1], triangle[j], iso[triangle[j + 1]]);
                        if (triangle[j] != iso[triangle[j]] && iso[triangle[j]] != iso[triangle[j + 1]])
                            addElement(triangle[j], iso[triangle[j]], iso[triangle[j + 1]]);
                    }
                    else
                    {
                        if (triangle[j] != iso[triangle[j]])
                            addElement(triangle[j + 1], triangle[j], iso[triangle[j]]);
                        if (triangle[j + 1] != iso[triangle[j + 1]] && iso[triangle[j]] != iso[triangle[j + 1]])
                            addElement(triangle[j + 1], iso[triangle[j]], iso[triangle[j + 1]]);
                    }
                }
            }
        }
    }

    // the end.
    std::cout << "Создана сетка треугольных элементов для функционального объекта: узлов - " << nodesCount() << ", элементов - " << elementsCount() << "." << std::endl;
}

}
