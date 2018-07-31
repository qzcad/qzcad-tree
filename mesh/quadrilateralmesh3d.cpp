#include "quadrilateralmesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <iostream>
#include <algorithm>
#include <set>

#include "consoleprogress.h"

namespace msh {

QuadrilateralMesh3D::QuadrilateralMesh3D() : Mesh3D(NULL)
{
}

QuadrilateralMesh3D::QuadrilateralMesh3D(const QuadrilateralMesh3D &mesh) : Mesh3D(&mesh)
{
    element_ = mesh.element_;
}

QuadrilateralMesh3D::QuadrilateralMesh3D(const QuadrilateralMesh3D *mesh) : Mesh3D(mesh)
{
    element_ = mesh->element_;
}

void QuadrilateralMesh3D::cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length)
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
            addElement(p0, p1, p2, p3);
        }
    }
    // размеры области
    xMin_ = zMin_ = -radius;
    xMax_ = zMax_ = radius;
    yMin_ = 0.0;
    yMax_ = length;
}

void QuadrilateralMesh3D::cylinderDomain(const UInteger &rCount, const UInteger &lCount, std::function<double (double)> radius, const double &length)
{
    clear();
    double hphi = 2.0 * M_PI / static_cast<double>(rCount);
    double hl = length / static_cast<double>(lCount);
    double phi = 0.0;
    double r_max = 0.0;
    // формирование массива узлов
    for (UInteger i = 0; i < rCount; i++)
    {
        double l = 0.0;
        for (UInteger j = 0; j <= lCount; j++)
        {
            double r = radius(l);
            Point3D point(r * cos(phi), l, r * sin(phi));
            if (j == 0 || j == lCount)
                pushNode(point, CHARACTER);
            else
                pushNode(point, BORDER);
            l += hl;
            if (r > r_max) r_max = r;
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
            addElement(p0, p1, p2, p3);
        }
    }
    // размеры области
    xMin_ = zMin_ = -r_max;
    xMax_ = zMax_ = r_max;
    yMin_ = 0.0;
    yMax_ = length;
}

void QuadrilateralMesh3D::coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length)
{
    clear();
    double hxi = 1.0 / static_cast<double>(rCount);
    double heta = 1.0 / static_cast<double>(lCount);
    double xi = 0.0;
    double radius;
    double phi;
    double p = (bottom_radius < top_radius) ? bottom_radius / top_radius : 2.0 - top_radius / bottom_radius;
    double q = 2.0;
    // формирование массива узлов
    for (UInteger i = 0; i < rCount; i++)
    {
        double eta = 0.0;

        for (UInteger j = 0; j <= lCount; j++)
        {
            double s = p * eta + (1.0 - p) * (1.0 - tanh(q * (1.0 - eta)) / tanh(q)); // формула растяжения Эйземана, Флетчер, Вычислительные методы в динамике жидкостей, том 2, стр. 123
            if (j == 0) s = 0.0;
            else if (j == lCount) s = 1.0;
            radius = bottom_radius + s * (top_radius - bottom_radius);
            phi = 2.0 * M_PI * xi;
            Point3D point(radius * cos(phi), s * length, radius * sin(phi));
            if (j == 0 || j == lCount)
                pushNode(point, CHARACTER);
            else
                pushNode(point, BORDER);
            eta += heta;
        }
        xi += hxi;
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
            addElement(p0, p1, p2, p3);
        }
    }
    // размеры области
    xMin_ = zMin_ = -(bottom_radius > top_radius ? bottom_radius : top_radius);
    xMax_ = zMax_ = (bottom_radius > top_radius ? bottom_radius : top_radius);
    yMin_ = 0.0;
    yMax_ = length;
}

void QuadrilateralMesh3D::backgroundGrid(const HexahedralMesh3D *mesh, std::function<double (double, double, double)> func, double level, int smooth, int optimize)
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
            if (value - level <= -epsilon_)
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
            Point3D p[4];
            p[0] = mesh->point3d(f[0]);
            p[1] = mesh->point3d(f[1]);
            p[2] = mesh->point3d(f[2]);
            p[3] = mesh->point3d(f[3]);
            double d = p[0].distanceTo(p[1]);
            bool isInner = false;
            // check an other element with these vertices
            for (ElementPointer el1: inner)
            {
                if (el1 != el && el1->in(f[0]) && el1->in(f[1]) && el1->in(f[2]) && el1->in(f[3]))
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
                UInteger index3 = addNode(p[3], BORDER);
                addElement(index0, index1, index2, index3);
                AdjacentSet a0 = node_[index0].adjacent;
                AdjacentSet a1 = node_[index1].adjacent;
                AdjacentSet a2 = node_[index2].adjacent;
                AdjacentSet a3 = node_[index3].adjacent;
                std::vector<UInteger> common01;
                std::vector<UInteger> common12;
                std::vector<UInteger> common23;
                std::vector<UInteger> common30;
                set_intersection(a0.begin(), a0.end(), a1.begin(), a1.end(), std::back_inserter(common01));
                set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), std::back_inserter(common12));
                set_intersection(a2.begin(), a2.end(), a3.begin(), a3.end(), std::back_inserter(common23));
                set_intersection(a3.begin(), a3.end(), a0.begin(), a0.end(), std::back_inserter(common30));
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
                if (common23.size() > 2)
                {
                    isWeak = true;
                    weak0 = f[2];
                    weak1 = f[3];
                    break;
                }
                if (common30.size() > 2)
                {
                    isWeak = true;
                    weak0 = f[3];
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
                if (ell->in(weak0) && ell->in(weak1))
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
            Quadrilateral q = element_[elnum];
            int index = q.index(i);
            Point3D prev = node_[q[index - 1]].point;
            Point3D next = node_[q[index + 1]].point;
            n = n + normal3(point, prev, next);
        }
        normals[i] = n.normalized();
    }
    for (short j = 0; j < 4; j++)
    {
        for (UInteger i = 0; i != node_.size(); ++i)
        {
            AdjacentSet adjacent = node_[i].adjacent;
            Point3D n = normals[i];
            AdjacentSet neighbours;
            for (UInteger elnum: adjacent)
            {
                Quadrilateral q = element_[elnum];
                int index = q.index(i);
                neighbours.insert(q[index + 1]);
                neighbours.insert(q[index - 1]);
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
    laplacianSmoothing(func, level, smooth);
    distlenSmoothing(func, level, optimize);
    refineTopology(func, level);
    distlenSmoothing(func, level, optimize);
    std::cout << "Surface quadrilateral mesh: nodes - " << nodesCount() << ", elements - " << elementsCount() << std::endl;
    updateDomain();
}

UInteger QuadrilateralMesh3D::elementsCount() const
{
    return element_.size();
}

ElementPointer QuadrilateralMesh3D::element(const UInteger &number) const
{
    ElementPointer elementPtr = &element_[number];
    return elementPtr;
}

double QuadrilateralMesh3D::surfaceArea() const
{
    double s = 0.0;
    for (UInteger i = 0; i < elementsCount(); i++)
        s += area(i);
    return s;
}

void QuadrilateralMesh3D::addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3)
{
    Quadrilateral quad(node0, node1, node2, node3);
    addElement(quad);
}

void QuadrilateralMesh3D::addElement(const Quadrilateral &quad)
{
    element_.push_back(quad);
    // обновление списка смежных узлов
    node_[quad[0]].adjacent.insert(element_.size() - 1);
    node_[quad[1]].adjacent.insert(element_.size() - 1);
    node_[quad[2]].adjacent.insert(element_.size() - 1);
    node_[quad[3]].adjacent.insert(element_.size() - 1);
}

void QuadrilateralMesh3D::addElement(const std::vector<UInteger> &nodes_ref)
{
    addElement(nodes_ref[0], nodes_ref[1], nodes_ref[2], nodes_ref[3]);
}

double QuadrilateralMesh3D::area(const UInteger &number) const
{
    Quadrilateral quad = element_[number];
    Point3D p0 = node_[quad[0]].point;
    Point3D p1 = node_[quad[1]].point;
    Point3D p2 = node_[quad[2]].point;
    Point3D p3 = node_[quad[3]].point;
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

void QuadrilateralMesh3D::add(const QuadrilateralMesh3D *mesh)
{
    std::vector<UInteger> nodesPointers(mesh->nodesCount());
    for (UInteger i = 0; i < mesh->nodesCount(); i++)
    {
        nodesPointers[i] = addNode(mesh->node_[i]);
    }
    for (UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        Quadrilateral q = mesh->element_[i];
        addElement(nodesPointers[q[0]], nodesPointers[q[1]], nodesPointers[q[2]], nodesPointers[q[3]]);
    }
    updateDomain();
}

void QuadrilateralMesh3D::translate(const double &x, const double &y, const double &z)
{
    for (UInteger i = 0; i < nodesCount(); i++)
    {
        Point3D p = node_[i].point;
        node_[i].point.set(p.x() + x, p.y() + y, p.z() + z);
    }
}

void QuadrilateralMesh3D::clearElements()
{
    element_.clear();
}

void QuadrilateralMesh3D::laplacianSmoothing(std::function<double (double, double, double)> func, double level, int iter_num)
{
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
                Quadrilateral q = element_[elnum];
                int index = q.index(nnode);
                neighbours.insert(q[index + 1]);
                neighbours.insert(q[index + 2]);
                neighbours.insert(q[index + 3]);
            }
            for (UInteger npointer: neighbours)
            {
                point = point + node_[npointer].point;
            }
            point.scale(1.0 /  static_cast<double>(neighbours.size()));
            for (UInteger elnum: adjacent)
            {
                Quadrilateral q = element_[elnum];
                int index = q.index(nnode);
                Point3D a(point, node_[q[index - 1]].point);
                Point3D b(point, node_[q[index + 1]].point);
                avr_len += 0.5 * (a.length() + b.length());
            }
            avr_len /=  static_cast<double>(adjacent.size());
            node_[nnode].point = findBorder(point, func, 0.1 * avr_len, level);
            ++progress;
        }
//        if ((iit + 1) % 4 == 0)
            flip();
    }
//    refineTopology(func, level);
}

void QuadrilateralMesh3D::distlenSmoothing(std::function<double (double, double, double)> func, double level, int iter_num)
{
    auto functor = [&](const AdjacentSet &adjasentset, const UInteger &nnode)
    {
        // alpha = 0.001 is optimal for gradient searching
        const double alpha = 0.001;
        const double beta = 1.0 - alpha;
        double F = 0.0; // функционал
        for (UInteger adj: adjasentset)
        {
            Quadrilateral q = element_[adj];
            int index = q.index(nnode);
            Point3D A = node_[q[index]].point;
            Point3D B = node_[q[index + 1]].point;
            Point3D C = node_[q[index - 1]].point;
            Point3D D = node_[q[index + 2]].point;
            double ab = A.distanceTo(B);
            double ac = A.distanceTo(C);
            double ad = A.distanceTo(D);
            Point3D center = (A + B + C + D) / 4.0;
            Point3D p = findBorder(center, func, min(ab, ac) * 0.25, level);
            double d2b = center.distanceTo(p);
            F += alpha * (ab*ab + ac*ac + ad*ad) + beta * d2b*d2b;

//            Point3D N = AB.product(AC);
//            Point3D Vx = AB.normalized();
//            Point3D Vz = N.normalized();
//            Point3D Vy = Vz.product(Vx);
//            Point3D a = (A - A).inCoordSystem(Vx, Vy, Vz);
//            Point3D b = (B - A).inCoordSystem(Vx, Vy, Vz);
//            Point3D c = (C - A).inCoordSystem(Vx, Vy, Vz);
////            Point3D d = (D - A).inCoordSystem(Vx, Vy, Vz);
//            double area = (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
//            Point3D center = 0.25 * (A + B + C + D);
//            Point3D p = findBorder(center, func, min(ab, ac) * 0.25, level);
//            double d2b = center.distanceTo(p);
//            F += alpha * exp(10.0 * area) + beta * d2b*d2b;
        }
//        for (UInteger adj: adjasentset)
//        {
//            Quadrilateral q = element_[adj];
//            for (int i = 0; i < 4; i++)
//            {
//                Point3D A = node_[q[i]].point;
//                Point3D B = node_[q[i + 1]].point;
//                Point3D C = node_[q[i - 1]].point;
//                Point3D AB(A, B);
//                Point3D AC(A, C);
//                double ab = AB.length();
//                double ac = AC.length();

////                Point3D N = AB.product(AC);
////                Point3D Vx = AB.normalized();
////                Point3D Vz = N.normalized();
////                Point3D Vy = Vz.product(Vx);
////                Point3D a = (A - A).inCoordSystem(Vx, Vy, Vz);
////                Point3D b = (B - A).inCoordSystem(Vx, Vy, Vz);
////                Point3D c = (C - A).inCoordSystem(Vx, Vy, Vz);
////                double area = (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
//                Point3D center = 0.3333333333333333 * (A + B + C);
//                Point3D p = findBorder(center, func, min(ab, ac) * 0.25, level);
//                double d2b = center.distanceTo(p);
////                F += alpha * exp(10.0 * area) + beta * d2b*d2b;

//                F += alpha * 0.5 * (ab*ab + ac*ac) + beta * d2b*d2b;
//            }
//        }
        return F;
    };

    std::cout << "Сглаживание функционала расстояния-длины: " << nodesCount() << " узлов, " << elementsCount() << " элементов." << std::endl;
    bool optimized = true;
    for (short iit = 0; iit < iter_num && optimized; iit++)
    {
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
//                Quadrilateral q = element_[*it];
//                int index = q.index(i);
//                Point3D B = node_[q[index + 1]].point;
//                Point3D C = node_[q[index - 1]].point;
////                Point3D D = node_[q[index + 2]].point;
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
//            }
            std::set<UInteger> neigbours;
            std::list<Point3D> points;
            for (UInteger elnum: adjasent)
            {
                Quadrilateral quad = element_[elnum];
                int index = quad.index(i);
                neigbours.insert(quad[index + 1]);
                neigbours.insert(quad[index + 2]);
                neigbours.insert(quad[index + 3]);
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
                    Point3D pp = step * (p - A) + A;
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
        if ((iit + 1) % 4 == 0)
            flip();
    }
    //refineTopology(func, level);
}

double QuadrilateralMesh3D::minAngle(const UInteger &elNum) const
{
    Quadrilateral quad = element_[elNum];
    Point3D A = node_[quad[0]].point;
    Point3D B = node_[quad[1]].point;
    Point3D C = node_[quad[2]].point;
    Point3D D = node_[quad[3]].point;
    double a = A.angle(D, B);
    double b = B.angle(A, C);
    double c = C.angle(B, D);
    double d = D.angle(C, A);
    return std::max(std::max(a, b), std::max(c, d));
}

double QuadrilateralMesh3D::maxAngle(const UInteger &elNum) const
{
    Quadrilateral quad = element_[elNum];
    Point3D A = node_[quad[0]].point;
    Point3D B = node_[quad[1]].point;
    Point3D C = node_[quad[2]].point;
    Point3D D = node_[quad[3]].point;
    return maxAngle(A, B, C, D);
}

double QuadrilateralMesh3D::maxAngle(const Point3D &A, const Point3D &B, const Point3D &C, const Point3D &D) const
{
    double a = A.angle(D, B);
    double b = B.angle(A, C);
    double c = C.angle(B, D);
    double d = D.angle(C, A);
    return std::max(std::max(a, b), std::max(c, d));
}

void QuadrilateralMesh3D::flip()
{
    std::cout << "Flip routine ";
    bool were_flips = true;
    int iic = 0;
    while (were_flips && iic < 40)
    {
        were_flips = false;
        ++iic;
        std::cout << ' ';
        for (UInteger quad_index = 0; quad_index < element_.size(); quad_index++)
        {
            Quadrilateral quad = element_[quad_index];
            for (int i = 0; i < 4; i++)
            {
                UInteger index0 = quad[i + 0];
                UInteger index1 = quad[i + 1];
                UInteger index2 = quad[i + 2];
                UInteger index3 = quad[i + 3];
                AdjacentSet a0 = node_[index0].adjacent;
                AdjacentSet a1 = node_[index1].adjacent;
//                AdjacentSet a2 = node_[index2].adjacent;
//                AdjacentSet a3 = node_[index3].adjacent;
                std::vector<UInteger> common;
                set_intersection(a0.begin(), a0.end(), a1.begin(), a1.end(), std::back_inserter(common));
                if (common.size() > 2)
                {
//                    std::cout << "Bad mesh" << std::endl;
                    were_flips = false;
//                    return;
                }
                if (common.size() == 2)
                {
                    Point3D p0 = node_[index0].point;
                    Point3D p1 = node_[index1].point;
                    Point3D p2 = node_[index2].point;
                    Point3D p3 = node_[index3].point;
                    UInteger index_of_flip_node2;
                    UInteger index_of_flip_node3;
                    int subindex2 = -1, subindex3 = -1;
                    Point3D f2, f3;
                    UInteger flip_quad_index;
//                    AdjacentSet af2, af3;
                    if (common[0] != quad_index)
                    {
                        flip_quad_index = common[0];
                    }
                    else
                    {
                        flip_quad_index = common[1];
                    }
                    if (element_[flip_quad_index][1] == index0 && element_[flip_quad_index][0] == index1)
                    {
                        subindex2 = 2;
                        subindex3 = 3;
                        index_of_flip_node2 = element_[flip_quad_index][subindex2];
                        index_of_flip_node3 = element_[flip_quad_index][subindex3];
                    }
                    else if (element_[flip_quad_index][2] == index0 && element_[flip_quad_index][1] == index1)
                    {
                        subindex2 = 3;
                        subindex3 = 0;
                        index_of_flip_node2 = element_[flip_quad_index][subindex2];
                        index_of_flip_node3 = element_[flip_quad_index][subindex3];
                    }
                    else if (element_[flip_quad_index][3] == index0 && element_[flip_quad_index][2] == index1)
                    {
                        subindex2 = 0;
                        subindex3 = 1;
                        index_of_flip_node2 = element_[flip_quad_index][subindex2];
                        index_of_flip_node3 = element_[flip_quad_index][subindex3];
                    }
                    else if (element_[flip_quad_index][0] == index0 && element_[flip_quad_index][3] == index1)
                    {
                        subindex2 = 1;
                        subindex3 = 2;
                        index_of_flip_node2 = element_[flip_quad_index][subindex2];
                        index_of_flip_node3 = element_[flip_quad_index][subindex3];
                    }
                    f2 = node_[index_of_flip_node2].point;
                    f3 = node_[index_of_flip_node3].point;
                    double current_angle = std::max(maxAngle(p0, p1, p2, p3), maxAngle(f2, f3, p0, p1));
                    double angle2 = std::max(maxAngle(f3, p1, p2, p3), maxAngle(f2, f3, p3, p0));
                    double angle3 = std::max(maxAngle(p0, f2, p2, p3), maxAngle(f2, f3, p1, p2));
//                    node_[index_of_flip_node].adjacent.sort();
//                    af = node_[index_of_flip_node].adjacent;
//                    std::vector<UInteger> a2_af_common;
//                    set_intersection(a2.begin(), a2.end(), af.begin(), af.end(), std::back_inserter(a2_af_common));
                    double tp0 = 0.0, tq0 = 0.0;
                    double tp1 = 0.0, tq1 = 0.0;
                    double tp2 = 0.0, tq2 = 0.0;
                    const double delta = 0.1;
                    if (angle3 < current_angle && angle3 < angle2 && (isSkew(p0, p1, p2, f2, tp0, tq0) && tp0 > delta && tp0 < (1.0 - delta) && tq0 > delta && tq0 < (1.0 - delta)) &&
                            (isSkew(p0, p1, p3, f2, tp1, tq1) && tp1 > delta && tp1 < (1.0 - delta) && tq1 > delta && tq1 < (1.0 - delta)) &&
                            (isSkew(p0, p1, p2, f3, tp2, tq2) && tp2 > delta && tp2 < (1.0 - delta) && tq2 > delta && tq2 < (1.0 - delta)))
                    {
                        std::cout << '.';
                        node_[index0].adjacent.erase(flip_quad_index);
                        node_[index1].adjacent.erase(quad_index);
                        element_[quad_index][i + 1] = index_of_flip_node2;
                        node_[index_of_flip_node2].adjacent.insert(quad_index);
                        element_[flip_quad_index][subindex2 - 1] = index2;
                        node_[index2].adjacent.insert(flip_quad_index);
                        were_flips = true;
                        break;
                    }
                    if (angle2 < current_angle && angle2 < angle3 && (isSkew(p0, p1, p3, f3, tp0, tq0) && tp0 > delta && tp0 < (1.0 - delta) && tq0 > delta && tq0 < (1.0 - delta)) &&
                            (isSkew(p0, p1, p2, f3, tp1, tq1) && tp1 > delta && tp1 < (1.0 - delta) && tq1 > delta && tq1 < (1.0 - delta)) &&
                            (isSkew(p0, p1, p3, f2, tp2, tq2) && tp2 > delta && tp2 < (1.0 - delta) && tq2 > delta && tq2 < (1.0 - delta)))
                    {
                        std::cout << '.';
                        node_[index1].adjacent.erase(flip_quad_index);
                        node_[index0].adjacent.erase(quad_index);
                        element_[quad_index][i + 0] = index_of_flip_node3;
                        node_[index_of_flip_node3].adjacent.insert(quad_index);
                        element_[flip_quad_index][subindex3 + 1] = index3;
                        node_[index3].adjacent.insert(flip_quad_index);
                        were_flips = true;
                        break;
                    }
                }
            }
        }
    }
    std::cout << iic << std::endl;
}

void QuadrilateralMesh3D::refineTopology(std::function<double(double, double, double)> func, double level)
{
    const double alpha = 2.44;
    std::list<std::vector<Point3D>> emesh;
    for (Quadrilateral quad: element_)
    {
        Point3D p0 = point3d(quad[0]);
        Point3D p1 = point3d(quad[1]);
        Point3D p2 = point3d(quad[2]);
        Point3D p3 = point3d(quad[3]);
        if (maxAngle(p0, p1, p2, p3) >= alpha)
        {
            if ((p0.angle(p3, p1) >= alpha && p2.angle(p1, p3) >= alpha) || (p1.angle(p0, p2) >= alpha && p3.angle(p2, p0) >= alpha))
            {
                std::vector<Point3D> t(3);
                Point3D c = 0.25 * (p0 + p1 + p2 + p3);
                t[0] = p0; t[1] = p1; t[2] = c;
                emesh.push_back(t);
                t[0] = p1; t[1] = p2; t[2] = c;
                emesh.push_back(t);
                t[0] = p2; t[1] = p3; t[2] = c;
                emesh.push_back(t);
                t[0] = p3; t[1] = p0; t[2] = c;
                emesh.push_back(t);
            }
            else if (p1.angle(p0, p2) >= alpha && p3.angle(p2, p0) >= alpha)
            {
                std::vector<Point3D> t0(3), t1(3);
                t0[0] = p0; t0[1] = p1; t0[2] = p3;
                t1[0] = p1; t1[1] = p2; t1[2] = p3;
                emesh.push_back(t0);
                emesh.push_back(t1);
            }
            else
            {
                std::vector<Point3D> t0(3), t1(3);
                t0[0] = p0; t0[1] = p1; t0[2] = p2;
                t1[0] = p0; t1[1] = p2; t1[2] = p3;
                emesh.push_back(t0);
                emesh.push_back(t1);
            }
        }
        else
        {
            std::vector<Point3D> t(4);
            t[0] = p0; t[1] = p1; t[2] = p2; t[3] = p3;
            emesh.push_back(t);
        }
    }
    clear();
    std::cout<<emesh.size();
    for (std::vector<Point3D> pol: emesh)
    {
        if (pol.size() == 3)
        {
            Point3D p0 = pol[0];
            Point3D p1 = pol[1];
            Point3D p2 = pol[2];
            Point3D c = 1.0 / 3.0 * (p0 + p1 + p2);
//            double h = c.distanceTo(p0);
            Point3D c01 = 0.5 * (p0 + p1);
            Point3D c12 = 0.5 * (p1 + p2);
            Point3D c20 = 0.5 * (p2 + p0);
//            Point3D c01 = findBorder(0.5 * (p0 + p1), func, h, level);
//            Point3D c12 = findBorder(0.5 * (p1 + p2), func, h, level);
//            Point3D c20 = findBorder(0.5 * (p2 + p0), func, h, level);
//            c = findBorder(c, func, h, level);
            addElement(addNode(p0, BORDER), addNode(c01, BORDER), addNode(c, BORDER), addNode(c20, BORDER));
            addElement(addNode(p1, BORDER), addNode(c12, BORDER), addNode(c, BORDER), addNode(c01, BORDER));
            addElement(addNode(p2, BORDER), addNode(c20, BORDER), addNode(c, BORDER), addNode(c12, BORDER));
        }
        else
        {
            Point3D p0 = pol[0];
            Point3D p1 = pol[1];
            Point3D p2 = pol[2];
            Point3D p3 = pol[3];
            Point3D c = 0.25 * (p0 + p1 + p2 + p3);
            Point3D c01 = 0.5 * (p0 + p1);
            Point3D c12 = 0.5 * (p1 + p2);
            Point3D c23 = 0.5 * (p2 + p3);
            Point3D c30 = 0.5 * (p3 + p0);
//            double h = c.distanceTo(p0);
//            Point3D c01 = findBorder(0.5 * (p0 + p1), func, h, level);
//            Point3D c12 = findBorder(0.5 * (p1 + p2), func, h, level);
//            Point3D c23 = findBorder(0.5 * (p2 + p3), func, h, level);
//            Point3D c30 = findBorder(0.5 * (p3 + p0), func, h, level);
//            c = findBorder(c, func, h, level);
            addElement(addNode(p0, BORDER), addNode(c01, BORDER), addNode(c, BORDER), addNode(c30, BORDER));
            addElement(addNode(p1, BORDER), addNode(c12, BORDER), addNode(c, BORDER), addNode(c01, BORDER));
            addElement(addNode(p2, BORDER), addNode(c23, BORDER), addNode(c, BORDER), addNode(c12, BORDER));
            addElement(addNode(p3, BORDER), addNode(c30, BORDER), addNode(c, BORDER), addNode(c23, BORDER));
        }
    }
}
}
