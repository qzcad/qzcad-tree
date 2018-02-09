#include "quadrilateralmesh3d.h"
#undef __STRICT_ANSI__
#include <math.h>
#include <iostream>

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
    double hphi = 2.0 * M_PI / (double)rCount;
    double hl = length / (double)lCount;
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
    double hxi = 1.0 / (double)rCount;
    double heta = 1.0 / (double)lCount;
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
    std::set<UInteger> border;
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
        else if (code != (1<<m) - 1)
        {
            for (int j = 0; j < m; j++)
            {
                if (values[j] - level >= epsilon_)
                    border.insert(el->vertexNode(j));
            }
        }
        ++progress;
    }
    // delete all isolated elements
    for (std::list<ElementPointer>::iterator it = inner.begin(); it != inner.end();)
    {
        short inner_faces_count = 0;
        ElementPointer el = *it;
        for (std::list<ElementPointer>::iterator iit = inner.begin(); iit != inner.end(); iit++)
        {
            ElementPointer ell = *iit;
            if (it != iit)
            {
                short common_nodes_count = 0;
                for (int j = 0; j < 8; j++)
                {
                    if (ell->in(el->vertexNode(j)))
                        ++common_nodes_count;
                }
                if (common_nodes_count == 4)
                {
                    ++inner_faces_count;
                }
                if (inner_faces_count == 3)
                    break;
            }
        }
        if (inner_faces_count < 3)
        {
            it = inner.erase(it);
            it = inner.begin();
            std::cout << '.' << std::ends;
        }
        else
        {
            it++;
        }
    }
    // loop a list of background elements and add border faces into the mesh
    for (ElementPointer el: inner)
    {
        int m = el->facesCount();
        for (int j = 0; j < m; j++)
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
            if (!isInner/* && border.find(f[0]) != border.end() && border.find(f[1]) != border.end() && border.find(f[2]) != border.end() && border.find(f[3]) != border.end()*/)
                addElement(addNode(p[0], BORDER), addNode(p[1], BORDER), addNode(p[2], BORDER), addNode(p[3], BORDER));
            if (d < h)
                h = d;
        }
    }
    // loop nodes of the mesh and find correct position
    h *= 0.25;
    const double dh = 0.5 * h;
    UInteger corrected_count = node_.size();
    while (corrected_count > 0)
    {
        corrected_count = 0;
        for (UInteger i = 0; i != node_.size(); ++i)
        {
            Point3D point = node_[i].point;
            double f = func(point.x(), point.y(), point.z()) - level;
            if (fabs(f) >= epsilon_)
            {
                Point3D g = grad(func, point, dh).normalized();
                Point3D p = (f < 0.0) ? (point + h * g) : (point - h * g);
                double fp = func(p.x(), p.y(), p.z()) - level;
                if (signbit(f) != signbit(fp))
                {
                    p = binary(point, p, func, level);
                }
                node_[i].point = p;
                corrected_count++;
            }
            else
            {
                AdjacentSet adjacent = node_[i].adjacent;
                Point3D center(0.0, 0.0, 0.0);
                AdjacentSet neighbours;
                for (UInteger elnum: adjacent)
                {
                    Quadrilateral q = element_[elnum];
                    int index = q.index(i);
                    neighbours.insert(q[index + 1]);
                    neighbours.insert(q[index + 2]);
                    neighbours.insert(q[index + 3]);
                }
                for (UInteger npointer: neighbours)
                {
                    center = center + node_[npointer].point;
                }
                center.scale(1.0 / (double)neighbours.size());
                node_[i].point = findBorder(center, func, h, level);
            }
        }
        std::cout << corrected_count << " " << std::ends;
    }
    /*
    progress.restart(node_.size());
    std::vector<Point3D> surface(node_.size());
    for (UInteger i = 0; i != node_.size(); ++i)
    {
//        AdjacentSet adjacent = node_[i].adjacent;
        Point3D point = node_[i].point;
//        Point3D n(0.0, 0.0, 0.0);
//        double avr_len = 0.0;
//        for (UInteger elnum: adjacent)
//        {
//            Quadrilateral q = element_[elnum];
//            int index = q.index(i);
//            Point3D prev = node_[q[index - 1]].point;
//            Point3D next = node_[q[index + 1]].point;
//            Point3D a(prev, point);
//            Point3D b(point, next);
//            avr_len += 0.5 * (a.length() + b.length());
//            n = n + normal3(point, prev, next);
//        }
//        avr_len /= (double)adjacent.size();
//        surface[i] = findBorder(point, n.normalized(), func, sqrt(3.0) * avr_len, level);
        surface[i] = findBorder(point, func, 0.25 * h, level);
        ++progress;
    }
    for (UInteger i = 0; i != node_.size(); ++i) node_[i].point = surface[i];*/
    laplacianSmoothing(func, level, smooth);
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
            point.scale(1.0 / (double)neighbours.size());
            for (UInteger elnum: adjacent)
            {
                Quadrilateral q = element_[elnum];
                int index = q.index(nnode);
                Point3D a(point, node_[q[index - 1]].point);
                Point3D b(point, node_[q[index + 1]].point);
                avr_len += 0.5 * (a.length() + b.length());
            }
            avr_len /= (double)adjacent.size();
            node_[nnode].point = findBorder(point, func, 0.1 * avr_len, level);
            ++progress;
        }
    }
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
            Point3D center = 0.25 * (A + B + C + D);
            Point3D p = findBorder(center, func, (B - A).length() * 0.25, level);
            double d2b = center.distanceTo(p);
//            double d2b = distToBorder(A, B, C, func, 0.33333, 0.33333, level);
            Point3D AB(A, B);
            Point3D AC(A, C);
            double ab = AB.length();
            double ac = AC.length();
            F += alpha * 0.5 * (ab*ab + ac*ac) + beta * d2b*d2b;
        }
        return F;
    };

    std::cout << "Сглаживание функционала расстояния-длины: " << nodesCount() << " узлов, " << elementsCount() << " элементов." << std::endl;
    bool optimized = true;
    for (short iit = 0; iit < iter_num && optimized; iit++)
    {
        ConsoleProgress progress(nodesCount());
        for (UInteger i = 0; i < nodesCount(); i++)
        {
            Node3D node = node_[i];
            AdjacentSet adjasent = node.adjacent;
            Point3D point = node.point;
            double f_current = functor(adjasent, i);

            for (std::set<UInteger>::iterator it = adjasent.begin(); it != adjasent.end(); ++it)
            {
                Quadrilateral q = element_[*it];
                int index = q.index(i);
                Point3D B = node_[q[index + 1]].point;
                Point3D C = node_[q[index - 1]].point;
                double step = 0.1;
                int iic = 0;
                do
                {
                    node_[i].point = findBorder(node_[i].point, B, C, func, step, 0.0, level);
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
                    iic++;
                } while (step >= 0.001 && iic < 10);
                /*if (iic == 3)*/step = 0.1;/*else step=0.001;*/
                iic = 0;
                do
                {
                    node_[i].point = findBorder(node_[i].point, B, C, func, 0.0, step, level);
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
                    iic++;
                } while (step >= 0.001 && iic < 10);
                iic = 0;
                do
                {
                    node_[i].point = findBorder(node_[i].point, B, C, func, step, step, level);
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
                    iic++;
                } while (step >= 0.001 && iic < 5);
            }
            ++progress;
        }
    }
}
}
