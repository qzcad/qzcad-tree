#ifndef MESHLIST_H
#define MESHLIST_H

#include <functional>
#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>

#include "point3d.h"
#include "nodetype.h"

namespace msh
{

class ElementList;

class NodeList3d
{
public:
    Point3D point;
    NodeType type;
    std::set<ElementList *> incident;

    NodeList3d(const Point3D &p, const NodeType &t)
    {
        point = p;
        type = t;
    }
};

class ElementList
{
public:
    std::vector<NodeList3d *> vertices;

    ElementList(int vcount=3)
    {
        vertices.resize(vcount);
    }

    ElementList(NodeList3d *n0, NodeList3d *n1, NodeList3d *n2, NodeList3d *n3)
    {
        vertices.resize(4);
        vertices[0] = n0;
        vertices[1] = n1;
        vertices[2] = n2;
        vertices[3] = n3;
    }

    ElementList(const std::vector<NodeList3d *> &v)
    {
        vertices = v;
    }
};

typedef std::list<NodeList3d *> NodesList;
typedef std::list<ElementList *> ElementsList;

class MeshList
{
public:
    MeshList();

    NodeList3d *addNode(const Point3D &p, const NodeType &t, double epsilon=1.0E-7)
    {
        for (NodesList::reverse_iterator i = nodes.rbegin(); i != nodes.rend(); i++)
        {
            NodeList3d *nptr = (*i);
            if (p.isEqualTo(nptr->point, epsilon))
            {
                return nptr;
            }
        }
        NodeList3d *n = new NodeList3d(p, t);
        nodes.push_back(n);
        return n;
    }

    void addElement(NodeList3d *n0, NodeList3d *n1, NodeList3d *n2, NodeList3d *n3)
    {
        ElementList *e = new ElementList(n0, n1, n2, n3);
        n0->incident.insert(e);
        n1->incident.insert(e);
        n2->incident.insert(e);
        n3->incident.insert(e);
        elements.push_back(e);
    }

    void addElement(const std::vector<NodeList3d *> v)
    {
        ElementList *e = new ElementList(v);
        for (std::vector<NodeList3d *>::iterator i = e->vertices.begin(); i != e->vertices.end(); i++)
        {
            (*i)->incident.insert(e);
        }
        elements.push_back(e);
    }

    void delNode(NodeList3d *n)
    {
        if (n->incident.empty())
        {
//            nodes.remove(n);
            NodeList3d *nul = NULL;
            std::replace(nodes.begin(), nodes.end(), n, nul);
            delete n;
        }
        else
        {
            std::set<ElementList *> Incident = n->incident;
            for (std::set<ElementList *>::iterator i = Incident.begin(); i != Incident.end(); i++)
            {
                ElementList *e = (*i);
                delElement(e);
            }
        }
    }

//    std::vector<NodeList3d *>::iterator delNode(std::vector<NodeList3d *>::iterator it)
//    {
////        if ((*it)->incident.empty())
////        {
////            delete (*it);
////            return nodes.erase(it)
////        }
//    }

    void delElement(ElementList *e)
    {
        for (std::vector<NodeList3d *>::iterator i = e->vertices.begin(); i != e->vertices.end(); i++)
        {
            NodeList3d *n = *i;
            n->incident.erase(e);
            if (n->incident.empty())
                delNode(n);
        }
        elements.remove(e);
        delete e;
    }

    void clearFuncNodes(std::function<double(double, double, double)> func, long maxCount)
    {
        long c = nodes.size();
        for (NodesList::iterator i = nodes.begin(); i != nodes.end() && c - nodes.size() < maxCount;)
        {
            NodeList3d *n = (*i);
            if (n != NULL)
            {
                Point3D p = n->point;
                if (func(p.x(), p.y(), p.z()) < 0.0)
                {
                    delNode(n);
//                    std::cout << n;
//                    if (*i == NULL) std::cout << "-";
//                    else  std::cout << "!";
                    i = nodes.erase(i);
                }
                else {
                    ++i;
                }
            }
            else
            {
                i = nodes.erase(i);
//                std::cout << ".";
            }
        }
    }

    ~MeshList()
    {
        for (NodesList::iterator i = nodes.begin(); i != nodes.end(); i++)
        {
            delete (*i);
        }
        for (ElementsList::iterator i = elements.begin(); i != elements.end(); i++)
        {
            delete (*i);
        }
    }
public:
    NodesList nodes;
    ElementsList elements;
};

}


#endif // MESHLIST_H
