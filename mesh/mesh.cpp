#include "mesh.h"
#include <iostream>
#include <algorithm>
#include <limits>

namespace msh {

double Mesh::epsilon_ = static_cast<double>(std::numeric_limits<double>::epsilon()) * 1000.0; // инициализация точности

Mesh::Mesh(const Mesh *mesh)
{
    if (mesh != NULL)
    {
        layer_ = mesh->layer_;
        data_ = mesh->data_;
    }
}

bool Mesh::isBorderElement(ElementPointer element) const
{
    for (int i = 0; i < element->verticesCount(); i++)
    {
        NodeType t = nodeType(element->vertexNode(i));
        if (t == BORDER || t == CHARACTER) return true;
    }
    return false;
}

bool Mesh::isBorderFace(const UIntegerVector &face) const
{
    for (UIntegerVector::size_type k = 0; k < face.size(); k++)
        if (nodeType(face[k]) == INNER)
            return false;
    return true;
}

int Mesh::layer(const msh::UInteger &number) const
{
    return layer_[number];
}

void Mesh::setLayer(const UInteger &number, const int &l)
{
    layer_[number] = l;
}

void Mesh::pushLayer(const int &l)
{
    layer_.push_back(l);
}

void Mesh::clearLayers()
{
    layer_.clear();
}

UInteger Mesh::sizeOfLayers() const
{
    return layer_.size();
}
double Mesh::epsilon()
{
    return epsilon_;
}

void Mesh::setEpsilon(double epsilon)
{
    epsilon_ = epsilon;
}

UInteger Mesh::dataVectorsCount() const
{
    return data_.size();
}

void Mesh::addDataVector(const NamedDoubleVector &vec)
{
    data_.push_back(vec);
}

void Mesh::addDataVector(const std::string &name, const std::vector<double> &values)
{
    NamedDoubleVector vec(name, values);
    data_.push_back(vec);
}

const NamedDoubleVector &Mesh::data(const UInteger &i) const
{
    return data_[i];
}

void Mesh::clearDataVectors()
{
    data_.clear();
}

void Mesh::clear()
{
    clearDataVectors();
    clearLayers();
    clearNodes();
    clearElements();
}

void Mesh::printDataExtremums()
{
    for (UInteger i = 0; i < data_.size(); i++)
    {
        NamedDoubleVector ndv = data_[i];
        std::cout << ndv.min() << "\t<=\t" << ndv.name() << "\t<=\t" << ndv.max() << std::endl;
    }
}

int Mesh::facePower(const UIntegerVector &face) const
{
    AdjacentSet s0 = adjacent(face[0]);
    AdjacentSet s1 = adjacent(face[1]);
    AdjacentSet intersect;
    std::set_intersection(s0.begin(), s0.end(), s1.begin(), s1.end(), std::inserter(intersect, intersect.begin()));
    for (UIntegerVector::size_type i = 2; i < face.size(); i++)
    {
        AdjacentSet current = adjacent(face[i]);
        AdjacentSet local_intersect;
        std::set_intersection(intersect.begin(), intersect.end(), current.begin(), current.end(), std::inserter(local_intersect, local_intersect.begin()));
        intersect = local_intersect;
    }
    return intersect.size();
}

}


