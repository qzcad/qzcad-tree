#include "mesh.h"
#include <iostream>

namespace msh {

double Mesh::epsilon_ = 1.0E-10; // инициализация точности

Mesh::Mesh(const Mesh *mesh)
{
    if (mesh != NULL)
    {
        layer_ = mesh->layer_;
        data_ = mesh->data_;
    }
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

}


