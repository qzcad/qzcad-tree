#include "mesh.h"

namespace msh {

double Mesh::epsilon_ = 1.0E-10; // инициализация точности

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

}


