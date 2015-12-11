#include "nameddoublevector.h"

namespace msh {

NamedDoubleVector::NamedDoubleVector(const std::string &name, const std::vector<double> &vector)
{
    name_ = name;
    vector_ = vector;
}

std::vector<double> NamedDoubleVector::vector() const
{
    return vector_;
}

void NamedDoubleVector::setVector(const std::vector<double> &vector)
{
    vector_ = vector;
}

std::string NamedDoubleVector::name() const
{
    return name_;
}

void NamedDoubleVector::setName(const std::string &name)
{
    name_ = name;
}

}
