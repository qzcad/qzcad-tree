#include "nameddoublevector.h"
#include <algorithm>

namespace msh {

NamedDoubleVector::NamedDoubleVector()
{
    name_ = "undef";
}

NamedDoubleVector::NamedDoubleVector(const std::string &name, const std::vector<double> &vector)
{
    name_ = name;
    vector_ = vector;
}

NamedDoubleVector::NamedDoubleVector(const NamedDoubleVector &ndv)
{
    name_ = ndv.name_;
    vector_ = ndv.vector_;
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

NamedDoubleVector::size_type NamedDoubleVector::size() const
{
    return vector_.size();
}

double &NamedDoubleVector::operator[](const NamedDoubleVector::size_type &i)
{
    return vector_[i];
}

const double &NamedDoubleVector::operator[](const NamedDoubleVector::size_type &i) const
{
    return vector_[i];
}

double NamedDoubleVector::min() const
{
    return *std::min_element(vector_.begin(), vector_.end());
}

double NamedDoubleVector::min(NamedDoubleVector::size_type &index) const
{
    std::vector<double>::const_iterator result = std::min_element(vector_.begin(), vector_.end());
    index = std::distance(vector_.begin(), result);
    return *result;
}

double NamedDoubleVector::max() const
{
    return *std::max_element(vector_.begin(), vector_.end());
}

double NamedDoubleVector::max(NamedDoubleVector::size_type &index) const
{
    std::vector<double>::const_iterator result = std::max_element(vector_.begin(), vector_.end());
    index = std::distance(vector_.begin(), result);
    return *result;
}

}
