#include "namedfloatingvector.h"
#include <algorithm>

NamedFloatingVector::NamedFloatingVector()
{
    name_ = "";
}

NamedFloatingVector::NamedFloatingVector(const QString &name, const std::vector<msh::Floating> &data)
{
    name_ = name;
    data_ = data;
}

NamedFloatingVector::NamedFloatingVector(const NamedFloatingVector &nfv)
{
    name_ = nfv.name_;
    data_ = nfv.data_;
}

NamedFloatingVector::size_type NamedFloatingVector::size()
{
    return data_.size();
}

msh::Floating NamedFloatingVector::data(const NamedFloatingVector::size_type &i)
{
    return data_[i];
}

msh::Floating &NamedFloatingVector::operator[](const NamedFloatingVector::size_type &i)
{
    return data_[i];
}

msh::Floating NamedFloatingVector::min() const
{
    return *std::min_element(data_.begin(), data_.end());
}

msh::Floating NamedFloatingVector::max() const
{
    return *std::max_element(data_.begin(), data_.end());
}

QString NamedFloatingVector::name() const
{
    return name_;
}

void NamedFloatingVector::setName(const QString &name)
{
    name_ = name;
}
std::vector<msh::Floating> NamedFloatingVector::data() const
{
    return data_;
}

void NamedFloatingVector::setData(const std::vector<msh::Floating> &data)
{
    data_ = data;
}


