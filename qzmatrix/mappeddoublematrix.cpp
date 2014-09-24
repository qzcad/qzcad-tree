#include <iostream>
#include "mappeddoublematrix.h"

namespace mtx {
MappedDoubleMatrix::MappedDoubleMatrix()
{
    size_ = 0;
    data_ = 0;
}

MappedDoubleMatrix::MappedDoubleMatrix(size_type size)
{
    alloc(size);
}

MappedDoubleMatrix::MappedDoubleMatrix(const MappedDoubleMatrix &mdm)
{
    alloc(mdm.size_);
    for (size_type i = 0; i < size_; i++)
    {
        data_[i]= mdm.data_[i];
    }
}

MappedDoubleMatrix::~MappedDoubleMatrix()
{
    clear();
}

size_type MappedDoubleMatrix::size() const
{
    return size_;
}

double MappedDoubleMatrix::data(size_type i, size_type j) const
{
    typename MappedDoubleVector::iterator it = data_[i].find(j);
    if (it != data_[i].end())
    {
        return it->second;
    }
    return 0.0;
}

void MappedDoubleMatrix::resize(size_type size)
{
    clear();
    alloc(size);
}

MappedDoubleVector &MappedDoubleMatrix::operator [](size_type i)
{
    return data_[i];
}

reference MappedDoubleMatrix::operator ()(size_type i, size_type j)
{
    return data_[i][j];
}

MappedDoubleMatrix &MappedDoubleMatrix::operator =(const MappedDoubleMatrix &mdm)
{
    if (this != &mdm)
    {
        clear();
        alloc(mdm.size_);
        for (size_type i = 0; i < size_; i++)
        {
            data_[i]= mdm.data_[i];
        }
    }
    return *this;
}

void MappedDoubleMatrix::print(char separator) const
{
    for (size_type i = 0; i < size_; i++)
    {
        std::cout << "[" << separator;
        for (size_type j = 0; j < size_; j++)
        {
            std::cout << data(i, j) << separator;
        }
        std::cout << "]" << std::endl;
    }
}

void MappedDoubleMatrix::alloc(size_type size)
{
    size_ = size;
    data_ = new MappedDoubleVector[size_];
}

void MappedDoubleMatrix::clear()
{
    if (data_ != 0)
    {
        delete []data_;
        data_ = 0;
        size_ = 0;
    }
}
}

