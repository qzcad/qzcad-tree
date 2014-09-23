#include <iostream>
#include <float.h>
#include <math.h>
#include "doublevector.h"

namespace mtx {
DoubleVector::DoubleVector()
{
    size_ = 0;
    data_ = 0;
}

DoubleVector::DoubleVector(size_type size)
{
    size_ = size;
    alloc();
}

DoubleVector::DoubleVector(size_type size, const_reference initValue)
{
    size_ = size;
    alloc(initValue);
}

DoubleVector::DoubleVector(const DoubleVector &dv)
{
    size_ = dv.size_;
    alloc();
    for (size_type i = 0; i < size_; i++) data_[i] = dv.data_[i];
}

DoubleVector::~DoubleVector()
{
    clear();
}

pointer DoubleVector::begin()
{
    return data_;
}

const_pointer DoubleVector::begin() const
{
    const_pointer cptr = data_;
    return cptr;
}

pointer DoubleVector::end()
{
    return data_ + size_;
}

const_pointer DoubleVector::end() const
{
    const_pointer cptr = data_ + size_;
    return cptr;
}

size_type DoubleVector::size() const
{
    return size_;
}

double DoubleVector::data(size_type i) const
{
    return data_[i];
}

void DoubleVector::set(const_reference initValue)
{
    for (pointer p = begin(); p != end(); p++) *p = initValue;
}

void DoubleVector::resize(size_type size)
{
    clear();
    size_ = size;
    alloc();
}

double DoubleVector::min() const
{
    double minValue = DBL_MAX;
    for (const_pointer p = begin(); p != end(); p++)
    {
        if (*p < minValue) minValue = *p;
    }
    return minValue;
}

double DoubleVector::max() const
{
    double maxValue = -DBL_MAX;
    for (const_pointer p = begin(); p != end(); p++)
    {
        if (*p > maxValue) maxValue = *p;
    }
    return maxValue;
}

double DoubleVector::norm_2() const
{
    double sum = 0.0;
    for (const_pointer p = begin(); p != end(); p++)
    {
        double val = *p;
        sum += val * val;
    }
    return sqrt(sum);
}

reference DoubleVector::operator [](size_type i)
{
    return data_[i];
}

DoubleVector &DoubleVector::operator=(const DoubleVector &dv)
{
    if (this != &dv)
    {
        resize(dv.size_);
        for (size_type i = 0; i < size_; i++) data_[i] = dv.data_[i];
    }
    return *this;
}

void DoubleVector::print(char separator) const
{
    std::cout << "[" << separator;
    for (const_pointer p = begin(); p != end(); p++)
    {
        std::cout << *p << separator;
    }
    std::cout << "]" << std::endl;
}

void DoubleVector::alloc()
{
    if (size_ > 0)
    {
        data_ = new double[size_];
    }
}

void DoubleVector::alloc(const_reference initValue)
{
    alloc();
    set(initValue);
}

void DoubleVector::clear()
{
    if (data_ != 0)
    {
        delete []data_;
        size_ = 0;
        data_ = 0;
    }
}

}
