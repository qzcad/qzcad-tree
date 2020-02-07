#include <iostream>
#include <float.h>
#include <math.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

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

void DoubleVector::resize(size_type size, const_reference initValue)
{
    resize(size);
    set(initValue);
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

double DoubleVector::norm_1() const
{
    double sum = 0.0;
    for (const_pointer p = begin(); p != end(); p++)
    {
        sum += fabs(*p);
    }
    return sum;
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

double DoubleVector::norm_inf() const
{
    double n = 0.0;
    for (const_pointer p = begin(); p != end(); p++)
    {
        double val = fabs(*p);
        if (val > n)
            n = val;
    }
    return n;
}

DoubleVector DoubleVector::dotProduct(const DoubleVector &dv) const
{
    DoubleVector result(size_);
    const_pointer pdv = dv.begin();
    pointer rdv = result.begin();
    for (const_pointer p = begin(); p != end(); p++, pdv++, rdv++)
    {
        (*rdv) = (*p) * (*pdv);
    }
    return result;
}

reference DoubleVector::operator ()(size_type i)
{
    return data_[i];
}

reference DoubleVector::operator [](size_type i)
{
    return data_[i];
}

const_reference DoubleVector::operator [](size_type i) const
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

double DoubleVector::operator *(const DoubleVector &dv) const
{
    register double sum = 0.0;
//#ifdef WITH_OPENMP
//    omp_set_num_threads(omp_get_max_threads()); // использовать максимальное количество потоков
//#pragma omp parallel for reduction(+:sum)
//#endif
//    for (size_type i = 0; i < size_; i++)
//    {
//        sum += data_[i] * dv.data_[i];
//    }
    const_pointer pdv = dv.begin();
    for (const_pointer p = begin(); p != end(); p++, pdv++)
    {
        sum += (*p) * (*pdv);
    }

    return sum;
}

void DoubleVector::scale(const double &d)
{
//    for (size_type i = 0; i < size_; i++)
//    {
//        data_[i] *= d;
//    }
    for (pointer p = begin(); p != end(); p++) (*p) *= d;
}

DoubleVector &DoubleVector::operator +=(const DoubleVector &vec)
{
//    for (size_type i = 0; i < size_; i++)
//    {
//        data_[i] += vec.data_[i];
//    }
    const_pointer pdv = vec.begin();
    for (pointer p = begin(); p != end(); p++, pdv++)
    {
        (*p) += (*pdv);
    }
    return *this;
}

DoubleVector &DoubleVector::operator -=(const DoubleVector &vec)
{
//    for (size_type i = 0; i < size_; i++)
//    {
//        data_[i] -= vec.data_[i];
//    }
    const_pointer pdv = vec.begin();
    for (pointer p = begin(); p != end(); p++, pdv++)
    {
        (*p) -= (*pdv);
    }
    return *this;
}

DoubleVector operator *(const double &d, const DoubleVector &vec)
{
    DoubleVector m(vec.size_);
    for (size_type i = 0; i < vec.size_; i++)
    {
        m.data_[i] = d * vec.data_[i];
    }
    return m;
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

std::vector<double> DoubleVector::to_std()
{
    std::vector<double> vec;
    vec.assign(data_, data_ + size_);
    return vec;
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
