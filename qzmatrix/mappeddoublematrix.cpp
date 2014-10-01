#include <iostream>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif
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

DoubleVector MappedDoubleMatrix::conjugateGradient(const DoubleVector &B, double epsilon, unsigned int niter, bool printMessages, unsigned int messageStep) const
{
    DoubleVector X(size_, 0.0); // начальное приближение - вектор нулей
    DoubleVector resid(size_); // невязка
    DoubleVector direction; // направление поиска
//    DoubleVector resid_old; // невязка на предыдущей итерации
    DoubleVector temp(size_); // ременное хранилище для обмена данными
    double resid_norm; // норма невязки
    double alpha;
    double beta;

    double resid_resid, resid_resid_new;

    residual(X, B, resid);

    direction = resid;

    resid_norm = resid.norm_2();

    if (printMessages) std::cout << "Начальная невязка: " << resid_norm << std::endl;
    if (resid_norm > epsilon)
    {
        resid_resid = resid * resid;
        for (unsigned int i = 0; i < niter; i++)
        {
            product(direction, temp);

            alpha = (resid_resid) / (direction * temp);
            X += alpha * direction;
            resid -= alpha * temp;
            resid_resid_new = resid * resid;
            resid_norm = sqrt(resid_resid_new);
            if (resid_norm <= epsilon)
            {
                if (printMessages)
                    std::cout << "Решение найдено. Итераций: " << i << ", невязка: " << resid_norm << std::endl;
                break;
            }
            if (printMessages && (i % messageStep == 0))
                std::cout << i << ", невязка: " << resid_norm << std::endl;

            beta = (resid_resid_new) / (resid_resid);
            // d = r + d*beta
            direction.scale(beta);
            direction += resid;
            //
            resid_resid = resid_resid_new;
        }
    }
    return X;
}

DoubleVector MappedDoubleMatrix::product(const DoubleVector &dv) const
{
    DoubleVector mul(size_);
#ifdef WITH_OPENMP
    omp_set_num_threads(omp_get_max_threads()); // использовать максимальное количество потоков
#pragma omp parallel for
#endif
    for (size_type i = 0; i < size_; i++)
    {
        double sum = 0.0;
        for (MappedDoubleVector::iterator it = data_[i].begin(); it != data_[i].end(); it++)
        {
            sum += it->second * dv[it->first];
        }
        mul[i] = sum;
    }
    return mul;
}

void MappedDoubleMatrix::product(const DoubleVector &dv, DoubleVector &res) const
{
#ifdef WITH_OPENMP
    omp_set_num_threads(omp_get_max_threads()); // использовать максимальное количество потоков
#pragma omp parallel for
#endif
    for (size_type i = 0; i < size_; i++)
    {
        double sum = 0.0;
        for (MappedDoubleVector::iterator it = data_[i].begin(); it != data_[i].end(); it++)
        {
            sum += it->second * dv[it->first];
        }
        res[i] = sum;
    }
}

void MappedDoubleMatrix::zeroRow(size_type i)
{
    for (MappedDoubleVector::iterator it = data_[i].begin(); it != data_[i].end(); it++)
        it->second = 0.0;
}

void MappedDoubleMatrix::zeroCol(size_type j)
{
    for (size_type i = 0; i < size_; i++)
    {
        typename MappedDoubleVector::iterator it = data_[i].find(j);
        if (it != data_[i].end())
        {
            it->second = 0.0;
        }
    }
}

void MappedDoubleMatrix::zeroSym(size_type i)
{
    for (MappedDoubleVector::iterator it = data_[i].begin(); it != data_[i].end(); it++)
    {
        it->second = 0.0;
        data_[it->first][i] = 0.0;
    }
}

DoubleVector operator *(const MappedDoubleMatrix &mdm, const DoubleVector dv)
{
    size_type size = mdm.size();
    DoubleVector mul(size);
#ifdef WITH_OPENMP
    omp_set_num_threads(omp_get_max_threads()); // использовать максимальное количество потоков
#pragma omp parallel for
#endif
    for (size_type i = 0; i < size; i++)
    {
        double sum = 0.0;
        for (MappedDoubleVector::iterator it = mdm.data_[i].begin(); it != mdm.data_[i].end(); it++)
        {
            sum += it->second * dv[it->first];
        }
        mul[i] = sum;
    }
    return mul;
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

void MappedDoubleMatrix::residual(const DoubleVector &X, const DoubleVector &B, DoubleVector &R) const
{
    for (size_type i = 0; i < size_; i++)
    {
        double t = B[i];
        for (MappedDoubleVector::iterator it = data_[i].begin(); it != data_[i].end(); it++)
        {
            t -= it->second * X[it->first];
        }
        R[i] = t;
    }
}
}

