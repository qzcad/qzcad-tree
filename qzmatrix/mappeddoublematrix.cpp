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

size_type MappedDoubleMatrix::size(size_type i) const
{
    return data_[i].size();
}

double MappedDoubleMatrix::data(size_type i, size_type j) const
{
    typename MappedDoubleVector::iterator it = data_[i].find(j);
    if (it != data_[i].end())
    {
        return it->second;
    }
    return 0;
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
//            std::cout << direction.norm_2() << "    " << temp.norm_2() << std::endl;
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
        double sum = 0;
        for (MappedDoubleVector::const_iterator it = data_[i].begin(); it != data_[i].end(); it++)
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
        double sum = 0;
        for (MappedDoubleVector::const_iterator it = data_[i].begin(); it != data_[i].end(); it++)
        {
            sum += it->second * dv[it->first];
        }
        res[i] = sum;
    }
}

void MappedDoubleMatrix::zeroRow(size_type i)
{
//    for (MappedDoubleVector::iterator it = data_[i].begin(); it != data_[i].end(); it++)
//        it->second = 0;
    data_[i].clear();
}

void MappedDoubleMatrix::zeroCol(size_type j)
{
    for (size_type i = 0; i < size_; i++)
    {
//        typename MappedDoubleVector::iterator it = data_[i].find(j);
//        if (it != data_[i].end())
//        {
//            it->second = 0;
//        }
        data_[i].erase(j);
    }
}

void MappedDoubleMatrix::zeroSym(size_type i)
{
    zeroCol(i);
    zeroRow(i);
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
        double sum = 0;
        for (MappedDoubleVector::const_iterator it = mdm.data_[i].begin(); it != mdm.data_[i].end(); it++)
        {
            sum += it->second * dv[it->first];
        }
        mul[i] = sum;
    }
    return mul;
}

void MappedDoubleMatrix::clear(size_type i)
{
    data_[i].clear();
}

MappedDoubleVector::const_iterator MappedDoubleMatrix::begin(size_type i) const
{
    return data_[i].begin();
}

MappedDoubleVector::const_iterator MappedDoubleMatrix::end(size_type i) const
{
    return data_[i].end();
}

DoubleVector MappedDoubleMatrix::cholesky(DoubleVector &B)
{
    DoubleVector L(size_ + size_ * (size_ - 1) / 2, 0.0);
    DoubleVector y(size_, 0.0);
    DoubleVector x(size_, 0.0);
//#ifdef WITH_OPENMP
//    omp_set_num_threads(omp_get_max_threads()); // использовать максимальное количество потоков
//#endif
    std::cout << size_ << "x" << size_ << " LL Init...";
    for (size_type i = 0; i < size_; i++)
    {
        L[i + i * (i + 1) / 2] = data_[i][i];
        for (MappedDoubleVector::const_iterator it = data_[i].begin(); it != data_[i].end(); it++)
        {
            size_type j = it->first;
            if (j < i)
            {
                L[j + i * (i + 1) / 2] = it->second;
            }
            else break;
        }
    }
    std::cout << "done" << std::endl;
    std::cout << "L...";
    size_type zeros = 0;
    for (register size_type i = 0; i < size_; i++)
    {
        for (register size_type j = 0; j < i; j++)
        {
            register double sij = L[j + i * (i + 1) / 2];
//#ifdef WITH_OPENMP
//#pragma omp parallel for reduction(+:sij)
//#endif
            for (register size_type k = 0; k < j; k++)
                sij -= L[k + i * (i + 1) / 2] * L[k + j * (j + 1) / 2];
            L[j + i * (i + 1) / 2] = sij / L[j + j * (j + 1) / 2];
            if (sij == 0.0) zeros++;
        }
        double sii = L[i + i * (i + 1) / 2];
//#ifdef WITH_OPENMP
//#pragma omp parallel for reduction(+:sii)
//#endif
        for (register size_type k = 0; k < i; k++)
            sii -= L[k + i * (i + 1) / 2] * L[k + i * (i + 1) / 2];
        L[i + i * (i + 1) / 2] = sqrt(sii);
//        std::cout << "sii = " << sii << std::endl;
        if (i % 500 == 0) std::cout << i << " ";
    }
    std::cout << "done: zeros = " << zeros << std::endl;
//    for (size_type i = 0; i < size_; i++)
//    {
//        for (size_type j = 0; j < size_; j++)
//        {
//            if (j <= i)
//                std::cout << L[j + i * (i + 1) / 2] << "\t";
//            else
//                std::cout << L[i + j * (j + 1) / 2] << "\t";
//        }
//        std::cout << std::endl;
//    }
    y[0] = B[0] / L[0];
    for (register size_type i = 1; i < size_; i++)
    {
        register double s = 0.0;
        for (register size_type j = 0; j < i; j++)
            s += L[j + i * (i + 1) / 2] * y[j];
        y[i] = (B[i] - s) / L[i + i * (i + 1) / 2];
    }
    x[size_ - 1] = y[size_ - 1] / L[size_ - 1 + (size_ - 1) * size_ / 2];
    for (register size_type ii = 1; ii < size_; ii++)
    {
        size_type i = size_ - 1 - ii;
        register double s = 0.0;
        for (register size_type j = i + 1; j < size_; j++)
            s += L[i + j * (j + 1) / 2] * x[j];
        x[i] = (y[i] - s) / L[i + i * (i + 1) / 2];
    }
    return x;
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

