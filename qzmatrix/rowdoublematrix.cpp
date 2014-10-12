#include <omp.h>
#include <math.h>
#include <iostream>
#include "rowdoublematrix.h"

namespace mtx {

RowDoubleMatrix::RowDoubleMatrix()
{
    size_ = 0;
    data_ = 0;
    row_size_ = 0;
}

RowDoubleMatrix::RowDoubleMatrix(MappedDoubleMatrix &M)
{
    size_ = 0;
    data_ = 0;
    row_size_ = 0;
    assign(M);
}

RowDoubleMatrix::~RowDoubleMatrix()
{
    clear();
}

void RowDoubleMatrix::assign(MappedDoubleMatrix &M)
{
    clear();
    size_ = M.size();
    row_size_ = new size_type[size_];
    data_ = new ColumnValuePair* [size_];
    for (size_type i = 0; i < size_; i++)
    {
        size_type rs = M.size(i);
        size_type j = 0;
        row_size_[i] = rs;
        data_[i] = new ColumnValuePair[rs];

        for (MappedDoubleVector::const_iterator it = M.begin(i); it != M.end(i); it++)
        {
            data_[i][j].column = it->first;
            data_[i][j].value = it->second;
            j++;
        }
        M.clear(i);
    }
}

void RowDoubleMatrix::product(const DoubleVector &dv, DoubleVector &res) const
{
#ifdef WITH_OPENMP
    omp_set_num_threads(omp_get_max_threads()); // использовать максимальное количество потоков
#pragma omp parallel for
#endif
    for (size_type i = 0; i < size_; i++)
    {
        double sum = 0;
        for (size_type j = 0; j < row_size_[i]; j++)
        {
            ColumnValuePair pair = data_[i][j];
            sum += pair.value * dv[pair.column];
        }
        res[i] = sum;
    }
}

DoubleVector RowDoubleMatrix::conjugateGradient(const DoubleVector &B, double epsilon, unsigned int niter, bool printMessages, unsigned int messageStep) const
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

void RowDoubleMatrix::clear()
{
    if (size_ > 0)
    {
        for (size_type i = 0; i < size_; i++)
        {
            if (row_size_[i] > 0)
                delete []data_[i];
        }
        delete []data_;
        delete []row_size_;
    }
}

void RowDoubleMatrix::residual(const DoubleVector &X, const DoubleVector &B, DoubleVector &R) const
{
    for (size_type i = 0; i < size_; i++)
    {
        double t = B[i];
        for (size_type j = 0; j < row_size_[i]; j++)
        {
            ColumnValuePair pair = data_[i][j];
            t -= pair.value * X[pair.column];
        }
        R[i] = t;
    }
}

}

