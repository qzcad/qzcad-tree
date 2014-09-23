#include "doublematrix.h"

namespace mtx {
DoubleMatrix::DoubleMatrix()
{
    data_ = 0;
    rowCount_ = 0;
}

DoubleMatrix::DoubleMatrix(size_type size)
{
    alloc(size, size);
}

DoubleMatrix::DoubleMatrix(size_type size, const_reference initValue)
{
    alloc(size, size);
    set(initValue);
}

DoubleMatrix::DoubleMatrix(size_type rows, size_type cols)
{
    alloc(rows, cols);
}

DoubleMatrix::DoubleMatrix(size_type rows, size_type cols, const_reference initValue)
{
    alloc(rows, cols);
    set(initValue);
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix &dm)
{
    alloc(dm.rowCount(), dm.colCount());
    for (size_type i = 0; i < rowCount_; i++)
    {
        data_[i] = dm.data_[i];
    }
}

DoubleMatrix::~DoubleMatrix()
{
    clear();
}

size_type DoubleMatrix::rowCount() const
{
    return rowCount_;
}

size_type DoubleMatrix::colCount() const
{
    if (rowCount_ == 0)
        return 0;
    return data_[0].size();
}

void DoubleMatrix::set(const_reference initValue)
{
    for (size_type i = 0; i < rowCount_; i++)
    {
        data_[i].set(initValue);
    }
}

void DoubleMatrix::setRow(size_type i, const_reference value)
{
    data_[i].set(value);
}

void DoubleMatrix::setCol(size_type j, const_reference value)
{
    for (size_type i = 0; i < rowCount_; i++)
    {
        data_[i][j] = value;
    }
}

double DoubleMatrix::data(size_type i, size_type j) const
{
    return data_[i][j];
}

void DoubleMatrix::resize(size_type rows, size_type cols)
{
    clear();
    alloc(rows, cols);
}

DoubleVector &DoubleMatrix::operator [](size_type i)
{
    return data_[i];
}

reference DoubleMatrix::operator ()(size_type i, size_type j)
{
    return data_[i][j];
}

DoubleMatrix &DoubleMatrix::operator =(const DoubleMatrix &dm)
{
    if (this != &dm)
    {
        alloc(dm.rowCount(), dm.colCount());
        for (size_type i = 0; i < rowCount_; i++)
        {
            data_[i] = dm.data_[i];
        }
    }
    return *this;
}

void DoubleMatrix::print(char separator) const
{
    for (size_type i = 0; i < rowCount_; i++)
    {
        data_[i].print(separator);
    }
}

DoubleMatrix DoubleMatrix::transpose()
{
    DoubleMatrix trans(colCount(), rowCount());
    for (size_type i = 0; i < rowCount(); i++)
    {
        for (size_type j = 0; j < colCount(); j++)
        {
            trans.data_[j][i] = data_[i][j];
        }
    }
    return trans;
}

DoubleVector operator *(const DoubleMatrix &a, const DoubleVector &b)
{
    DoubleVector mul(a.rowCount(), 0.0);
    for (size_type i = 0; i < a.rowCount(); i++)
    {
        for (size_type k = 0; k < a.colCount(); k++)
        {
            mul[i] += a.data_[i][k] * b.data(k);
        }
    }
    return mul;
}

DoubleMatrix operator *(const DoubleMatrix &a, const DoubleMatrix &b)
{
    DoubleMatrix mul(a.rowCount(), b.colCount(), 0.0);
    for (size_type i = 0; i < a.rowCount(); i++)
    {
        for (size_type j = 0; j < b.colCount(); j++)
        {
            for (size_type k = 0; k < a.colCount(); k++)
            {
                mul.data_[i][j] += a.data_[i][k] * b.data_[k][j];
            }
        }
    }
    return mul;
}

void DoubleMatrix::alloc(size_type rows, size_type cols)
{
    rowCount_ = rows;
    if (rows > 0 && cols > 0)
    {
        data_ = new DoubleVector[rowCount_];
        for (size_type i = 0; i < rowCount_; i++)
        {
            data_[i].resize(cols);
        }
    }
}

void DoubleMatrix::clear()
{
    if (data_ != 0)
    {
        delete []data_;
        rowCount_ = 0;
        data_ = 0;
    }
}
}


