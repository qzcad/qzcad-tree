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

double DoubleMatrix::det2x2() const
{
    return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
}

double DoubleMatrix::det3x3() const
{
    return data_[0][0] * (data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1]) -
            data_[0][1] * (data_[1][0] * data_[2][2] - data_[1][2] * data_[2][0]) +
            data_[0][2] * (data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0]);
}

DoubleMatrix DoubleMatrix::inverted2x2() const
{
    DoubleMatrix inv(2);
    double det = det2x2();
    inv.data_[0][0] = data_[1][1] / det;
    inv.data_[1][1] = data_[0][0] / det;
    inv.data_[0][1] = -data_[0][1] / det;
    inv.data_[1][0] = -data_[1][0] / det;
    return inv;
}

DoubleMatrix DoubleMatrix::inverted3x3() const
{
    double a = data_[0][0], b = data_[0][1], c = data_[0][2];
    double d = data_[1][0], e = data_[1][1], f = data_[1][2];
    double g = data_[2][0], h = data_[2][1], i = data_[2][2];
    double A = e*i - f*h,       D = -(b*i - c*h),   G = b*f - c*e;
    double B = -(d*i - f*g),    E = a*i - c*g,      H = -(a*f - c*d);
    double C = d*h - e*g,       F = -(a*h - b*g),   I = a*e - b*d;
    double det = det3x3();
    DoubleMatrix inv(3);
    inv.data_[0][0] = A / det;  inv.data_[0][1] = D / det;  inv.data_[0][2] = G / det;
    inv.data_[1][0] = B / det;  inv.data_[1][1] = E / det;  inv.data_[1][2] = H / det;
    inv.data_[2][0] = C / det;  inv.data_[2][1] = F / det;  inv.data_[2][2] = I / det;
    return inv;
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


