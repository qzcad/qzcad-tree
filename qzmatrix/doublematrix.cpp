#include <math.h>
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

const DoubleVector &DoubleMatrix::operator [](size_type i) const
{
    return data_[i];
}

reference DoubleMatrix::operator ()(size_type i, size_type j)
{
    return data_[i][j];
}

const_reference DoubleMatrix::operator ()(size_type i, size_type j) const
{
    return data_[i][j];
}

DoubleMatrix &DoubleMatrix::operator =(const DoubleMatrix &dm)
{
    if (this != &dm)
    {
        resize(dm.rowCount(), dm.colCount()); // очистка и запрос памяти
        for (size_type i = 0; i < rowCount_; i++)
        {
            data_[i] = dm.data_[i];
        }
    }
    return *this;
}

DoubleMatrix &DoubleMatrix::operator +=(const DoubleMatrix &dm)
{
    for (size_type i = 0; i < rowCount(); i++)
    {
        for (size_type j = 0; j < colCount(); j++)
        {
            data_[i][j] += dm.data_[i][j];
        }
    }
    return *this;
}

DoubleMatrix operator *(const double &d, const DoubleMatrix &dm)
{
    DoubleMatrix m(dm);
    for (size_type i = 0; i < m.rowCount(); i++)
    {
        for (size_type j = 0; j < m.colCount(); j++)
        {
            m.data_[i][j] *= d;
        }
    }
    return m;
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

DoubleMatrix DoubleMatrix::inverted() const
{
    size_type dimension = rowCount_;
    DoubleMatrix augmentedmatrix(dimension, 2 * dimension, 0.0);
    size_type i, j, k, temp;
    DoubleMatrix inv(dimension, 0.0);

    for (i = 0; i < dimension; i++)
        for (j = dimension; j < 2 * dimension; j++)
            if (i == j%dimension)
                augmentedmatrix[i][j] = 1.0;
            else
                augmentedmatrix[i][j] = 0.0;

    for (i = 0; i < dimension; i++)
        for (j = 0; j < dimension; j++)
            augmentedmatrix[i][j] = data_[i][j];

    for (j = 0; j < dimension; j++)
    {
        temp = j;
        /* finding maximum jth column element in last (dimension-j) rows */
        for (i = j + 1; i < dimension; i++)
            if (augmentedmatrix[i][j] > augmentedmatrix[temp][j])
                temp = i;

        if (fabs(augmentedmatrix[temp][j]) < 1.0E-150)
        {
            return inv;
        }
        /* swapping row which has maximum jth column element */
        if (temp != j)
            for (k = 0; k < 2 * dimension; k++)
            {
                double temporary = augmentedmatrix[j][k];
                augmentedmatrix[j][k] = augmentedmatrix[temp][k];
                augmentedmatrix[temp][k] = temporary;
            }
        /* performing row operations to form required identity matrix out
           of the input matrix */
        for (i = 0; i < dimension; i++)
            if (i != j)
            {
                double r = augmentedmatrix[i][j];
                for (k = 0; k < 2 * dimension; k++)
                    augmentedmatrix[i][k] -= augmentedmatrix[j][k] *
                            r / augmentedmatrix[j][j];
            } else
            {
                double r = augmentedmatrix[i][j];
                for (k = 0; k < 2 * dimension; k++)
                    augmentedmatrix[i][k] /= r;
            }
    }


    for (i = 0; i < dimension; i++)
        for (j = dimension; j < 2 * dimension; j++)
            inv[i][j-dimension] = augmentedmatrix[i][j];
    return inv;

}


DoubleMatrix operator +(const DoubleMatrix &a, const DoubleMatrix &b)
{
    DoubleMatrix sum(a.rowCount(), a.colCount());
    for (size_type i = 0; i < a.rowCount(); i++)
    {
        for (size_type j = 0; j < a.colCount(); j++)
            sum.data_[i][j] = a.data_[i][j] + b.data_[i][j];
    }
    return sum;
}

DoubleMatrix operator -(const DoubleMatrix &a, const DoubleMatrix &b)
{
    DoubleMatrix c(a.rowCount(), a.colCount());
    for (size_type i = 0; i < a.rowCount(); i++)
    {
        for (size_type j = 0; j < a.colCount(); j++)
            c.data_[i][j] = a.data_[i][j] - b.data_[i][j];
    }
    return c;
}

DoubleVector operator *(const DoubleMatrix &a, const DoubleVector &b)
{
    DoubleVector mul(a.rowCount());
    for (size_type i = 0; i < a.rowCount(); i++)
    {
        double sum = 0.0;
        for (size_type k = 0; k < a.colCount(); k++)
        {
            sum += a.data_[i][k] * b.data(k);
        }
        mul[i] = sum;
    }
    return mul;
}

DoubleMatrix operator *(const DoubleMatrix &a, const DoubleMatrix &b)
{
    DoubleMatrix mul(a.rowCount(), b.colCount());
    for (size_type i = 0; i < a.rowCount(); i++)
    {
        for (size_type j = 0; j < b.colCount(); j++)
        {
            double sum = 0.0;
            for (size_type k = 0; k < a.colCount(); k++)
            {
                sum += a.data_[i][k] * b.data_[k][j];
            }
            mul.data_[i][j] = sum;
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


