#include "fem.h"

#include <iostream>
#include <math.h>

Fem::Fem(Mesh *mesh)
{
    mesh_ = mesh;
}

Mesh *Fem::mesh()
{
    return mesh_;
}

UInteger Fem::freedom() const
{
    return freedom_;
}

std::vector<double> Fem::nodeVector(UInteger num) const
{
    UInteger nodesCount = mesh_->nodesCount();
    std::vector<double> v(nodesCount);
    UInteger d = nodesCount * num;
    for (UInteger i = 0; i < nodesCount; i++)
    {
        v[i] = nodeValues_[i + d];
    }
    return v;
}

UInteger Fem::elementVectorsCount() const
{
    return elementVectorsCount_;
}

std::vector<double> Fem::elementVector(UInteger num) const
{
    UInteger elementsCount = mesh_->elementsCount();
    std::vector<double> v(elementsCount);
    UInteger d = elementsCount * num;
    for (UInteger i = 0; i < elementsCount; i++)
    {
        v[i] = elementValues_[i + d];
    }
    return v;
}

void Fem::printNodeValuesExtremums() const
{
    UInteger nodesCount = mesh_->nodesCount();
    for (UInteger f = 0; f < freedom_; f++)
    {
        UInteger d = nodesCount * f;
        double max = nodeValues_[0 + d];
        double min = nodeValues_[0 + d];
        for (UInteger i = 1; i < nodesCount; i++)
        {
            double v = nodeValues_[i + d];
            if (max < v) max = v;
            if (min > v) min = v;
        }
        std::cout << f << "\t:\t" << min << " <= " << nodeVectorName(f) << " <= " << max << std::endl;
    }
}

void Fem::printElementValuesExtremums() const
{
    UInteger elementsCount = mesh_->elementsCount();
    for (UInteger f = 0; f < elementVectorsCount_; f++)
    {
        UInteger d = elementsCount * f;
        double max = elementValues_[0 + d];
        double min = elementValues_[0 + d];
        for (UInteger i = 1; i < elementsCount; i++)
        {
            double v = elementValues_[i + d];
            if (max < v) max = v;
            if (min > v) min = v;
        }
        std::cout << f << "\t:\t" << min << " <= " << elementVectorName(f) << " <= " << max << std::endl;
    }
}

void Fem::quadrature(int count, DoubleVector &point, DoubleVector &weight)
{
    point.resize(count);
    weight.resize(count);
    switch (count)
    {
    case 1:
        point(0) = 0.0;
        weight(0) = 2.0;
        break;
    case 2:
        point(0) = -1.0 / sqrt(3.0);
        point(1) = -point(0);
        weight(0) = 1.0;
        weight(1) = weight(0);
        break;
    case 3:
        point(0) = -sqrt(3.0 / 5.0);
        point(1) = 0.0;
        point(2) = -point(0);
        weight(0) = 5.0 / 9.0;
        weight(1) = 8.0 / 9.0;
        weight(2) = weight(0);
        break;
    case 4:
        point(0) = -sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
        point(1) = -sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
        point(2) = -point(1);
        point(3) = -point(0);
        weight(0) = (18.0 - sqrt(30.0)) / 36.0;
        weight(1) = (18.0 + sqrt(30.0)) / 36.0;
        weight(2) = weight(1);
        weight(3) = weight(0);
        break;
    default:
        point(0) = -(1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
        point(1) = -(1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
        point(2) = 0.0;
        point(3) = -point(1);
        point(4) = -point(0);
        weight(0) = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
        weight(1) = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
        weight(2) = 128.0 / 255.0;
        weight(3) = weight(1);
        weight(4) = weight(0);
    }
}

void Fem::quadrature(int count, DoubleVector &xi, DoubleVector &eta, DoubleVector &weight)
{
    xi.resize(count);
    eta.resize(count);
    weight.resize(count);
    switch (count)
    {
    case 1:
        xi(0) = 1.0 / 3.0;  eta(0) = 1.0 / 3.0; weight(0) = 1.0;
        break;
    case 3://////////////////////////////////////////////////////////////
        xi(0) = 1.0 / 6.0;  eta(0) = 1.0 / 6.0; weight(0) = 1.0 / 3.0;
        //
        xi(1) = 2.0 / 3.0;  eta(1) = 1.0 / 6.0; weight(1) = 1.0 / 3.0;
        //
        xi(2) = 1.0 / 6.0;  eta(2) = 2.0 / 3.0; weight(2) = 1.0 / 3.0;
        break;
    case 4: /////////////////////////////////////////////////////////////
        xi(0) = 1.0 / 3.0;  eta(0) = 1.0 / 3.0; weight(0) = -27.0 / 48.0;
        //
        xi(1) = 1.0 / 5.0;  eta(1) = 3.0 / 5.0; weight(1) = 25.0 / 48.0;
        //
        xi(2) = 1.0 / 5.0;  eta(2) = 1.0 / 5.0; weight(2) = 25.0 / 48.0;
        //
        xi(3) = 3.0 / 5.0;  eta(3) = 1.0 / 5.0; weight(3) = 25.0 / 48.0;
        break;
    default:
        xi(0) = 1.0 / 3.0;  eta(0) = 1.0 / 3.0; weight(0) = 1.0;
    }
}

void Fem::setInitialNodalValue(MappedDoubleMatrix &global, DoubleVector &force, const UInteger &rowNumber, const double value)
{
    for (UInteger j = 0; j < global.size(); j++)
    {
        if (rowNumber != j && global.data(rowNumber, j) != 0)
        { // см. Зенкевич, стр. 485
            force(j) = force(j) - global.data(rowNumber, j) * value;
        }
    } // for j
    global.zeroSym(rowNumber);
    force(rowNumber) = value;
    global(rowNumber, rowNumber) = 1.0;
}
