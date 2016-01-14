#include "fem.h"

#include <iostream>
#include <math.h>
#include "rowdoublematrix.h"

Fem::Fem(Mesh *mesh, UInteger freedom_value)
{
    mesh_ = mesh;
    freedom_ = freedom_value;
    dimension_ = mesh->nodesCount() * freedom_;
    std::cout << "<== The Finite Element Analysis Engine ==>" << std::endl;
    std::cout << "Mesh: " << mesh->nodesCount() << " nodes, " << mesh->elementsCount() << " elements." << std::endl;
    std::cout << "Freedom: " << freedom_ << std::endl;
    std::cout << "Dimension: " << dimension_ << std::endl;
}

Mesh *Fem::mesh()
{
    return mesh_;
}

UInteger Fem::freedom() const
{
    return freedom_;
}

void Fem::reportNodeValues(Fem::PointFilterFunc func)
{
    UInteger nodesCount = mesh_->nodesCount();
    std::vector<UInteger> idx;
    for (UInteger i = 0; i < nodesCount; i++)
        if (func(mesh_->node(i)) == true) idx.push_back(i);

    std::cout << "x:\t";
    for (UInteger i = 0; i < idx.size(); i++)
        std::cout << mesh_->node(idx[i])->x() << "\t";

    std::cout << std::endl << "y:\t";
    for (UInteger i = 0; i < idx.size(); i++)
        std::cout << mesh_->node(idx[i])->y() << "\t";

    std::cout << std::endl << "z:\t";
    for (UInteger i = 0; i < idx.size(); i++)
        std::cout << mesh_->node(idx[i])->z() << "\t";

    for (UInteger f = 0; f < mesh_->dataVectorsCount(); f++)
    {
        std::cout << std::endl << mesh_->data(f).name() << "\t";
        std::vector<double> values = mesh_->data(f).vector();
        for (UInteger i = 0; i < idx.size(); i++)
            std::cout << values[idx[i]] << "\t";
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
        weight(2) = 128.0 / 225.0;
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
    weight.scale(0.5); // площадь единичного треугольника равна 0,5
}

void Fem::assembly(ElementPointer element, const DoubleMatrix &local, MappedDoubleMatrix &global)
{
    UInteger index_i = 0;
    UInteger index_j = 0;
    for (UInteger i = 0; i < local.rowCount(); i++)
    {
        index_i = freedom_ * element->vertexNode(i / freedom_) + (i % freedom_);

        for (UInteger j = i; j < local.colCount(); j++)
        {
            index_j = freedom_ * element->vertexNode(j / freedom_) + (j % freedom_);
            global(index_i, index_j) += local(i, j);
            if (index_i != index_j) global(index_j, index_i) = global(index_i, index_j);
        } // for j
    } // for i
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

DoubleVector Fem::solve(MappedDoubleMatrix &global, DoubleVector &force, bool cg)
{
    // решение СЛАУ
    std::cout << "Linear Equations (СЛАУ)..." << std::endl;
    if(cg)
    {
        RowDoubleMatrix rdm(global);
        return rdm.preconditionedConjugateGradient(force);
    }
    return global.cholesky(force);
}
