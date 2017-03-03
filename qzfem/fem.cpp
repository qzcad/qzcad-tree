#include "fem.h"

#include <iostream>
#include <math.h>
#include "rowdoublematrix.h"
#include "consoleprogress.h"

Fem::Fem(Mesh *mesh, UInteger freedom, const std::list<FemCondition *> &conditions) :
    mesh_(mesh),
    freedom_ (freedom),
    dimension_ (mesh->nodesCount() * freedom),
    global_(mesh->nodesCount() * freedom),
    force_(mesh->nodesCount() * freedom, 0.0),
    conditions_(conditions)
{
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

void Fem::solve()
{
    mesh_->clearDataVectors();
    buildGlobalMatrix();
    buildGlobalVector();
    processInitialValues();
    DoubleVector solution = solveLinearSystem();
    processSolution(solution);
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

void Fem::assembly(ElementPointer element, const DoubleMatrix &local)
{
    UInteger index_i = 0;
    UInteger index_j = 0;
    for (UInteger i = 0; i < local.rowCount(); i++)
    {
        index_i = freedom_ * element->vertexNode(i / freedom_) + (i % freedom_);

        for (UInteger j = i; j < local.colCount(); j++)
        {
            index_j = freedom_ * element->vertexNode(j / freedom_) + (j % freedom_);
            global_(index_i, index_j) += local(i, j);
            if (index_i != index_j) global_(index_j, index_i) = global_(index_i, index_j);
        } // for j
    } // for i
}

void Fem::setInitialNodalValue(const UInteger &rowNumber, const double &value)
{
    for (UInteger j = 0; j < global_.size(); j++)
    {
        if (rowNumber != j && global_.data(rowNumber, j) != 0)
        { // см. Зенкевич, стр. 485
            force_(j) = force_(j) - global_.data(rowNumber, j) * value;
        }
    } // for j
    global_.zeroSym(rowNumber);
    force_(rowNumber) = value;
    global_(rowNumber, rowNumber) = 1.0;
}

void Fem::processInitialValues()
{
    //учет учет граничных условий
    std::cout << "Boundary Conditions...";
    for (std::list<FemCondition *>::iterator condition = conditions_.begin(); condition != conditions_.end(); condition++)
    {
        if ((*condition)->type() == FemCondition::INITIAL_VALUE)
        {
            ConsoleProgress progressBar(mesh_->nodesCount());
            for (UInteger i = 0; i < mesh_->nodesCount(); i++)
            {
                PointPointer point = mesh_->node(i);
                if ((*condition)->isApplied(point))
                {
                    FemCondition::FemDirection dir = (*condition)->direction();

                    if (dir == FemCondition::ALL || dir == FemCondition::FIRST)
                        setInitialNodalValue(freedom_ * i, (*condition)->value(point));

                    if ((dir == FemCondition::ALL || dir == FemCondition::SECOND) && freedom_ >= 2)
                        setInitialNodalValue(freedom_ * i + 1UL, (*condition)->value(point));

                    if ((dir == FemCondition::ALL || dir == FemCondition::THIRD) && freedom_ >= 3)
                        setInitialNodalValue(freedom_ * i + 2UL, (*condition)->value(point));

                    if ((dir == FemCondition::ALL || dir == FemCondition::FOURTH) && freedom_ >= 4)
                        setInitialNodalValue(freedom_ * i + 3UL, (*condition)->value(point));

                    if ((dir == FemCondition::ALL || dir == FemCondition::FIFTH) && freedom_ >= 5)
                        setInitialNodalValue(freedom_ * i + 4UL, (*condition)->value(point));

                    if ((dir == FemCondition::ALL || dir == FemCondition::SIXTH) && freedom_ >= 6)
                        setInitialNodalValue(freedom_ * i + 5UL, (*condition)->value(point));
                }
                ++progressBar;
            } // for i
        } // if
    } // iterator
}

DoubleVector Fem::solveLinearSystem(bool cg)
{
    // решение СЛАУ
    std::cout << "Linear Equations (СЛАУ)..." << std::endl;
    if(cg)
    {
        RowDoubleMatrix rdm(global_);
        return rdm.preconditionedConjugateGradient(force_);
    }
    return global_.cholesky(force_);
}
