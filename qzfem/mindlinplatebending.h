#ifndef MINDLINPLATEBENDING_H
#define MINDLINPLATEBENDING_H

#include <functional>

#include "fem2d.h"
#include "elasticmatrix.h"

#include "quadrilateralmesh2d.h"
#include "trianglemesh2d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

class MindlinPlateBending : public Fem2D
{
public:
    typedef std::function<double(double, double)> DistributedValueFunc;
    MindlinPlateBending(QuadrilateralMesh2D *mesh, double thickness, BoundaryConditionFunction fixFunc, DistributedValueFunc distributed);
};

#endif // MINDLINPLATEBENDING_H
