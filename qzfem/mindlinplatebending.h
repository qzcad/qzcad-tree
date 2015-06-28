#ifndef MINDLINPLATEBENDING_H
#define MINDLINPLATEBENDING_H

#include <functional>

#include "fem2d.h"
#include "elasticmatrix.h"
#include "femcondition.h"

#include "quadrilateralmesh2d.h"
#include "trianglemesh2d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

class MindlinPlateBending : public Fem2D
{
public:
    MindlinPlateBending(Mesh2D *mesh,
                        double thickness,
                        const ElasticMatrix &elasticMatrix,
                        std::list<FemCondition *> conditions);
};

#endif // MINDLINPLATEBENDING_H
