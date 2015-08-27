#ifndef MINDLINSHELLBENDING_H
#define MINDLINSHELLBENDING_H

#include <functional>
#include <list>

#include "fem2d.h"
#include "elasticmatrix.h"
#include "femcondition.h"

#include "quadrilateralmesh3d.h"
#include "trianglemesh3d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

class MindlinShellBending : public Fem2D
{
public:
    MindlinShellBending(Mesh3D * mesh,
                        double thickness,
                        const ElasticMatrix &elasticMatrix,
                        std::list<FemCondition *> conditions);
    MindlinShellBending(Mesh3D * mesh,
                        const std::vector<double> &thickness,
                        const std::vector<ElasticMatrix> &elasticMatrix,
                        std::list<FemCondition *> conditions);
private:
    DoubleMatrix cosinuses(const Point3D &A, const Point3D &B, const Point3D &C);
};

#endif // MINDLINSHELLBENDING_H
