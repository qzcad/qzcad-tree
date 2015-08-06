#ifndef MINDLINSHELLBENDING_H
#define MINDLINSHELLBENDING_H

#include <functional>

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
    MindlinShellBending(Mesh3D * mesh);
};

#endif // MINDLINSHELLBENDING_H
