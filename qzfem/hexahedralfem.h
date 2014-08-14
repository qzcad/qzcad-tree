#ifndef HEXAHEDRALFEM_H
#define HEXAHEDRALFEM_H

#include "hexahedralmesh3d.h"
#include "femcondition3d.h"
#include "mechanicalparameters.h"
#include "floatingvector.h"

using namespace msh;

class HexahedralFEM
{

public:
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters &parameters, FEMCondition3DPointer pressure, std::vector<FEMCondition3DPointer> boundary);
private:
    FloatingVector displacement;
};

#endif // HEXAHEDRALFEM_H
