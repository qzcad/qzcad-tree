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
    /**
     * @brief Установить значения перемещения в значения узлах
     * @param mesh Сетка, в которой устанавливается значение в узлах
     * @param direction Направление перемещения (0, 1 или 2)
     */
    void setNodeDisplacement(HexahedralMesh3D* mesh, const UInteger &direction);
    /**
     * @brief Получить вектор перемещений для заданного узла
     * @param i Номер узла
     * @param nodesCount Колчисетво узлов в сетке
     * @return Вектор перемещений
     */
    Point3D getDisplacemementVector(const UInteger &i, const UInteger &nodesCount);
private:
    FloatingVector displacement;
};

#endif // HEXAHEDRALFEM_H
