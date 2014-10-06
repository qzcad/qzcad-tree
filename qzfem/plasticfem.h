/**
  * @author Сергей Чопоров
  * @date 06/10/2014
  * @version 1.0.0
  * */
#ifndef PLASTICFEM_H
#define PLASTICFEM_H
#include <vector>
#include "hexahedralfem.h"
/**
 * @brief Класс для расчета конечно-элементных моделей в упруго-пластичной/пластичной постановке
 */
class PlasticFem
{
public:
    /**
     * @brief Конструктор для упруго-пластичного расчета. Шастигранные конечные элементы
     * @param mesh Указатель на сетку шестигранных элементов
     * @param epsilon Массив деформаций
     * @param sigma Соответствующий массиву деформаций массив напряжений
     * @param nu Коэффициент Пуассона
     * @param boundaryForces Массив сил (шаг)
     * @param boundaryConditions Массив граничных условий
     */
    PlasticFem(HexahedralMesh3D* mesh, const std::vector<double> &epsilon, const std::vector<double> &sigma, const double &nu, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions);
};

#endif // PLASTICFEM_H
