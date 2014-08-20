#ifndef HEXAHEDRALFEM_H
#define HEXAHEDRALFEM_H

#include "hexahedralmesh3d.h"
#include "femcondition3d.h"
#include "mechanicalparameters.h"
#include "floatingvector.h"
#include "floatingmatrix.h"
#include "globalmatrix.h"

using namespace msh;

class HexahedralFEM
{

public:
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку шестигранных элементов
     * @param parameters Значения механических параметров (модуль Юнга и коэффициент Пуассона)
     * @param boundaryForce Поверхностная нагрузка
     * @param boundaryConditions Массмв граничных условий
     */
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters &parameters, FEMCondition3DPointer boundaryForce, std::vector<FEMCondition3DPointer> boundaryConditions);
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку шестигранных элементов
     * @param parameters Значения механических параметров (модуль Юнга и коэффициент Пуассона)
     * @param boundaryForces Массив поверхностных нагрузок
     * @param boundaryConditions Массмв граничных условий
     */
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters &parameters, std::vector<FEMCondition3DPointer> boundaryForces, std::vector<FEMCondition3DPointer> boundaryConditions);
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
protected:
    /**
     * @brief Процедура построения матрицы упргости
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @param D Матрица упругости (матрица 6*6, перезаписуется)
     */
    void buildElasticMatrix(const Floating &E, const Floating &nu, FloatingMatrix &D);
    /**
     * @brief Процедура построения глобальной матрицы жесткости
     * @param mesh Указатель на сетку
     * @param D Матрица упругости
     * @param globalMatrix Глобальная матрица жесткости (результат)
     */
    void assebly(HexahedralMesh3D* mesh, const FloatingMatrix &D, GlobalMatrix &globalMatrix);
    /**
     * @brief Процедура учета поверхностной нагрузокуи
     * @param mesh Указатель на сетку
     * @param boundaryForce Параметры поверхностной нагрузки
     * @param force Вектор сил (должен быть инициализирован нулями, результат процедуры)
     */
    void processForce(HexahedralMesh3D* mesh, FEMCondition3DPointer boundaryForce, FloatingVector &force);
    /**
     * @brief Процедура учета граничных условий
     * @param mesh Укзатель на сетку
     * @param boundaryConditions Массив граничных условий
     * @param globalMatrix Глобальная матрица жесткости (модифицируется)
     * @param force Вектор сил (модифицируется)
     */
    void processBoundaryConditions(HexahedralMesh3D* mesh, std::vector<FEMCondition3DPointer> boundaryConditions, GlobalMatrix &globalMatrix, FloatingVector &force);
    /**
     * @brief Процедура для поиска и вывода на экран экстремальных значений перемещений
     * @param nodesCount Количество узлов в сетке
     */
    void printDisplacementExtremum(const UInteger &nodesCount);
private:
    FloatingVector displacement;
};

#endif // HEXAHEDRALFEM_H
