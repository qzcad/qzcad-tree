#ifndef HEXAHEDRALFEM_H
#define HEXAHEDRALFEM_H

#include "hexahedralmesh3d.h"
#include "femcondition3d.h"
#include "mechanicalparameters3d.h"
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
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters3D &parameters, FEMCondition3DPointer boundaryForce, std::vector<FEMCondition3DPointer> boundaryConditions);
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку шестигранных элементов
     * @param parameters Значения механических параметров (модуль Юнга и коэффициент Пуассона)
     * @param boundaryForces Массив поверхностных нагрузок
     * @param boundaryConditions Массмв граничных условий
     */
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters3D &parameters, std::vector<FEMCondition3DPointer> boundaryForces, std::vector<FEMCondition3DPointer> boundaryConditions);
    /**
     * @brief Конструктор для расчета многослойных конструкций
     * @param mesh Указатель на сетку шестигранных элементов
     * @param parameters Значения механических параметров (модуль Юнга и коэффициент Пуассона) - массив, послойное представление
     * @param boundaryForces Массив поверхностных нагрузок
     * @param boundaryConditions Массмв граничных условий
     */
    HexahedralFEM(HexahedralMesh3D* mesh, std::vector<MechanicalParameters3D> parameters, std::vector<FEMCondition3DPointer> boundaryForces, std::vector<FEMCondition3DPointer> boundaryConditions);
    /**
     * @brief Установить значения перемещения в значения узлах
     * @param mesh Сетка, в которой устанавливается значение в узлах
     * @param direction Направление перемещения (0, 1 или 2)
     */
    void setNodeDisplacement(HexahedralMesh3D* mesh, const UInteger &direction);
    /**
     * @brief Установить напряжения SigmaX в значения на элементе
     * @param mesh Указатель на сетку
     */
    void setElementSigmaX(HexahedralMesh3D * mesh);
    /**
     * @brief Установить напряжения SigmaY в значения на элементе
     * @param mesh Указатель на сетку
     */
    void setElementSigmaY(HexahedralMesh3D * mesh);
    /**
     * @brief Установить напряжения SigmaZ в значения на элементе
     * @param mesh Указатель на сетку
     */
    void setElementSigmaZ(HexahedralMesh3D * mesh);
    /**
     * @brief Установить напряжения TauXY в значения на элементе
     * @param mesh Указатель на сетку
     */
    void setElementTauXY(HexahedralMesh3D * mesh);
    /**
     * @brief Установить напряжения TauYZ в значения на элементе
     * @param mesh Указатель на сетку
     */
    void setElementTauYZ(HexahedralMesh3D * mesh);
    /**
     * @brief Установить напряжения TauZX в значения на элементе
     * @param mesh Указатель на сетку
     */
    void setElementTauZX(HexahedralMesh3D * mesh);
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
     * @param params Механические параметры материала
     * @param D Матрица упругости (матрица 6*6, перезаписуется)
     */
    void buildElasticMatrix(const MechanicalParameters3D &params, FloatingMatrix &D);
    /**
     * @brief Процедура построени матрицы упргости для многослойных конструкций
     * @param params Массив механических параметров материала
     * @param D Массив матриц упругости
     */
    void buildElasticMatrix(std::vector<MechanicalParameters3D> params, FloatingMatrix D[]);
    /**
     * @brief Процедура построения глобальной матрицы жесткости
     * @param mesh Указатель на сетку
     * @param D Матрица упругости
     * @param globalMatrix Глобальная матрица жесткости (результат)
     */
    void assebly(HexahedralMesh3D* mesh, const FloatingMatrix &D, GlobalMatrix &globalMatrix);
    /**
     * @brief Процедура построения глобальной матрицы жесткости для многослойных конструкций
     * @param mesh Указатель на сетку
     * @param D Массив матриц упругости
     * @param globalMatrix Глобальная матрица жесткости (результат)
     */
    void assebly(HexahedralMesh3D* mesh, FloatingMatrix D[], GlobalMatrix &globalMatrix);
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
    /**
     * @brief Процедура вычисления напряжений
     * @param mesh Указатель на сетку
     * @param D Матрица упругости
     */
    void recoverStress(HexahedralMesh3D* mesh, const FloatingMatrix &D);
private:
    FloatingVector displacement;
    std::vector<Floating> sigmaX;
    std::vector<Floating> sigmaY;
    std::vector<Floating> sigmaZ;
    std::vector<Floating> tauXY;
    std::vector<Floating> tauYZ;
    std::vector<Floating> tauZX;
};

#endif // HEXAHEDRALFEM_H
