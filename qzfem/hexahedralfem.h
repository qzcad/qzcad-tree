/**
  * @author Сергей Чопоров
  * @date 30/07/2014
  * @version 1.0.5
  * */
#ifndef HEXAHEDRALFEM_H
#define HEXAHEDRALFEM_H

#include "hexahedralmesh3d.h"
#include "femcondition3d.h"
#include "mechanicalparameters3d.h"
#include "floatingvector.h"
#include "floatingmatrix.h"
#include "globalmatrix.h"

using namespace msh;
/**
 * @brief Класс для определения перемещений и напраяжений при помощи МКЭ на шестигранных элементах
 */
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
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters3D &parameters, FEMCondition3DPointer boundaryForce, const std::vector<FEMCondition3DPointer> &boundaryConditions);
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку шестигранных элементов
     * @param parameters Значения механических параметров (модуль Юнга и коэффициент Пуассона)
     * @param boundaryForces Массив поверхностных нагрузок
     * @param boundaryConditions Массмв граничных условий
     */
    HexahedralFEM(HexahedralMesh3D* mesh, const MechanicalParameters3D &parameters, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions);
    /**
     * @brief Конструктор для расчета многослойных конструкций
     * @param mesh Указатель на сетку шестигранных элементов
     * @param parameters Значения механических параметров (модуль Юнга и коэффициент Пуассона) - массив, послойное представление
     * @param boundaryForces Массив поверхностных нагрузок
     * @param boundaryConditions Массмв граничных условий
     */
    HexahedralFEM(HexahedralMesh3D* mesh, const std::vector<MechanicalParameters3D> &parameters, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions);
    /**
     * @brief Получить значения перемещений в первом направлении
     * @return Значения перемещений в первом направлении
     */
    std::vector<double> u() const;
    /**
     * @brief Получить значения перемещений во втором направлении
     * @return Значения пермещений во втором направлении
     */
    std::vector<double> v() const;
    /**
     * @brief Получить значения перемещений в третьем направлении
     * @return Значения пермещений в третьем направлении
     */
    std::vector<double> w() const;
    /// Методы для получения значений напряжений на элементе
    /// @{
    std::vector<double> sigmaX() const;
    std::vector<double> sigmaY() const;
    std::vector<double> sigmaZ() const;
    std::vector<double> tauXY() const;
    std::vector<double> tauYZ() const;
    std::vector<double> tauZX() const;
    /// @}

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
    void buildElasticMatrix(const std::vector<MechanicalParameters3D> &params, FloatingMatrix D[]);
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
    void processBoundaryConditions(HexahedralMesh3D* mesh, const std::vector<FEMCondition3DPointer> &boundaryConditions, GlobalMatrix &globalMatrix, FloatingVector &force);
    /**
     * @brief Процедура для поиска и вывода на экран экстремальных значений перемещений
     * @param nodesCount Количество узлов в сетке
     */
    void printDisplacementExtremum();
    /**
     * @brief Процедура вычисления напряжений
     * @param mesh Указатель на сетку
     * @param D Матрица упругости
     */
    void recoverStress(HexahedralMesh3D* mesh, const FloatingMatrix &D);
    void recoverStress(HexahedralMesh3D* mesh, const FloatingMatrix D[]);
    /**
     * @brief Метод извлечения компонент перемещения из массива
     * @param displacement Массив пермещений
     * @param nodesCount Количество узлов в сетке
     */
    void displacementToUVW(const FloatingVector &displacement, const UInteger &nodesCount);
private:
    std::vector<double> u_; //!< Перемещения в первом направлении (x)
    std::vector<double> v_; //!< Перемещения во втором направлении (y)
    std::vector<double> w_; //!< Перемещения в третьем направлении (z)
    std::vector<double> sigmaX_; //!< Нормальные компоненты напряжения, параллельные первому направлению (x)
    std::vector<double> sigmaY_; //!< Нормальные компоненты напряжения, параллельные второму направлению (y)
    std::vector<double> sigmaZ_; //!< Нормальные компоненты напряжения, параллельные третьему направлению (z)
    std::vector<double> tauXY_; //!< Касательные компоненты напряжения, в плоскости 1-2 (x-y)
    std::vector<double> tauYZ_; //!< Касательные компоненты напряжения, в плоскости 2-3 (y-z)
    std::vector<double> tauZX_; //!< Касательные компоненты напряжения, в плоскости 3-1 (z-x)
};

#endif // HEXAHEDRALFEM_H
