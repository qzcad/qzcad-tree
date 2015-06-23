/**
  * @author Сергей Чопоров
  * @date 17/06/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef PLANESTRESSSTRAIN_H
#define PLANESTRESSSTRAIN_H

#include <functional>
#include <list>

#include "fem2d.h"
#include "elasticmatrix.h"

#include "femcondition.h"

#include "quadrilateralmesh2d.h"
#include "trianglemesh2d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

/**
 * @brief Класс для решения задач плоского напряжения и плоской деформации МКЭ
 */
class PlaneStressStrain : public Fem2D
{
public:
    /**
     * @brief Конструктор для решения задачи исследования плоского напряженного или плоского деформированного состояния двумерного объекта
     * @param mesh Указатель на сетку
     * @param thickness Толщина объекта
     * @param elasticMatrix Матрица упругости
     * @param fixFunc Функция условий закрепеления (параметры - координаты; возвращает 0, если начальное значение применяется для обеих компонент перемещения, 1 - начальные значения при меняются только в первом направлении, 2 - зтолько во втором напралении)
     * @param boundaryValue Функция начальных значений (параметры - координаты; результат начальные перемещения, соответствующие узлу с указанными координатами)
     * @param nodalForce Функция узловых нагрузок
     * @param surfaceForce Функция поверхностных нагрузок
     * @param volumeForce Фукнция объемных нагрузок
     */
    PlaneStressStrain(QuadrilateralMesh2D *mesh,
                      double thickness,
                      const ElasticMatrix &elasticMatrix,
                      BoundaryConditionFunction fixFunc,
                      VectorFunction2D boundaryValue,
                      VectorFunction2D nodalForce,
                      VectorFunction2D surfaceForce,
                      VectorFunction2D volumeForce);
    /**
     * @brief Конструктор для решения задачи исследования плоского напряженного или плоского деформированного состояния двумерного объекта
     * @param mesh Указатель на сетку четырехугольников
     * @param thickness Толщина объекта
     * @param elasticMatrix Матрица упругости
     * @param conditions Список указателей на граничные условия и силы
     */
    PlaneStressStrain(QuadrilateralMesh2D *mesh,
                      double thickness,
                      const ElasticMatrix &elasticMatrix,
                      std::list<FemCondition *> conditions);
    /**
     * @brief Конструктор для решения задачи исследования плоского напряженного или плоского деформированного состояния двумерного объекта
     * @param mesh Указатель на сетку
     * @param thickness Толщина объекта
     * @param elasticMatrix Матрица упругости
     * @param fixFunc Функция условий закрепеления (параметры - координаты; возвращает 0, если начальное значение применяется для обеих компонент перемещения, 1 - начальные значения при меняются только в первом направлении, 2 - зтолько во втором напралении)
     * @param boundaryValue Функция начальных значений (параметры - координаты; результат начальные перемещения, соответствующие узлу с указанными координатами)
     * @param nodalForce Функция узловых нагрузок
     * @param surfaceForce Функция поверхностных нагрузок
     * @param volumeForce Фукнция объемных нагрузок
     */
    PlaneStressStrain(TriangleMesh2D *mesh,
                      double thickness,
                      const ElasticMatrix &elasticMatrix,
                      BoundaryConditionFunction fixFunc,
                      VectorFunction2D boundaryValue,
                      VectorFunction2D nodalForce,
                      VectorFunction2D surfaceForce,
                      VectorFunction2D volumeForce);
    /**
     * @brief Конструктор для решения задачи исследования плоского напряженного или плоского деформированного состояния двумерного объекта
     * @param mesh Указатель на сетку треугольников
     * @param thickness Толщина объекта
     * @param elasticMatrix Матрица упругости
     * @param conditions Список указателей на граничные условия и силы
     */
    PlaneStressStrain(TriangleMesh2D *mesh,
                      double thickness,
                      const ElasticMatrix &elasticMatrix,
                      std::list<FemCondition *> conditions);
    /**
     * @brief Метод возвращает название вектора узловых значений
     * @param num Номер вектора узловых значений
     * @return Название вектора узловых значений
     */
    virtual std::string nodeVectorName(UInteger num) const;
    /**
     * @brief Метод возвращает название вектора значений, определенных на элементе
     * @param num Номер вектора значений, определеных на элементе
     * @return Название вектора значений, определенных на элеменете
     */
    virtual std::string elementVectorName(UInteger num) const;
};

#endif // PLANESTRESSSTRAIN_H
