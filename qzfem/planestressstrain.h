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
     * @param mesh Указатель на сетку треугольников или четырехугольников
     * @param thickness Толщина объекта
     * @param elasticMatrix Матрица упругости
     * @param conditions Список указателей на граничные условия и силы
     */
    PlaneStressStrain(Mesh2D *mesh,
                      double thickness,
                      const ElasticMatrix &elasticMatrix,
                      std::list<FemCondition *> conditions);
};

#endif // PLANESTRESSSTRAIN_H
