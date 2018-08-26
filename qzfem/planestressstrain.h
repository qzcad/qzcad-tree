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
     * @param alphaT Коэффициент температурного напряжения
     */
    PlaneStressStrain(Mesh2D *mesh,
                      double thickness,
                      const DoubleMatrix &elasticMatrix,
                      const std::list<FemCondition *> &conditions,
                      double alphaT = 0.0);
    virtual void solve(std::function<double(double, double)> func, double delta=0.2, int maxiter=3);
protected:
    /**
     * @brief Метод для построения глобальной матрицы системы
     */
    virtual void buildGlobalMatrix();
    /**
     * @brief Метод для построения вектора системы
     */
    virtual void buildGlobalVector();
    /**
     * @brief Метод для обработки результтов решения
     * @param nodalValues
     */
    virtual void processSolution(const DoubleVector &displacement);
    DoubleVector adaptationVector(const DoubleVector &displacement);
protected:
    double thickness_; //!< Толщина объекта
    DoubleMatrix D_; //!< Матрица упругости
    double alpha_; //!< Коэффициент температурного напряжения
};

#endif // PLANESTRESSSTRAIN_H
