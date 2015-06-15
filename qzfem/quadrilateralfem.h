/**
  * @author Сергей Чопоров
  * @date 15/06/2015
  * @version 1.0.5
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef QUADRILATERALFEM_H
#define QUADRILATERALFEM_H

#include <functional>

#include "quadrilateralmesh2d.h"
#include "doublevector.h"
#include "doublematrix.h"
#include "mappeddoublematrix.h"

using namespace msh;
using namespace mtx;
/**
 * @brief Класс для решения задач МКЭ на базе четырехугольных конечных элементов
 */
class QuadrilateralFEM
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    QuadrilateralFEM();
//    void planeStress(QuadrilateralMesh2D *mesh, std::function<double(double, double)> func)
private:
};

#endif // QUADRILATERALFEM_H
