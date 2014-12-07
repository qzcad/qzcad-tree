/**
  * @author Сергей Чопоров
  * @date 14/08/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef FEMCONDITION3D_H
#define FEMCONDITION3D_H
#include "femcondition2d.h"
/**
 * @brief Граничные условия для трехмерных задач (виртуальный класс)
 */
class FEMCondition3D : public FEMCondition2D
{
public:
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~FEMCondition3D() {}
    /**
     * @brief Значение граничного условия в третьем направлении
     * @return Значение граничного условия во втором направлении
     */
    virtual double w() = 0;
    /**
     * @brief Проверить действие условия в третьем направлении
     * @return true, если действует во втором направлении
     */
    virtual bool isW() = 0;
};
/**
 * @brief Указатель на трехмерные граничные условия
 */
typedef FEMCondition3D* FEMCondition3DPointer;
#endif // FEMCONDITION3D_H
