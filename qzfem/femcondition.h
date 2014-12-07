/**
  * @author Сергей Чопоров
  * @date 14/08/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef FEMCONDITION_H
#define FEMCONDITION_H
#include "pointpointer.h"
/**
 * @brief Класс граничных условий (виртуальный класс)
 */
class FEMCondition
{
public:
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~FEMCondition() {}
    /**
     * @brief Проверить соответсвие точки граничным условиям
     * @param point казатель на точку, в которой необходимо проверить соответствие условиям
     * @return true, если точка соответствует граничному условию
     */
    virtual bool isApplied(msh::PointPointer point) = 0;
};
/**
 * @brief Указатель на граничные условия
 */
typedef FEMCondition* FEMConditionPointer;
#endif // FEMCONDITION_H
