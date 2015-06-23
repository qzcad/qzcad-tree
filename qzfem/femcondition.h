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
/**
 * @brief Новая реализация класса для передачи граничных условий
 */
class FemCondition
{
public:
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~FemCondition() {}
    /**
     * @brief Проверить соответсвие точки граничным условиям
     * @param point казатель на точку, в которой необходимо проверить соответствие условиям
     * @return true, если точка соответствует граничному условию
     */
    virtual bool isApplied(msh::PointPointer point) = 0;
    /**
     * @brief Вычислить значение граничного условия в точке
     * @param point казатель на точку, в которой необходимо проверить соответствие условиям
     * @return true, если точка соответствует граничному условию
     */
    virtual double value(msh::PointPointer point) = 0;

    enum FemConditionType
    {
        INITIAL_VALUE = 0,
        NODAL_FORCE = 1,
        SURFACE_FORCE = 2,
        VOLUME_FORCE = 3
    };
    /**
     * @brief Метод возвращает тип граничного условия
     * @return Тип граничного условия
     */
    virtual FemConditionType type() const = 0;

    enum FemDirection
    {
        ALL = 0,
        FIRST = 1,
        SECOND = 2,
        THIRD = 3,
        FOURTH = 4,
        FIFTH = 5,
        SIXTH = 6,
        SEVENTH = 7,
        EIGHTH = 8,
        NINETH = 9,
        TENTH = 10
    };

    virtual FemDirection direction() const = 0;
};

#endif // FEMCONDITION_H
